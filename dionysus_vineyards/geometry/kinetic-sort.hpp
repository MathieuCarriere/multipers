#include "utilities/log.h"
#include "utilities/counter.h"

#ifdef LOGGING
static rlog::RLogChannel* rlKineticSort =           DEF_CHANNEL("geometry/kinetic-sort", rlog::Log_Debug);
static rlog::RLogChannel* rlKineticSortAudit =      DEF_CHANNEL("geometry/kinetic-sort/audit", rlog::Log_Debug);
static rlog::RLogChannel* rlKineticSortSchedule =   DEF_CHANNEL("geometry/kinetic-sort/schedule", rlog::Log_Debug);
static rlog::RLogChannel* rlKineticSortProcess =    DEF_CHANNEL("geometry/kinetic-sort/process", rlog::Log_Debug);
#endif // LOGGING

#ifdef COUNTERS
static Counter*  cKineticSort =                     GetCounter("kinetic-sort");
static Counter*  cKineticSortSwap =                 GetCounter("kinetic-sort/swap");
#endif // COUNTERS


template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
KineticSort()
{}	

template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
KineticSort(ElementIterator b, ElementIterator e, Swap swap, Simulator* simulator, const TrajectoryExtractor& te):
    te_(te)
{
	initialize(b, e, swap, simulator);
}

template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
void
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
initialize(ElementIterator b, ElementIterator e, Swap swap, Simulator* simulator)
{
	swap_ = swap;
	for (ElementIterator cur = b; cur != e; ++cur)
		list_.push_back(Node(cur, simulator->null_key()));
	schedule_swaps(list_.begin(), list_.end(), simulator);
}

/// Adds elements in the range [f,l) of the underlying data structure to the management by the KineticSort.
/// pos must be a valid iterator, i.e., it cannot be end().
template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
void
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
insert(iterator pos, ElementIterator f, ElementIterator l, Simulator* simulator)
{
	iterator previous = pos; --previous;
	if (previous != list_.end()) simulator->remove(previous->swap_event_key);

	ElementIterator cur = boost::next(previous)->element;
	while(cur != pos->element)
		list_.insert(pos->element, Node(cur++, simulator->null_key()));

	if (previous != list_.end()) 
		schedule_swaps(previous, pos, simulator);
	else
		schedule_swaps(list_.begin(), pos, simulator);
}

/// Removes pos from the KineticSort structure. Assumes that the element of the underlying data structure 
/// that *pos refers to has already been removed.
template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
void
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
erase(iterator pos, Simulator* simulator)
{
	simulator->remove(pos->swap_event_key);
	iterator prev = pos; --prev;
	list_.erase(pos);
	schedule_swaps(prev, simulator);
}

template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
void
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
update_trajectory(iterator pos, Simulator* simulator)
{
	iterator prev = boost::prior(pos);
	if (prev != list_.end())
	{
		simulator->remove(prev->swap_event_key);
		schedule_swaps(prev, simulator);
	}

	if (boost::next(pos) != list_.end())
	{
		simulator->remove(pos->swap_event_key);
		schedule_swaps(pos, simulator);
	}
}


template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
void						
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
swap(iterator pos, Simulator* simulator)
{
	// TODO: AssertMsg(boost::next(pos) != list_.end(), "Cannot swap the last element");

    Count(cKineticSortSwap);
	swap_(pos->element, simulator);

	// Remove events
	iterator prev = boost::prior(pos);
	if (prev != list_.end())
		simulator->remove(prev->swap_event_key);
	iterator next = boost::next(pos);
	simulator->remove(next->swap_event_key);

	// Swap
	list_.splice(pos, list_, next);
	
	// update events
	next->swap_event_key = pos->swap_event_key;
	static_cast<SwapEvent*>(*(next->swap_event_key))->set_position(next);
	schedule_swaps(prev, simulator);
	schedule_swaps(pos, simulator);
	//audit(simulator);
}

template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
bool
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
audit(Simulator* simulator) const
{
	typedef 		typename Simulator::Function		        Function;
	typedef 		typename Simulator::Time					Time;
	
	Time t = simulator->audit_time();
	rLog(rlKineticSortAudit, "Auditing at %s", tostring(t).c_str());

	typename NodeList::const_iterator next = list_.begin();
	typename NodeList::const_iterator cur = next++;
	Function cur_trajectory = te_(cur->element);
	while (next != list_.end())
	{
		rLog(rlKineticSortAudit, "  %s", intostring(**(cur->swap_event_key)).c_str());

		Function next_trajectory = te_(next->element);
		rLog(rlKineticSortAudit, "  Auditing:   %s, %s", tostring(cur_trajectory).c_str(),
                                                         tostring(next_trajectory).c_str());
		rLog(rlKineticSortAudit, "  Difference: %s", tostring(next_trajectory - cur_trajectory).c_str());
		rLog(rlKineticSortAudit, "  Sign at:    %s, %s", tostring(t).c_str(),
                                                         tostring(FunctionKernel::sign_at(next_trajectory - cur_trajectory, t)).c_str());
		if (FunctionKernel::sign_at(next_trajectory - cur_trajectory, t) == -1)
		{
			// rError("Audit failed at %s, %s", tostring(*cur->element).c_str(), 
            //                                  tostring(*next->element).c_str());
			return false;
		}

		cur_trajectory = next_trajectory;
		cur = next++;
	}
	if (cur != list_.end()) rLog(rlKineticSortAudit, "  %s", intostring(**(cur->swap_event_key)).c_str());
	return true;
}

		
template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
void						
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
schedule_swaps(iterator b, iterator e, Simulator* simulator)
{
	typedef 		typename Simulator::Function		        Function;
	
	iterator next = b; 
	iterator cur = next++;
	Function cur_trajectory = te_(cur->element);
	while (next != e)
	{
		Function next_trajectory = te_(next->element);
		rLog(rlKineticSortSchedule, "Next trajectory: %s", tostring(next_trajectory).c_str());
		// TODO: add assertion that (next_trajectory - cur_trajectory)(s->curren_time()) > 0
		cur->swap_event_key = simulator->add(next_trajectory - cur_trajectory, SwapEvent(this, cur));
		cur = next++;
		cur_trajectory = next_trajectory;
	}
	if (cur != e) schedule_swaps(cur, simulator);
}

template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
void						
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::
schedule_swaps(iterator i, Simulator* simulator)
{
	typedef 		typename Simulator::Function		        Function;
	
	if (i == list_.end()) return;
	if (boost::next(i) == list_.end())
	{
		i->swap_event_key = simulator->add(SwapEvent(this, i));
		return;
	}

	iterator next = boost::next(i); 
	Function i_trajectory = te_(i->element);
	Function next_trajectory = te_(next->element);
	
	//std::cout << "Updating swaps for: " << i_trajectory << ", " << next_trajectory << std::endl;
	//std::cout << "Difference:         " << next_trajectory - i_trajectory << std::endl;

	i->swap_event_key = simulator->add(next_trajectory - i_trajectory, SwapEvent(this, i));
}

/* SwapEvent */
template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
class KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::SwapEvent: public Simulator::Event
{
	public:
		typedef						typename Simulator::Event					Parent;

									SwapEvent(const SwapEvent& e):
										sort_(e.sort_), pos_(e.pos_)			{}
									SwapEvent(KineticSort* sort, iterator pos):
										sort_(sort), pos_(pos)					{}

		virtual bool				process(Simulator* s) const;
		void						set_position(iterator i)					{ pos_ = i; }
		iterator					position() const							{ return pos_; }
		std::ostream&				operator<<(std::ostream& out) const;

	private:
		KineticSort*				sort_;
		iterator					pos_;
};

template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
bool
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::SwapEvent::
process(Simulator* s) const
{ 
	rLog(rlKineticSortProcess, "Swapping. Current time: %s", tostring(s->current_time()).c_str());
	sort_->swap(pos_, s); 
	return true; 
}

template<class ElementIterator_, class TrajectoryExtractor_, class Simulator_, class Swap_>
std::ostream&				
KineticSort<ElementIterator_, TrajectoryExtractor_, Simulator_, Swap_>::SwapEvent::
operator<<(std::ostream& out) const
{
	Parent::operator<<(out) << "SwapEvent at " << sort_->trajectory_extractor()(position()->element);
	return out;
}
