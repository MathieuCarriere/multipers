#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#include <utilities/eventqueue.h>
#include <utilities/indirect.h>

#include <limits>

/**
 * Simulator class. Keeps a queue of events. Infinity is reached if the Event
 * at the front of the queue has an empty root stack. Keeps track of current time,
 * Event addition, and processes events one by one. Degeneracies are handled by
 * assuming that the FunctionKernel::Function responsible for the event must be
 * positive before the Event occurs.
 *
 * \ingroup kinetic
 */
template<class FuncKernel_, template<class Event> class EventComparison_ = std::less>
class Simulator
{
    public:
        typedef                     FuncKernel_                                 FunctionKernel;
        typedef                     typename FunctionKernel::Function           Function;
        typedef                     typename FunctionKernel::RootStack          RootStack;
        typedef                     typename FunctionKernel::RootType           RootType;
        typedef                     RootType                                    Time;

        class Event;
        typedef                     EventComparison_<Event>                     EventComparison;

        class IndirectEventComparison;
        typedef                     EventQueue<Event*, IndirectEventComparison> EventQueueS;
        typedef                     typename EventQueueS::iterator              Key;
        typedef                     typename EventQueueS::const_iterator        const_Key;


                                    Simulator(Time start = FunctionKernel::root(0)):
                                        current_(start)                         {}
                                    ~Simulator()                                { for (Key cur = queue_.top(); cur != queue_.end(); ++cur) delete *cur; }


        template<class Event_>
        Key                         add(const Event_& e);
        template<class Event_>
        Key                         add(const Function& f, const Event_& e);
        void                        process();
        //void                        update(Key k, const Function& f);

        void                        remove(Key k)                               { Event* e = *k; queue_.remove(k); delete e; }
        Key                         null_key()                                  { return queue_.end(); }

        Time                        current_time() const                        { return current_; }
        Time                        audit_time() const;
        bool                        reached_infinity() const                    { return queue_.empty() || (*queue_.top())->root_stack().empty(); }
        Time                        next_event_time() const                     { return queue_.empty() ? std::numeric_limits<Time>::infinity():(*queue_.top())->root_stack().top(); }

        Event*                      top() const                                 { return *(queue_.top()); }
        unsigned                    size() const                                { return queue_.size(); }
        unsigned                    event_count() const                         { return count_; }

        std::ostream&               operator<<(std::ostream& out) const;

    private:
        Time                        current_;
        EventQueueS                 queue_;
        unsigned                    count_;
};


/**
 * Base class for events. Stores a root stack, subclasses need to define process().
 * Event with an empty root stack compares greater than any other Event,
 * pushing those events to the end of the queue.
 */
template<class FuncKernel_, template<class Event> class EventComparison_>
class Simulator<FuncKernel_, EventComparison_>::Event
{
    public:
        typedef                     FuncKernel_                                 FunctionKernel;
        typedef                     typename FunctionKernel::RootStack          RootStack;

        virtual                     ~Event()                                    {}

        /// process() is called when the event is at the top of the queue
        /// in the simulator.
        /// Returns true if the event needs to remain in the Simulator
        /// (top of the root_stack() will be used for new time).
        virtual bool                process(Simulator* s) const                 =0;
        RootStack&                  root_stack()                                { return root_stack_; }
        const RootStack&            root_stack() const                          { return root_stack_; }

        bool                        operator<(const Event& e) const
        {
            if (root_stack().empty())
                return false;
            else if (e.root_stack().empty())
                return true;
            else
                return root_stack().top() < e.root_stack().top();
        }

        virtual std::ostream&       operator<<(std::ostream& out) const
        {
            out << "Event with " << root_stack_.size() << " roots";
            if (!root_stack_.empty()) out << "; top root: " << root_stack_.top();
            out << ", ";
            return out;
        }

    private:
        RootStack                   root_stack_;
};

/**
 * Compares elements pointed at by its arguments using the provided Comparison_
 * (which must not take any arguments during construction).
 */
template<class FuncKernel_, template<class Event> class EventComparison_>
class Simulator<FuncKernel_, EventComparison_>::IndirectEventComparison:
    public std::binary_function<const typename EventComparison::first_argument_type*,
                                const typename EventComparison::second_argument_type*,
                                bool>
{
    public:
        typedef                     EventComparison                                     Comparison;
        typedef                     const typename Comparison::first_argument_type*     first_argument_type;
        typedef                     const typename Comparison::second_argument_type*    second_argument_type;

        bool                        operator()(first_argument_type e1, second_argument_type e2) const
        { return Comparison()(*e1, *e2); }
};


#include "simulator.hpp"

#endif // __SIMULATOR_H__
