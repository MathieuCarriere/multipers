#ifndef __EVENTQUEUE_H__
#define __EVENTQUEUE_H__

#include <list>
#include <vector>
#include <functional>

#include <utilities/boost.h>
#include <boost/iterator/iterator_adaptor.hpp>
namespace b = boost;


#include <utilities/log.h>
#ifdef LOGGING
static rlog::RLogChannel* rlEventQueue =             DEF_CHANNEL("utilities/eventqueue", rlog::Log_Debug);
#endif // LOGGING

#include <utilities/binaryheap.h>

#include <iostream>
#include <string>
#include <algorithm>

template<class Event_, class EventComparison_>
class EventQueue
{
    public:
        typedef                 Event_                                          Event;
        typedef                 EventComparison_                                EventComparison;

        struct                  HeapNode;
        struct                  CompareNodes;
        struct                  MarkedCompareNodes;
        struct                  UpdateNode;
        class                   iterator;
        class                   const_iterator;

        typedef                 std::list<HeapNode>                             QueueRepresentation;
        typedef                 std::vector<iterator>                           EventListHeap;

                                EventQueue()                {}

        const_iterator          top() const                 { AssertMsg(!empty(), "Queue must not be empty"); return heap_.front(); }
        iterator                top()                       { AssertMsg(!empty(), "Queue must not be empty"); return heap_.front(); }
        iterator                push(Event e);
        void                    pop();
        void                    remove(iterator i);
        void                    replace(iterator i, Event e);
        void                    promoted(iterator i);
        void                    demoted(iterator i);

        iterator                begin()                     { return queue_.begin(); }
        const_iterator          begin() const               { return queue_.begin(); }
        iterator                end()                       { return queue_.end(); }
        const_iterator          end() const                 { return queue_.end(); }
        bool                    empty() const               { return queue_.empty(); }
        size_t                  size() const                { return heap_.size(); }

        std::ostream&           print(std::ostream& out, const std::string& prefix) const;

    private:
        QueueRepresentation     queue_;
        EventListHeap           heap_;
};

template<class Event_, class EventComparison_>
struct EventQueue<Event_, EventComparison_>::HeapNode
{
                        HeapNode(const Event& e): e_(e)     {}

    size_t              heap_position_;
    Event               e_;
};

template<class Event_, class EventComparison_>
class EventQueue<Event_, EventComparison_>::iterator:
    public boost::iterator_adaptor<iterator,
                                   typename QueueRepresentation::iterator,
                                   Event>
{
    public:
        iterator():
            iterator::iterator_adaptor_(0)                              {}

        iterator(typename QueueRepresentation::iterator i):
            iterator::iterator_adaptor_(i)                              {}

        using iterator::iterator_adaptor_::base;

        typename iterator::iterator_adaptor_::reference
                        dereference() const                             { return base()->e_; }
};

template<class Event_, class EventComparison_>
class EventQueue<Event_, EventComparison_>::const_iterator:
    public boost::iterator_adaptor<const_iterator,
                                   typename QueueRepresentation::const_iterator,
                                   Event,
                                   b::use_default,
                                   const Event&>
{
    public:
        const_iterator():
            const_iterator::iterator_adaptor_(0)                        {}

        const_iterator(iterator i):
            const_iterator::iterator_adaptor_(i.base())                 {}

        const_iterator(typename QueueRepresentation::const_iterator i):
            const_iterator::iterator_adaptor_(i)                        {}

        using const_iterator::iterator_adaptor_::base;

        typename const_iterator::iterator_adaptor_::reference
                        dereference() const                             { return base()->e_; }
};

template<class Event_, class EventComparison_>
struct EventQueue<Event_, EventComparison_>::UpdateNode
{
    void operator()(iterator i, size_t pos) const                       { i.base()->heap_position_ = pos; }
};

template<class Event_, class EventComparison_>
struct EventQueue<Event_, EventComparison_>::CompareNodes
{
    bool operator()(iterator i, iterator j) const                       { return EventComparison()(*j, *i); }
};

template<class Event_, class EventComparison_>
struct EventQueue<Event_, EventComparison_>::MarkedCompareNodes
{
         MarkedCompareNodes(iterator i): i_(i)                          {}
    bool operator()(iterator i, iterator j) const
    {
        // i_ is less than everything else
        if (i == i_)
            return false;
        if (j == i_)
            return true;
        return EventComparison()(*j, *i);
    }
    iterator i_;
};


template<class Event_, class EventComparison_>
typename EventQueue<Event_, EventComparison_>::iterator
EventQueue<Event_, EventComparison_>::
push(Event e)
{
    queue_.push_back(e);
    iterator i = b::prior(queue_.end());
    heap_.push_back(i);
    std::push_heap(heap_.begin(), heap_.end(), CompareNodes(), UpdateNode());

    return i;
}

template<class Event_, class EventComparison_>
void
EventQueue<Event_, EventComparison_>::
pop()
{
    AssertMsg(!empty(), "Queue must not be empty");
    std::pop_heap(heap_.begin(), heap_.end(), CompareNodes(), UpdateNode());
    queue_.erase(heap_.back().base()); heap_.pop_back();
}

template<class Event_, class EventComparison_>
void
EventQueue<Event_, EventComparison_>::
remove(iterator i)
{
    MarkedCompareNodes mcmp(i);
    std::update_heap_pos(heap_.begin(), heap_.end(), heap_.begin() + i.base()->heap_position_, mcmp, UpdateNode());
    std::pop_heap(heap_.begin(), heap_.end(), mcmp, UpdateNode());
    AssertMsg(heap_.back() == i, "i should be in the back");
    queue_.erase(heap_.back().base());
    heap_.pop_back();
}

template<class Event_, class EventComparison_>
void
EventQueue<Event_, EventComparison_>::
replace(iterator i, Event e)
{

    if (EventComparison()(*i, e))
    {
        *i = e;
        promoted(i);
    } else
    {
        *i = e;
        demoted(i);
    }
}

template<class Event_, class EventComparison_>
void
EventQueue<Event_, EventComparison_>::
promoted(iterator i)
{
    std::update_heap_pos(heap_.begin(), heap_.end(), heap_.begin() + i.base()->heap_position_, CompareNodes(), UpdateNode());
}

template<class Event_, class EventComparison_>
void
EventQueue<Event_, EventComparison_>::
demoted(iterator i)
{
    // Bring to front
    MarkedCompareNodes mcmp(i);
    std::update_heap_pos(heap_.begin(), heap_.end(), heap_.begin() + i.base()->heap_position_, mcmp, UpdateNode());

    // Bring to back
    std::pop_heap(heap_.begin(), heap_.end(), mcmp, UpdateNode());

    // Find new place
    std::update_heap_pos(heap_.begin(), heap_.end(), heap_.begin() + i.base()->heap_position_, CompareNodes(), UpdateNode());
}

template<class Event_, class EventComparison_>
std::ostream&
EventQueue<Event_, EventComparison_>::
print(std::ostream& out, const std::string& prefix) const
{
    for (typename EventListHeap::const_iterator cur = heap_.begin(); cur != heap_.end(); ++cur)
        out << prefix << **cur << std::endl;
    return out;
}

#endif // __EVENTQUEUE_H__
