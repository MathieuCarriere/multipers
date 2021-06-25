#ifndef __CYCLES_H__
#define __CYCLES_H__

#include "chain.h"
#include "utilities/circular_list.h"

#if DEBUG_CONTAINERS
    #include <debug/vector>
    #include <debug/deque>
    using std::__debug::vector;
    using std::__debug::deque;
#else
    #include <vector>
    #include <deque>
    using std::vector;
    using std::deque;
#endif

template<class OrderIndex_ = int>
struct VectorChains
{
    typedef             OrderIndex_                                             OrderIndex;
    typedef             ChainWrapper<vector<OrderIndex> >                       Chain;
    
    template<class U> struct rebind
    { typedef           VectorChains<U>         other; };
};

template<class OrderIndex_ = int>
struct DequeChains
{
    typedef             OrderIndex_                                             OrderIndex;
    typedef             ChainWrapper<deque<OrderIndex> >                        Chain;

    template<class U> struct rebind
    { typedef           DequeChains<U>         other; };
};

template<class OrderIndex_ = int>
struct ListChains
{
    typedef             OrderIndex_                                             OrderIndex;
    typedef             ChainWrapper<List<OrderIndex> >                         Chain;
    
    template<class U> struct rebind
    { typedef           ListChains<U>           other; };
};

#endif // __CYCLES_H__
