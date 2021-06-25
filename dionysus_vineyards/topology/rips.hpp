#include <algorithm>
#include <utility>
#include <iostream>
#include <utilities/log.h>
#include <utilities/counter.h>
#include <utilities/indirect.h>
#include <utilities/boost.h>
#include <boost/iterator/counting_iterator.hpp>
#include <functional>

#ifdef LOGGING
static rlog::RLogChannel* rlRips =                  DEF_CHANNEL("rips/info", rlog::Log_Debug);
static rlog::RLogChannel* rlRipsDebug =             DEF_CHANNEL("rips/debug", rlog::Log_Debug);
#endif // LOGGING

#ifdef COUNTERS
static Counter*  cClique =                          GetCounter("rips/clique");
#endif // COUNTERS

template<class D, class S>
template<class Functor, class Iterator>
void
Rips<D,S>::
generate(Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    rLog(rlRipsDebug,       "Entered generate with %d indices", distances().size());

    WithinDistance neighbor(distances(), max);

    // current      = empty
    // candidates   = everything
    VertexContainer current;
    VertexContainer candidates(bg, end);
    bron_kerbosch(current, candidates, boost::prior(candidates.begin()), k, neighbor, f);
}

template<class D, class S>
template<class Functor, class Iterator>
void
Rips<D,S>::
vertex_cofaces(IndexType v, Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    WithinDistance neighbor(distances(), max);

    // current      = [v]
    // candidates   = everything - [v]
    VertexContainer current; current.push_back(v);
    VertexContainer candidates;
    for (Iterator cur = bg; cur != end; ++cur)
        if (*cur != v && neighbor(v, *cur))
            candidates.push_back(*cur);
    bron_kerbosch(current, candidates, boost::prior(candidates.begin()), k, neighbor, f);
}

template<class D, class S>
template<class Functor, class Iterator>
void
Rips<D,S>::
edge_cofaces(IndexType u, IndexType v, Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    rLog(rlRipsDebug,   "In edge_cofaces(%d, %d)", u, v);

    WithinDistance neighbor(distances(), max);
    AssertMsg(neighbor(u,v), "The new edge must be in the complex");

    // current      = [u,v]
    // candidates   = everything - [u,v]
    VertexContainer current; current.push_back(u); current.push_back(v);

    VertexContainer candidates;
    for (Iterator cur = bg; cur != end; ++cur)
        if (*cur != u && *cur != v && neighbor(v,*cur) && neighbor(u,*cur))
        {
            candidates.push_back(*cur);
            rLog(rlRipsDebug,   "  added candidate: %d", *cur);
        }

    bron_kerbosch(current, candidates, boost::prior(candidates.begin()), k, neighbor, f);
}

template<class D, class S>
template<class Functor, class Iterator>
void
Rips<D,S>::
cofaces(const Simplex& s, Dimension k, DistanceType max, const Functor& f, Iterator bg, Iterator end) const
{
    rLog(rlRipsDebug,   "In cofaces(%s)", tostring(s).c_str());

    WithinDistance neighbor(distances(), max);

    // current      = s.vertices()
    VertexContainer current(s.vertices().begin(), s.vertices().end());
    
    // candidates   = everything - s.vertices()     that is a neighbor() of every vertex in the simplex
    VertexContainer candidates;
    typedef difference_iterator<Iterator, 
                                typename VertexContainer::const_iterator, 
                                std::less<Vertex> >                     DifferenceIterator;
    for (DifferenceIterator cur =  DifferenceIterator(bg, end, s.vertices().begin(), s.vertices().end()); 
                            cur != DifferenceIterator(end, end, s.vertices().end(), s.vertices().end()); 
                            ++cur)
    {
        bool nghbr = true;
        for (typename VertexContainer::const_iterator v = s.vertices().begin(); v != s.vertices().end(); ++v)
            if (!neighbor(*v, *cur))        { nghbr = false; break; }

        if (nghbr)  
        {
            candidates.push_back(*cur);
            rLog(rlRipsDebug,   "  added candidate: %d", *cur);
        }
    }

    bron_kerbosch(current, candidates, boost::prior(candidates.begin()), k, neighbor, f, false);
}


template<class D, class S>
template<class Functor, class NeighborTest>
void
Rips<D,S>::
bron_kerbosch(VertexContainer&                          current,    
              const VertexContainer&                    candidates,     
              typename VertexContainer::const_iterator  excluded,
              Dimension                                 max_dim,    
              const NeighborTest&                       neighbor,       
              const Functor&                            functor,
              bool                                      check_initial) const
{
    rLog(rlRipsDebug,       "Entered bron_kerbosch");
    
    if (check_initial && !current.empty())
    {
        Simplex s(current);
        rLog(rlRipsDebug,   "Reporting simplex: %s", tostring(s).c_str());
        functor(s);
    }

    if (current.size() == static_cast<size_t>(max_dim) + 1) 
        return;

    rLog(rlRipsDebug,       "Traversing %d vertices", candidates.end() - boost::next(excluded));
    for (typename VertexContainer::const_iterator cur = boost::next(excluded); cur != candidates.end(); ++cur)
    {
        current.push_back(*cur);
        rLog(rlRipsDebug,   "  current.size() = %d, current.back() = %d", current.size(), current.back());

        VertexContainer new_candidates;
        for (typename VertexContainer::const_iterator ccur = candidates.begin(); ccur != cur; ++ccur)
            if (neighbor(*ccur, *cur))
                new_candidates.push_back(*ccur);
        size_t ex = new_candidates.size();
        for (typename VertexContainer::const_iterator ccur = boost::next(cur); ccur != candidates.end(); ++ccur)
            if (neighbor(*ccur, *cur))
                new_candidates.push_back(*ccur);
        excluded  = new_candidates.begin() + (ex - 1);

        bron_kerbosch(current, new_candidates, excluded, max_dim, neighbor, functor);
        current.pop_back();
    }
}

template<class Distances_, class Simplex_>
typename Rips<Distances_, Simplex_>::DistanceType
Rips<Distances_, Simplex_>::
distance(const Simplex& s1, const Simplex& s2) const
{
    DistanceType mx = 0;
    for (typename Simplex::VertexContainer::const_iterator      a = s1.vertices().begin();   a != s1.vertices().end();    ++a)
        for (typename Simplex::VertexContainer::const_iterator  b = s2.vertices().begin();   b != s2.vertices().end();    ++b)
            mx = std::max(mx, distances_(*a,*b));
    return mx;
}

template<class Distances_, class Simplex_>
typename Rips<Distances_, Simplex_>::DistanceType
Rips<Distances_, Simplex_>::
max_distance() const
{
    DistanceType mx = 0;
    for (IndexType a = distances_.begin(); a != distances_.end(); ++a)
        for (IndexType b = boost::next(a); b != distances_.end(); ++b)
            mx = std::max(mx, distances_(a,b));
    return mx;
}

template<class Distances_, class Simplex_>
typename Rips<Distances_, Simplex_>::DistanceType
Rips<Distances_, Simplex_>::Evaluator::
operator()(const Simplex& s) const
{
    DistanceType mx = 0;
    for (typename Simplex::VertexContainer::const_iterator      a = s.vertices().begin();   a != s.vertices().end();    ++a)
        for (typename Simplex::VertexContainer::const_iterator  b = boost::next(a);         b != s.vertices().end();    ++b)
            mx = std::max(mx, distances_(*a,*b));
    return mx;
}
