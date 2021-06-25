#include <utilities/log.h>

#ifdef LOGGING
static rlog::RLogChannel* rlTranspositions =    DEF_CHANNEL("topology/persistence/transpositions", rlog::Log_Debug);
#endif // LOGGING

#ifdef COUNTERS
static Counter*  cTransposition =               GetCounter("persistence/transposition");
static Counter*  cTranspositionDiffDim =        GetCounter("persistence/transposition/diffdim");
static Counter*  cTranspositionCase12 =         GetCounter("persistence/transposition/case/1/2");
static Counter*  cTranspositionCase12s =        GetCounter("persistence/transposition/case/1/2/special");
static Counter*  cTranspositionCase112 =        GetCounter("persistence/transposition/case/1/1/2");
static Counter*  cTranspositionCase111 =        GetCounter("persistence/transposition/case/1/1/1");
static Counter*  cTranspositionCase22 =         GetCounter("persistence/transposition/case/2/2");
static Counter*  cTranspositionCase212 =        GetCounter("persistence/transposition/case/2/1/2");
static Counter*  cTranspositionCase211 =        GetCounter("persistence/transposition/case/2/1/1");
static Counter*  cTranspositionCase32 =         GetCounter("persistence/transposition/case/3/2");
static Counter*  cTranspositionCase31 =         GetCounter("persistence/transposition/case/3/1");
static Counter*  cTranspositionCase4 =          GetCounter("persistence/transposition/case/4");
#endif // COUNTERS


/* Trails */

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::
DynamicPersistenceTrails():
    ccmp_(consistent_order())
{}

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
template<class Filtration>
DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::
DynamicPersistenceTrails(const Filtration& f):
    Parent(f), ccmp_(consistent_order())
{}
        
template<class D, class CT, class OT, class E, class Cmp, class CCmp>
void
DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::
pair_simplices()
{ 
    PairingTrailsVisitor visitor(order(), ccmp_, size());
    Parent::pair_simplices(begin(), end(), true, visitor);
}

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
template<class DimensionFunctor, class Visitor>
bool
DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::
transpose(iterator i, const DimensionFunctor& dimension, Visitor visitor)
{
#if LOGGING
    typename Traits::OutputMap outmap(order());
#endif

    Count(cTransposition);
    typedef                 typename Element::Trail::iterator           TrailIterator;

    visitor.transpose(i);
    
    iterator i_prev = i++;

    if (dimension(i_prev) != dimension(i))
    {
        swap(i_prev, i);
        rLog(rlTranspositions, "Different dimension");
        Count(cTranspositionDiffDim);
        return false;
    }
    
    bool si = i_prev->sign(), sii = i->sign();
    if (si && sii)
    {
        rLog(rlTranspositions, "Trail prev: %s", i_prev->trail.tostring(outmap).c_str());

        // Case 1
        if (trail_remove_if_contains(i_prev, index(i)))
            rLog(rlTranspositions, "Case 1, U[i,i+1] = 1");

        iterator k = iterator_to(i_prev->pair);
        iterator l = iterator_to(i->pair);
        
        // rLog(rlTranspositions, "(i_prev, k), (i, l): (%s, %s), (%s, %s)", 
        //                         outmap(i_prev).c_str(), outmap(k).c_str(),
        //                         outmap(i).c_str(),      outmap(l).c_str());

        // Explicit treatment of unpaired simplex
        if (l == i)
        {
            swap(i_prev, i);
            rLog(rlTranspositions, "Case 1.2 --- unpaired");
            rLog(rlTranspositions, "%s", outmap(i_prev).c_str());
            Count(cTranspositionCase12);
            return false;
        } else if (k == i_prev)
        {
            if (!(l->cycle.contains(index(i_prev))))
            {
                // Case 1.2
                swap(i_prev, i);
                rLog(rlTranspositions, "Case 1.2 --- unpaired");
                rLog(rlTranspositions, outmap(i_prev).c_str());
                Count(cTranspositionCase12);
                return false;
            } else
            {
                // Case 1.2 --- special version (plain swap, but pairing switches)
                swap(i_prev, i);
                pairing_switch(i_prev, i);
                visitor.switched(i, Case12);
                rLog(rlTranspositions, "Case 1.2 --- unpaired (pairing switch)");
                rLog(rlTranspositions, outmap(i_prev).c_str());
                Count(cTranspositionCase12s);
                return true;
            }
        }
        
        rLog(rlTranspositions, "l cycle: %s", l->cycle.tostring(outmap).c_str());
        if (!(l->cycle.contains(index(i_prev))))
        {
            // Case 1.2
            rLog(rlTranspositions, "k is in l: %d", (bool) l->trail.contains(index(k)));       // if true, a special update would be needed to maintain lazy decomposition
            swap(i_prev, i);
            rLog(rlTranspositions, "Case 1.2");
            Count(cTranspositionCase12);
            return false;
        } else
        {
            // Case 1.1
            if (std::not2(order_comparison())(index(k),index(l)))
            {
                // Case 1.1.1
                swap(i_prev, i);
                cycle_add(l, k->cycle);               // Add column k to l
                trail_add(k, l->trail);               // Add row l to k
                rLog(rlTranspositions, "Case 1.1.1");
                Count(cTranspositionCase111);
                return false;
            } else
            {
                // Case 1.1.2
                swap(i_prev, i);
                cycle_add(k, l->cycle);               // Add column l to k
                trail_add(l, k->trail);               // Add row k to l
                pairing_switch(i_prev, i);
                visitor.switched(i, Case112);
                rLog(rlTranspositions, "Case 1.1.2");
                Count(cTranspositionCase112);
                return true;
            }
        }
    } else if (!si && !sii)
    {
        // Case 2
        if (!(i_prev->trail.contains(index(i))))
        {
            // Case 2.2
            swap(i_prev, i);
            rLog(rlTranspositions, "Case 2.2");
            Count(cTranspositionCase22);
            return false;
        } else
        {
            // Case 2.1
            iterator low_i =    iterator_to(i_prev->pair);
            iterator low_ii =   iterator_to(i->pair);
            trail_add(i_prev, i->trail);                   // Add row i to i_prev
            cycle_add(i, i_prev->cycle);                   // Add column i_prev to i
            swap(i_prev, i);    
            if (std::not2(order_comparison())(index(low_ii), index(low_i)))
            {
                // Case 2.1.2
                cycle_add(i_prev, i->cycle);               // Add column i to i_prev (after transposition)
                trail_add(i, i_prev->trail);               // Add row i to i_prev
                pairing_switch(i_prev, i);
                visitor.switched(i, Case212);
                rLog(rlTranspositions, "Case 2.1.2");
                Count(cTranspositionCase212);
                return true;
            } 
            
            // Case 2.1.1
            rLog(rlTranspositions, "Case 2.1.1");
            Count(cTranspositionCase211);
            return false;
        }
    } else if (!si && sii)
    {
        // Case 3
        if (!(i_prev->trail.contains(index(i))))
        {
            // Case 3.2
            swap(i_prev, i);
            rLog(rlTranspositions, "Case 3.2");
            Count(cTranspositionCase32);
            return false;
        } else
        {
            // Case 3.1
            trail_add(i_prev, i->trail);                   // Add row i to i_prev
            cycle_add(i, i_prev->cycle);                   // Add column i_prev to i
            swap(i_prev, i);
            cycle_add(i_prev, i->cycle);                   // Add column i_prev to i (after transposition)
            trail_add(i, i_prev->trail);                   // Add row i to i_prev
            pairing_switch(i_prev, i);
            visitor.switched(i, Case31);
            rLog(rlTranspositions, "Case 3.1");
            Count(cTranspositionCase31);
            return true;
        }
    } else if (si && !sii)
    {
        // Case 4
        if (trail_remove_if_contains(i_prev, index(i)))
            rLog(rlTranspositions, "Case 4, U[i,i+1] = 1");
        swap(i_prev, i);
        rLog(rlTranspositions, "Case 4");
        Count(cTranspositionCase4);
        return false;
    }
    
    return false; // to avoid compiler complaints; we should never reach this point
}


template<class D, class CT, class OT, class E, class Cmp, class CCmp>
template<class Iter>
void
DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::
rearrange(Iter i)
{ 
    order().rearrange(i); 
    consistent_order().rearrange(i);

    // Resort the cycles
    Cycle z;
    for(iterator i = begin(); i != end(); ++i)
    {
        Parent::swap_cycle(i, z);
        z.sort(ccmp_);
        Parent::swap_cycle(i, z);
    }
}

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
void
DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::
swap(iterator i, iterator j)
{
    order().relocate(i,j);
}

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
void
DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::
pairing_switch(iterator i, iterator j)
{
    OrderIndex i_pair = i->pair;
    OrderIndex j_pair = j->pair;

    // rLog(rlTranspositions, "  (%x %x %x) (%x %x %x)", &*i, i->pair, i->pair->pair, &*j, j->pair, j->pair->pair);
    if (i_pair == index(i))
        set_pair(j, j);
    else
    {
        set_pair(j, i_pair);
        set_pair(i_pair, j);
    }

    if (j_pair == index(j))
        set_pair(i, i);
    else
    {
        set_pair(i, j_pair);
        set_pair(j_pair, i);
    }
    // rLog(rlTranspositions, "  (%x %x %x) (%x %x %x)", &*i, i->pair, i->pair->pair, &*j, j->pair, j->pair->pair);
}

// Helper classes
template<class D, class CT, class OT, class E, class Cmp, class CCmp>
struct DynamicPersistenceTrails<D,CT,OT,E,Cmp,CCmp>::TrailRemover: 
    public std::unary_function<Element&, void>
{
                                TrailRemover(OrderIndex i):
                                    i_(i)                                       {}
    
    void                        operator()(Element& e)                          { result = e.trail.remove_if_contains(i_); }
    
    OrderIndex                  i_;
    bool                        result;
};


/* Chains */

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
DynamicPersistenceChains<D,CT,OT,E,Cmp,CCmp>::
DynamicPersistenceChains():
    ccmp_(order().template get<consistency>())
{}

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
template<class Filtration>
DynamicPersistenceChains<D,CT,OT,E,Cmp,CCmp>::
DynamicPersistenceChains(const Filtration& f):
    Parent(f), ccmp_(order().template get<consistency>())
{}

template<class D, class CT, class OT, class E, class Cmp, class CCmp>
void
DynamicPersistenceChains<D,CT,OT,E,Cmp,CCmp>::
pair_simplices()
{ 
    PairingChainsVisitor visitor(order(), ccmp_, size());
    Parent::pair_simplices(begin(), end(), true, visitor);
}
