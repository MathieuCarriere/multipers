#include <utilities/log.h>
#include <utilities/containers.h>
#include <utilities/property-maps.h>

#include <boost/utility/enable_if.hpp>
#include <utilities/boost.h>
#include <boost/foreach.hpp>

#include <boost/foreach.hpp>

#ifdef LOGGING
static rlog::RLogChannel* rlPersistence =                   DEF_CHANNEL("topology/persistence", rlog::Log_Debug);
#endif // LOGGING

#ifdef COUNTERS
static Counter*  cPersistencePair =                         GetCounter("persistence/pair");
static Counter*  cPersistencePairBoundaries =               GetCounter("persistence/pair/boundaries");
static Counter*  cPersistencePairCycleLength =              GetCounter("persistence/pair/cyclelength");
#endif // COUNTERS

template<class D, class CT, class OT, class E, class Cmp>
template<class Filtration>
void
StaticPersistence<D, CT, OT, E, Cmp>::
initialize(const Filtration& filtration)
{ 
    order_.assign(filtration.size(), OrderElement());
    rLog(rlPersistence, "Initializing persistence");

    OffsetMap<typename Filtration::Index, iterator>                         om(filtration.begin(), begin());
    for (typename Filtration::Index cur = filtration.begin(); cur != filtration.end(); ++cur)
    {
        Cycle z;   
        BOOST_FOREACH(const typename Filtration::Simplex& s, std::make_pair(cur->boundary_begin(), cur->boundary_end()))
            z.push_back(index(om[filtration.find(s)]));
        z.sort(ocmp_); 

        iterator ocur = om[cur];
        swap_cycle(ocur, z);
        set_pair(ocur,   ocur);
    }
}
                    
template<class D, class CT, class OT, class E, class Cmp>
void 
StaticPersistence<D, CT, OT, E, Cmp>::
pair_simplices(bool progress)
{ 
    if (progress) 
    {
        PairVisitor visitor(size());
        pair_simplices<PairVisitor>(begin(), end(), false, visitor); 
    }
    else
    {
        PairVisitorNoProgress visitor;
        pair_simplices<PairVisitorNoProgress>(begin(), end(), false, visitor);
    }
}

template<class D, class CT, class OT, class E, class Cmp>
template<class Visitor>
void 
StaticPersistence<D, CT, OT, E, Cmp>::
pair_simplices(iterator bg, iterator end, bool store_negative, const Visitor& visitor)
{
#if LOGGING
    typename ContainerTraits::OutputMap outmap(order_);
#endif

    // FIXME: need sane output for logging
    rLog(rlPersistence, "Entered: pair_simplices");
    for (iterator j = bg; j != end; ++j)
    {
        visitor.init(j);
        rLog(rlPersistence, "Simplex %s", outmap(j).c_str());

        Cycle z;
        swap_cycle(j, z);
        rLog(rlPersistence, "  has boundary: %s", z.tostring(outmap).c_str());

        // Sparsify the cycle by removing the negative elements
        if (!store_negative)
        {
            typename OrderElement::Cycle zz;
            BOOST_FOREACH(OrderIndex i, z)
                if (i->sign())           // positive
                    zz.push_back(i);
            z.swap(zz);
        }
        // --------------------------
        
        CountNum(cPersistencePairBoundaries, z.size());
        Count(cPersistencePair);

        while(!z.empty())
        {
            OrderIndex i = z.top(ocmp_);            // take the youngest element with respect to the OrderComparison
            rLog(rlPersistence, "  %s: %s", outmap(i).c_str(), outmap(i->pair).c_str());
            // TODO: is this even a meaningful assert?
            AssertMsg(!ocmp_(i, index(j)), 
                      "Simplices in the cycle must precede current simplex: (%s in cycle of %s)",
                      outmap(i).c_str(), outmap(j).c_str());

            // i is not paired, so we pair j with i
            if (iterator_to(i->pair) == iterator_to(i))
            {
                rLog(rlPersistence, "  Pairing %s and %s with cycle %s", 
                                   outmap(i).c_str(), outmap(j).c_str(), 
                                   z.tostring(outmap).c_str());
             
                set_pair(i, j);
                swap_cycle(j, z);
                set_pair(j, i);
                
                CountNum(cPersistencePairCycleLength,   j->cycle.size());
                CountBy (cPersistencePairCycleLength,   j->cycle.size());
                break;
            }

            // update element
            z.add(i->pair->cycle, ocmp_);
            visitor.update(j, iterator_to(i));
            rLog(rlPersistence, "    new cycle: %s", z.tostring(outmap).c_str());
        }
        // if z was empty, so is (already) j->cycle, so nothing to do
        visitor.finished(j);
        rLog(rlPersistence, "Finished with %s: %s", 
                            outmap(j).c_str(), outmap(j->pair).c_str());
    }
}

template<class D, class CT, class OT, class E, class Cmp>
template<class Filtration>
class StaticPersistence<D, CT, OT, E, Cmp>::SimplexMap
{
    public:
        typedef                             OrderIndex                              key_type;
        // typedef                             iterator                                key_type;
        typedef                             const typename Filtration::Simplex&     value_type;


                                            SimplexMap(const StaticPersistence& persistence,
                                                       const Filtration&        filtration):
                                                persistence_(persistence), 
                                                filtration_(filtration)             {}

        value_type                          operator[](OrderIndex k) const          { return (*this)[persistence_.iterator_to(k)]; }
        value_type                          operator[](iterator i) const            { return filtration_.simplex(filtration_.begin() + (i - persistence_.begin())); } 

    private:
        const StaticPersistence&            persistence_;
        const Filtration&                   filtration_;
};

/* Private */
