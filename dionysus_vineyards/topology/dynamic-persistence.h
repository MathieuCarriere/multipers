#ifndef __DYNAMIC_PERSISTENCE_H__
#define __DYNAMIC_PERSISTENCE_H__

#include "static-persistence.h"
#include <utilities/types.h>

#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

#ifdef COUNTERS
static Counter*  cTrailLength =             GetCounter("persistence/pair/traillength");     // the size of matrix U in RU decomposition
static Counter*  cChainLength =             GetCounter("persistence/pair/chainlength");     // the size of matrix V in R=DV decomposition
#endif // COUNTERS

template<class Data_, class ChainTraits_>
struct TrailData: public PairCycleData<Data_, ChainTraits_, TrailData<Data_, ChainTraits_> >
{
    typedef     Data_                                                                   Data;

    typedef     PairCycleData<Data_, ChainTraits_, TrailData>                           Parent;
    typedef     TrailData<Data_, ChainTraits_>                                          Self;

    typedef     typename Parent::Index                                                  Index;
    typedef     typename Parent::Cycle                                                  Cycle;
    typedef     typename Parent::Chain                                                  Chain;
    typedef     Chain                                                                   Trail;

    // Modifiers
    template<class Cmp>
    void        trail_append(Index i, const Cmp& cmp)                                   { trail.append(i, cmp); }
    template<class Cmp>
    void        trail_add(const Trail& t, const Cmp& cmp)                               { trail.add(t, cmp); }

    template<class Cmp>
    void        cycle_add(const Cycle& z, const Cmp& cmp)                               { cycle.add(z, cmp); }

    using       Parent::cycle;
    Trail                                                                               trail;
};

/**
 * Class: DynamicPersistenceTrails
 * Derives from StaticPersistence and allows one to update persistence
 * after a transposition of two contiguous simplices in a filtration.
 * In addition to reduced cycles, it stores each OrderElement's trails,
 * i.e. in addition to matrix R, it stores matrix U in vineyard notation.
 *
 * Template parameters:
 *   Data_ -                auxilliary contents to store with each OrderElement
 *   OrderDescriptor_ -     class describing how the order is stored; it defaults to <VectorOrderDescriptor>
 *                          which serves as a prototypical class
 */
// TODO: perhaps Consistency should be wrapped into a ConsistencyDescriptor that somehow knows how to initialize it.
// That way one could provide a simple consistency descriptor that just stored some integers describing the original
// position, or one could provide consistency that is references into the complex
template<class Data_ =                  Empty<>,
         class ChainTraits_ =           VectorChains<>,
         class ContainerTraits_ =       OrderConsistencyContainer<>,
         class Element_ =               TrailData<Data_, ChainTraits_>,
         class Comparison_ =            ElementComparison<typename ContainerTraits_::template rebind<Element_>::other::Container,
                                                          std::greater<typename ContainerTraits_::template rebind<Element_>::other::Container::iterator> >,
         class ConsistencyComparison_ = ElementComparison<typename ContainerTraits_::template rebind<Element_>::other::ConsistentContainer,
                                                          std::greater<typename ContainerTraits_::template rebind<Element_>::other::ConsistentContainer::iterator> >
        >
class DynamicPersistenceTrails:
    public StaticPersistence<Data_, ChainTraits_, ContainerTraits_, Element_, Comparison_>
{
    public:
        typedef         Data_                                                           Data;
        typedef         Element_                                                        Element;
        typedef         StaticPersistence<Data_, ChainTraits_,
                                          ContainerTraits_, Element_, Comparison_>      Parent;

        typedef         typename Parent::ContainerTraits                                Traits;
        typedef         typename Parent::Order                                          Order;
        typedef         typename Traits::ConsistentContainer                            Consistency;
        typedef         typename Parent::OrderComparison                                OrderComparison;
        typedef         typename Parent::OrderIndex                                     OrderIndex;
        typedef         ConsistencyComparison_                                          ConsistencyComparison;
        typedef         typename Parent::iterator                                       iterator;

        typedef         typename Element::Trail                                         Trail;
        typedef         typename Element::Cycle                                         Cycle;

         /* Constructor: DynamicPersistenceTrails() */
                                        DynamicPersistenceTrails();

        /**
         * Constructor: DynamicPersistenceTrails()
         * TODO: write a description
         *
         * Template parameters:
         *   Filtration -           filtration of the complex whose persistence we are computing
         */
        template<class Filtration>      DynamicPersistenceTrails(const Filtration& f);

        template<class Filtration>
        void                            initialize(const Filtration& f)                 { Parent::initialize(f); }

        void                            pair_simplices();

        // Function: transpose(i)
        // Tranpose i and the next element.
        // Returns: true iff the pairing switched.
        template<class DimensionFunctor, class Visitor>
        bool                            transpose(iterator i, const DimensionFunctor& dimension, Visitor visitor = Visitor());

        template<class DimensionFunctor>
        bool                            transpose(iterator i, const DimensionFunctor& dimension)    { return transpose(i,dimension,TranspositionVisitor()); }

        using                           Parent::begin;
        using                           Parent::end;
        using                           Parent::iterator_to;
        using                           Parent::index;
        using                           Parent::size;
        using                           Parent::order_comparison;

        template<class Iter>
        void                            rearrange(Iter i);

        // Struct: TranspositionVisitor
        //
        // For example, a VineardVisitor could implement this archetype.
        struct TranspositionVisitor
        {
            // Function: transpose(i)
            // This function is called before transposition is processed
            // (at the very beginning of <transpose(i, visitor)>). It is meant to update any structures
            // that may need to be updated, but perhaps it has other uses as well.
            void                        transpose(iterator i) const                     {}

            // Function: switched(i, type)
            // This function is called after the transposition if the switch in pairing has occured.
            // `i` is the index of the preceding simplex after the transposition.
            // `type` indicates the <SwitchType>.
            void                        switched(iterator i, SwitchType type) const     {}
        };

    protected:
        using                           Parent::order;
        using                           Parent::set_pair;
        using                           Parent::swap_cycle;

        Consistency&                    consistent_order()                              { return order().template get<consistency>(); }
        const Consistency&              consistent_order() const                        { return order().template get<consistency>(); }

        bool                            trail_remove_if_contains
                                            (iterator i, OrderIndex j)                  { TrailRemover rm(j); order().modify(i, rm); return rm.result; }
        void                            cycle_add(iterator i, const Cycle& z)           { order().modify(i, boost::bind(&Element::template cycle_add<ConsistencyComparison>, bl::_1, boost::ref(z), ccmp_)); }      // i->cycle_add(z, ccmp_)
        void                            trail_add(iterator i, const Trail& t)           { order().modify(i, boost::bind(&Element::template trail_add<ConsistencyComparison>, bl::_1, boost::ref(t), ccmp_)); }      // i->trail_add(t, ccmp_)

    private:
        void                            swap(iterator i, iterator j);
        void                            pairing_switch(iterator i, iterator j);

        struct PairingTrailsVisitor: public Parent::PairVisitor
        {
                                        PairingTrailsVisitor(Order& order, ConsistencyComparison ccmp, unsigned size):
                                            Parent::PairVisitor(size), order_(order), ccmp_(ccmp)   {}

            void                        init(iterator i) const                          { order_.modify(i,                                  boost::bind(&Element::template trail_append<ConsistencyComparison>, bl::_1, &*i, ccmp_)); Count(cTrailLength); }        // i->trail_append(&*i, ccmp)
            void                        update(iterator j, iterator i) const            { order_.modify(order_.iterator_to(*(i->pair)),     boost::bind(&Element::template trail_append<ConsistencyComparison>, bl::_1, &*j, ccmp_)); Count(cTrailLength); }        // i->pair->trail_append(&*j, ccmp)
            void                        finished(iterator i) const                      { Parent::PairVisitor::finished(i); }

            Order&                      order_;
            ConsistencyComparison       ccmp_;
        };

        struct TrailRemover;

        ConsistencyComparison           ccmp_;
};

/* Chains */
template<class Data_, class ChainTraits_>
struct ChainData: public PairCycleData<Data_, ChainTraits_, ChainData<Data_, ChainTraits_> >
{
    typedef     Data_                                                                   Data;

    typedef     PairCycleData<Data_, ChainTraits_, ChainData>                           Parent;
    typedef     ChainData<Data_, ChainTraits_>                                          Self;

    typedef     typename Parent::Index                                                  Index;
    typedef     typename Parent::Cycle                                                  Cycle;
    typedef     typename Parent::Chain                                                  Chain;
    typedef     Chain                                                                   Trail;

    // Modifiers
    template<class Cmp>
    void        chain_append(Index i, const Cmp& cmp)                                   { chain.append(i, cmp); }
    template<class Cmp>
    void        chain_add(const Chain& c, const Cmp& cmp)                               { chain.add(c, cmp); }

    template<class Cmp>
    void        cycle_add(const Cycle& z, const Cmp& cmp)                               { cycle.add(z, cmp); }

    using       Parent::cycle;
    Chain                                                                               chain;
};

/**
 * Class: DynamicPersistenceChains
 *
 * TODO: below comment is incorrect; nothing dynamic about this yet.
 * Derives from StaticPersistence and allows one to update persistence
 * after a transposition of two contiguous simplices in a filtration.
 * In addition to reduced cycles, it stores each OrderElement's chains,
 * i.e. in addition to matrix R, it stores matrix V in vineyard notation.
 *
 * Template parameters:
 *   Data_ -                auxilliary contents to store with each OrderElement
 *   OrderDescriptor_ -     class describing how the order is stored; it defaults to <VectorOrderDescriptor>
 *                          which serves as a prototypical class
 */
template<class Data_ =                  Empty<>,
         class ChainTraits_ =           VectorChains<>,
         class ContainerTraits_ =       OrderConsistencyContainer<>,
         class Element_ =               ChainData<Data_, ChainTraits_>,
         class Comparison_ =            ElementComparison<typename ContainerTraits_::template rebind<Element_>::other::Container,
                                                          std::greater<typename ContainerTraits_::template rebind<Element_>::other::Container::iterator> >,
         class ConsistencyComparison_ = ElementComparison<typename ContainerTraits_::template rebind<Element_>::other::ConsistentContainer,
                                                          std::greater<typename ContainerTraits_::template rebind<Element_>::other::ConsistentContainer::iterator> >
        >
class DynamicPersistenceChains:
    public StaticPersistence<Data_, ChainTraits_, ContainerTraits_, Element_, Comparison_>
{
    public:
        typedef         Data_                                                           Data;
        typedef         Element_                                                        Element;
        typedef         StaticPersistence<Data_, ChainTraits_,
                                          ContainerTraits_, Element_, Comparison_>      Parent;

        typedef         typename Parent::ContainerTraits                                Traits;
        typedef         typename Parent::Order                                          Order;

        typedef         typename Parent::OrderComparison                                OrderComparison;
        typedef         typename Parent::OrderIndex                                     OrderIndex;
        typedef         ConsistencyComparison_                                          ConsistencyComparison;
        typedef         typename Parent::iterator                                       iterator;

        typedef         typename Element::Chain                                         Chain;
        typedef         typename Element::Cycle                                         Cycle;

        /* Constructor: DynamicPersistenceChains() */
                                        DynamicPersistenceChains();

        /**
         * Constructor: DynamicPersistenceChains()
         * TODO: write a description
         *
         * Template parameters:
         *   Filtration -           filtration of the complex whose persistence we are computing
         */
        template<class Filtration>      DynamicPersistenceChains(const Filtration& f);

        template<class Filtration>
        void                            initialize(const Filtration& f)                 { Parent::initialize(f); }
        void                            pair_simplices();

        // Function: transpose(i)
        // Tranpose i and the next element.
        // Returns: true iff the pairing switched.
        // TODO
        //bool                            transpose(OrderIndex i)                         { return transpose(i, TranspositionVisitor()); }

        // TODO: the main missing piece to be dynamic
        //template<class Visitor>
        //bool                            transpose(OrderIndex i, Visitor& visitor = Visitor());

        using                           Parent::begin;
        using                           Parent::end;
        using                           Parent::iterator_to;
        using                           Parent::index;
        using                           Parent::size;

        // Struct: TranspositionVisitor
        //
        // For example, a VineardVisitor could implement this archetype.
        struct TranspositionVisitor
        {
            // Function: transpose(i)
            // This function is called before transposition is processed
            // (at the very beginning of <transpose(i, visitor)>). It is meant to update any structures
            // that may need to be updated, but perhaps it has other uses as well.
            void                        transpose(iterator i) const                     {}

            // Function: switched(i, type)
            // This function is called after the transposition if the switch in pairing has occured.
            // `i` is the index of the preceding simplex after the transposition.
            // `type` indicates the <SwitchType>.
            void                        switched(iterator i, SwitchType type) const     {}
        };

    protected:
        using                           Parent::order;
        using                           Parent::set_pair;
        using                           Parent::swap_cycle;

        void                            cycle_add(iterator i, const Cycle& z)           { order().modify(i, boost::bind(&Element::template cycle_add<ConsistencyComparison>, bl::_1, boost::ref(z), ccmp_)); }      // i->cycle_add(z, ccmp_)
        void                            chain_add(iterator i, const Chain& c)           { order().modify(i, boost::bind(&Element::template chain_add<ConsistencyComparison>, bl::_1, boost::ref(c), ccmp_)); }      // i->chain_add(c, ccmp_)

    private:
        void                            swap(OrderIndex i, OrderIndex j);
        void                            pairing_switch(OrderIndex i, OrderIndex j);

        struct PairingChainsVisitor: public Parent::PairVisitor
        {
                                        PairingChainsVisitor(Order& order, ConsistencyComparison ccmp, unsigned size):
                                            Parent::PairVisitor(size), order_(order), ccmp_(ccmp)       {}

            void                        init(iterator i) const                          { order_.modify(i,                  boost::bind(&Element::template chain_append<ConsistencyComparison>, bl::_1, &*i, ccmp_)); }                 // i->chain_append(&*i, ccmp)
            void                        update(iterator j, iterator i) const            { order_.modify(j,                  boost::bind(&Element::template chain_add<ConsistencyComparison>, bl::_1, i->pair->chain, ccmp_)); }         // j->chain.add(i->pair->chain, ccmp_)
            void                        finished(iterator i) const                      { Parent::PairVisitor::finished(i); CountBy(cChainLength, i->chain.size()); }

            Order&                      order_;
            ConsistencyComparison       ccmp_;
        };

        ConsistencyComparison           ccmp_;
};


#include "dynamic-persistence.hpp"

#endif  // __DYNAMIC_PERSISTENCE_H__
