#ifndef __STATIC_PERSISTENCE_H__
#define __STATIC_PERSISTENCE_H__

#include "order.h"
#include "cycles.h"
#include "filtration.h"

#include <boost/ref.hpp>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

#include <utilities/types.h>

#include <boost/progress.hpp>


// Element_ should derive from PairCycleData
template<class Data_, class ChainTraits_, class Element_ = use_default>
struct PairCycleData: public Data_
{
    typedef     Data_                                                                   Data;
    typedef     typename if_default<Element_, PairCycleData>::type                      Element;

    typedef     const Element*                                                          Index;
    typedef     typename ChainTraits_::template rebind<Index>::other                    ChainTraits;
    typedef     typename ChainTraits::Chain                                             Chain;
    typedef     Chain                                                                   Cycle;

                PairCycleData(Index p = Index(), const Cycle& z = Cycle(), const Data& d = Data()):
                    Data(d), pair(p), cycle(z)
                {}
    
    bool        sign() const                                                            { return cycle.empty(); }
    bool        unpaired() const                                                        { return pair == this; }
    
    void        swap_cycle(Cycle& z)                                                    { cycle.swap(z); }
    void        set_pair(Index i)                                                       { pair = i; }

    Index       pair;
    Cycle       cycle;
};

/**
 * Class: StaticPersistence
 * The class that encapsulates order and pairing information as well as 
 * implements methods to compute and maintain the pairing. Its most 
 * significant method is <pair_simplices()>.
 *
 * Template parameters:
 *   Data_ -                auxilliary contents to store with each OrderElement
 *   OrderDescriptor_ -     class describing how the order is stored; it defaults to <VectorOrderDescriptor> 
 *                          which serves as a prototypical class
 */
template<class Data_ =              Empty<>,
         class ChainTraits_ =       VectorChains<>,
         class ContainerTraits_ =   OrderContainer<>,
         class Element_ =           PairCycleData<Data_, ChainTraits_>,
         class Comparison_ =        ElementComparison<typename ContainerTraits_::template rebind<Element_>::other::Container,
                                                      std::greater<typename ContainerTraits_::template rebind<Element_>::other::Container::iterator> > >
class StaticPersistence
{
    public:
        // Typedef: Data
        // The data type stored in each order element
        typedef                         Data_                                                   Data;

        typedef                         Element_                                                Element;
        typedef                         typename ContainerTraits_::
                                                    template rebind<Element>::other             ContainerTraits;
        typedef                         typename ContainerTraits::Container                     Container;
        typedef                         Container                                               Order;
        typedef                         typename Element::Index                                 OrderIndex;
        typedef                         Element                                                 OrderElement;
        typedef                         typename Order::iterator                                iterator;

        typedef                         typename ChainTraits_::
                                                    template rebind<OrderIndex>::other          ChainTraits;
        typedef                         typename ChainTraits::Chain                             Chain;
        typedef                         Chain                                                   Cycle;
        
        typedef                         Comparison_                                             OrderComparison;

        /* Constructor: StaticPersistence() */
                                        StaticPersistence(): ocmp_(order_)                      {}

        /* Constructor: StaticPersistence()
         * TODO: write a description
         *
         * Template parameters:
         *   Filtration -           filtration of the complex whose persistence we are computing
         */
        template<class Filtration>      StaticPersistence(const Filtration& f): ocmp_(order_)   { initialize(f); }

        // Function: initialize(const Filtration& f)
        // Initialize the boundary map from the Filtration
        template<class Filtration>
        void                            initialize(const Filtration& f);
        
        // Function: pair_simplices()                                        
        // Compute persistence of the filtration
        void                            pair_simplices(bool progress = true);

        // Functions: Accessors
        //   begin() -              returns OrderIndex of the first element
        //   end() -                returns OrderIndex of one past the last element
        //   size() -               number of elements in the StaticPersistence
        iterator                        begin() const                                           { return order_.begin(); }
        iterator                        end() const                                             { return order_.end(); }
        iterator                        iterator_to(OrderIndex i) const                         { return order_.iterator_to(*i); }
        OrderIndex                      index(iterator i) const                                 { return &*i; }
        size_t                          size() const                                            { return order_.size(); }
        const OrderComparison&          order_comparison() const                                { return ocmp_; }

        // A map to extract simplices
        template<class Filtration>      class SimplexMap;
        template<class Filtration>
        SimplexMap<Filtration>          make_simplex_map(const Filtration& filtration) const    { return SimplexMap<Filtration>(*this, filtration); }

        class                           OrderModifier
        {
            public:
                                        OrderModifier(Order& order): order_(order)              {}
                template<class Functor> 
                void                    operator()(iterator i, const Functor& f)                { order_.modify(i, f); }

            private:
                Order&                  order_;
        };
        OrderModifier                   modifier()                                              { return OrderModifier(order()); }

        // Function: pair_simplices(bg, end)
        // Compute persistence of the simplices in filtration between bg and end
        template<class Visitor>
        void                            pair_simplices(iterator bg, iterator end, bool store_negative = false, const Visitor& visitor = Visitor());

        // Struct: PairVisitor
        // Acts as an archetype and if necessary a base class for visitors passed to <pair_simplices(bg, end, visitor)>.
        struct                          PairVisitor
        {
                                        PairVisitor(unsigned size): show_progress(size)         {}
            // Function: init(i)
            // Called after OrderElement pointed to by `i` has been initialized 
            // (its cycle is set to be its boundary, and pair is set to self, i.e. `i`)
            void                        init(iterator i) const                                  {}
            
            // Function: update(j, i)
            // Called after the cycle of `i` has been added to the cycle of `j`, 
            // this allows the derived class to perform the necessary updates 
            // (e.g., add `i`'s chain to `j`'s chain)
            void                        update(iterator j, iterator i) const                    {}

            // Function: finished(j)
            // Called after the processing of `j` is finished.
            void                        finished(iterator j) const                              { ++show_progress; }
            mutable boost::progress_display     
                                        show_progress;
        };
        
        struct                          PairVisitorNoProgress
        {
                                        PairVisitorNoProgress()                                 {}
            void                        init(iterator i) const                                  {}
            void                        update(iterator j, iterator i) const                    {}
            void                        finished(iterator j) const                              {}
        };

    protected:
        const Order&                    order() const                                           { return order_; }
        Order&                          order()                                                 { return order_; }

        void                            set_pair(iterator i,    iterator j)                     { set_pair(i, &*j); }
        void                            set_pair(iterator i,    OrderIndex j)                   { order_.modify(i, boost::bind(&OrderElement::set_pair, bl::_1, j)); }                  // i->set_pair(j)
        void                            set_pair(OrderIndex i,  iterator j)                     { set_pair(iterator_to(i), &*j); }
        void                            set_pair(OrderIndex i,  OrderIndex j)                   { set_pair(iterator_to(i), j); }
        void                            swap_cycle(iterator i,  Cycle& z)                       { order_.modify(i, boost::bind(&OrderElement::swap_cycle, bl::_1, boost::ref(z))); }    // i->swap_cycle(z)

    private:
        Order                           order_;
        OrderComparison                 ocmp_;
};

#include "static-persistence.hpp"

#endif // __STATIC_PERSISTENCE_H__
