#ifndef __ORDER_H__
#define __ORDER_H__

#include "utilities/types.h"
#include "utilities/indirect.h"
#include "utilities/property-maps.h"

#include <vector>
#include <list>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
namespace bmi = boost::multi_index;

//#include <iostream>
#include <sstream>
#include <string>


/* Tags */
// TODO: pollutes global namespace; find a place to localize
struct  order           {};
struct  consistency     {};

template<class Container>
class OffsetOutputMap
{
    public:
        typedef             const typename Container::value_type*   const_element_pointer;
        typedef             typename Container::iterator            iterator;
            
                            OffsetOutputMap(const Container& order):
                                om(order.begin(),0), order_(order)  {}

        // returns a string with (i - bg_)                                
        std::string         operator()(iterator i) const                       
        { 
            std::stringstream s; 
            s << om[i];
            return s.str();
        }
        
        std::string         operator()(const_element_pointer p) const
        { 
            return (*this)(order_.iterator_to(*p));
        }

    private:
        OffsetMap<typename Container::iterator, unsigned>           om;
        const Container&                                            order_;
};

template<class Element_ = Empty<> >
struct OrderContainer
{
    typedef             Element_                                                        Element;
    typedef             boost::multi_index_container<Element,
                                                     bmi::indexed_by<
                                                        bmi::random_access<bmi::tag<order> >             /* order index */
                                                     > 
                                                    >                                   Container;
    typedef             typename Container::template index<order>::type                 OrderedContainer;

    typedef             OffsetOutputMap<Container>                                      OutputMap;

    template<class U> struct rebind
    { typedef           OrderContainer<U>                                               other; };
};


template<class Container_, class Comparison_>
struct  ElementComparison: public std::binary_function<const typename Container_::value_type*,
                                                       const typename Container_::value_type*,
                                                       bool>
{
    typedef             Container_                                                      Container;
    typedef             Comparison_                                                     Comparison;
    typedef             typename Container::value_type                                  Element;

                        ElementComparison(const Container& container, 
                                          const Comparison& cmp = Comparison()): 
                            container_(container), cmp_(cmp)                            {}
    
    bool                operator()(const Element* a, const Element* b) const            { return cmp_(container_.iterator_to(*a), container_.iterator_to(*b)); }
    
    const Container&    container_;
    const Comparison&   cmp_;
};



template<class Element_ = Empty<> >
struct OrderConsistencyContainer
{
    typedef             Element_                                                        Element;
    typedef             boost::multi_index_container<Element,
                                                     bmi::indexed_by<
                                                        bmi::random_access<bmi::tag<order> >,            /* current index */
                                                        bmi::random_access<bmi::tag<consistency> >       /* original index */
                                                     > 
                                                    >                                   Container;
    
    typedef             typename Container::template index<order>::type                 OrderedContainer;
    typedef             typename Container::template index<consistency>::type           ConsistentContainer;
    
    typedef             OffsetOutputMap<Container>                                      OutputMap;
    
    template<class U> struct rebind
    { typedef           OrderConsistencyContainer<U>                                    other; };
};


#endif // __ORDER_H__
