#ifndef __PROPERTY_MAPS_H__
#define __PROPERTY_MAPS_H__

#include <boost/version.hpp>
#if BOOST_VERSION >= 104000
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <boost/iterator/iterator_traits.hpp>
#include <algorithm>
#include "utilities/log.h"


/* Associative Map */           // FIXME: this needs to have more thought put into it
template<class UniquePairAssociativeContainer_>
class AssociativeMap: public boost::associative_property_map<UniquePairAssociativeContainer_>
{
    public:
        typedef         boost::associative_property_map<UniquePairAssociativeContainer_>        Parent;
        
    public:
        // FIXME: takes begin, end, and initializes with the reverse
        AssociativeMap():
            container(), Parent(container)                                                      {}

    private:
        UniquePairAssociativeContainer_                                                         container;
};


/* Dereference Map */
template<class Iterator_>
class DereferenceMap
{
    public:
        typedef         Iterator_                                                               Iterator;
        typedef         Iterator                                                                key_type;
        typedef         typename boost::iterator_value<Iterator>::type                          value_type;
        typedef         boost::readable_property_map_tag                                        category;

    public:
        value_type      operator[](const key_type& k) const                                     { return *k; }
};


/* Binary Search Map */
template<class Query_, class Index_, class Comparison_ = std::less<Query_> >
class BinarySearchMap
{
    public:
        typedef         Query_                                                                  Query;
        typedef         Index_                                                                  Index;
        typedef         Comparison_                                                             Comparison;

        typedef         Query                                                                   key_type;
        typedef         Index                                                                   value_type;
        typedef         boost::readable_property_map_tag                                        category;

   public:
                        BinarySearchMap(Index bg, Index end, 
                                        Comparison cmp = Comparison()): 
                            bg_(bg), end_(end), cmp_(cmp)                                       {}

        value_type      operator[](const key_type& k) const                                     
        { 
            value_type res = std::lower_bound(bg_, end_, k, cmp_);
            AssertMsg(!cmp_(*res, k) && !cmp_(k, *res), "Query must always be found");
            return res;
        }

   private:
        Index           bg_;
        Index           end_;
        Comparison      cmp_;

};

/* Offset Map */
template<class From_, class To_>
struct OffsetMap
{
    typedef             From_                                                   From;
    typedef             To_                                                     To;
    typedef             From                                                    key_type;
    typedef             To                                                      value_type;

                        OffsetMap(From bg_from, To bg_to):
                            bg_from_(bg_from), bg_to_(bg_to)                    {}
                        
    To                  operator[](From i) const                                { return bg_to_ + (i - bg_from_); }

    From                from() const                                            { return bg_from_; }
    To                  to() const                                              { return bg_to_; }
    
    
    template<class NewFrom_> struct rebind_from
    { typedef           OffsetMap<NewFrom_, To_>                                other; };
    template<class NewTo_> struct rebind_to
    { typedef           OffsetMap<From_, NewTo_>                                other; };


    private:
                        From                                                    bg_from_;
                        To                                                      bg_to_;
};

template<class From_, class To_>
OffsetMap<From_, To_>
make_offset_map(From_ bg_from, To_ bg_to)
{ return OffsetMap<From_, To_>(bg_from, bg_to); }


template<class From_, class To_, class FromIndex_, class ToIndex_>
struct OffsetBeginMap
{
    typedef             From_                                                   From;
    typedef             To_                                                     To;
    typedef             FromIndex_                                              key_type;
    typedef             ToIndex_                                                value_type;

                        OffsetBeginMap(const From& from, const To& to):
                            from_(from), to_(to)                                {}
                        
    value_type          operator[](key_type i) const                            { return to_.begin() + (i - from_.begin()); }

    const From&         from() const                                            { return from_; }
    const To&           to() const                                              { return to_; }
    
    private:
                  const From&                                                   from_;
                  const To&                                                     to_;
};


/* ChainMap */
template<class Map1, class Map2>
class ChainMap
{
    public:
        typedef         typename Map1::key_type                                 key_type;
        typedef         typename Map2::value_type                               value_type;


                        ChainMap(const Map1& m1, const Map2& m2): 
                            m1_(m1), m2_(m2)                                    {}
        value_type      operator[](const key_type& k) const                     { return m2_[m1_[k]]; }

    private:
        const Map1&     m1_;
        const Map2&     m2_;
};

/* ThroughMap */
template<class Functor_, class Map_>
class ThroughMap
{
    public:
        typedef                 Map_                                            Map;
        typedef                 Functor_                                        Functor;

        typedef                 typename Functor::result_type                   result_type;
        typedef                 typename Map_::key_type                         first_argument_type;

                                ThroughMap(const Map&       map,
                                           const Functor&   functor):
                                    map_(map),
                                    functor_(functor)                           {}

        result_type             operator()(first_argument_type a) const         { return functor_(map_[a]); }

    private:
        const Map&              map_;
        const Functor&          functor_;
};

template<class Map, class Functor>
ThroughMap<Functor, Map>
evaluate_through_map(const Map& map, const Functor& functor)
{ return ThroughMap<Functor, Map>(map, functor); }

#endif // __PROPERTY_MAPS_H__
