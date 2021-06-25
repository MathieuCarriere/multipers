#ifndef __FILTRATION_H__
#define __FILTRATION_H__

#include <vector>
#include <iostream>

#include "complex-traits.h"

#include "utilities/indirect.h"
#include "utilities/property-maps.h"
#include "utilities/types.h"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/serialization.hpp>


namespace b = boost;
namespace bmi = boost::multi_index;


// Class: Filtration
//
// Filtration keeps track of the ordering of the simplices in a complex.
// The most significant function it provides is <boundary()> which converts
// the boundary of a simplex at a given index into a list of indices.
template<class Simplex_,
         class SimplexOrderIndex_ = bmi::ordered_unique<bmi::identity<Simplex_>,
                                                        typename Simplex_::VertexComparison> >
class Filtration
{
    private:
        struct                  order {};           // tag

    public:
        // Typedefs: Template parameters
        typedef                 Simplex_                                        Simplex;
        typedef                 SimplexOrderIndex_                              SimplexOrderIndex;

        typedef                 b::multi_index_container<Simplex,
                                                         bmi::indexed_by<SimplexOrderIndex,
                                                                         bmi::random_access<bmi::tag<order> >
                                                                        >
                                                        >                       Container;
        typedef                 typename Container::value_type                  value_type;

        // Typedefs: Complex and Order views
        typedef                 typename Container::template nth_index<0>::type Complex;
        typedef                 typename Container::template nth_index<1>::type Order;
        typedef                 typename Order::const_iterator                  Index;

                                Filtration()                                    {}

        // Constructor: Filtration(bg, end, cmp)
                                template<class ComplexIndex>
                                Filtration(ComplexIndex bg, ComplexIndex end):
                                    container_(bg, end)                         {}

        // Constructor: Filtration(bg, end, cmp)
                                template<class ComplexIndex, class Comparison>
                                Filtration(ComplexIndex bg, ComplexIndex end, const Comparison& cmp = Comparison()):
                                    container_(bg, end)                         { sort(cmp); }

        // Lookup
        const Simplex&          simplex(Index i) const                          { return *i; }
        Index                   find(const Simplex& s) const                    { return bmi::project<order>(container_, container_.find(s)); }

        // Modifiers
        template<class Comparison>
        void                    sort(const Comparison& cmp = Comparison())      { container_.template get<order>().sort(cmp); }
        void                    push_back(const Simplex& s)                     { container_.template get<order>().push_back(s); }
        void                    transpose(Index i)                              { container_.template get<order>().relocate(i, i+1); }
        void                    clear()                                         { container_.template get<order>().clear(); }
        template<class Iter>
        void                    rearrange(Iter i)                               { container_.template get<order>().rearrange(i); }

        Index                   begin() const                                   { return container_.template get<order>().begin(); }
        Index                   end() const                                     { return container_.template get<order>().end(); }
        size_t                  size() const                                    { return container_.size(); }

        std::ostream&           operator<<(std::ostream& out) const             { std::copy(begin(), end(), std::ostream_iterator<Simplex>(out, "\n")); return out; }

    private:
        Container               container_;

    private:
        // Serialization
        friend class                            boost::serialization::access;
        template<class Archive>
        void                                    serialize(Archive& ar, const unsigned int)
        { ar & boost::serialization::make_nvp("order", container_); }
};

template<class S, class SOI>
std::ostream&
operator<<(std::ostream& out, const Filtration<S,SOI>& f)                       { return f.operator<<(out); }


template<class Functor_, class Filtration_>
class ThroughFiltration
{
    public:
        typedef                 Filtration_                                     Filtration;
        typedef                 Functor_                                        Functor;

        typedef                 typename Functor::result_type                   result_type;
        typedef                 typename Filtration::Index                      first_argument_type;

                                ThroughFiltration(const Filtration& filtration,
                                                  const Functor&    functor):
                                    filtration_(filtration),
                                    functor_(functor)                           {}

        result_type             operator()(first_argument_type a) const         { return functor_(filtration_.simplex(a)); }

    private:
        const Filtration&       filtration_;
        const Functor&          functor_;
};

template<class Filtration, class Functor>
ThroughFiltration<Functor, Filtration>
evaluate_through_filtration(const Filtration& filtration, const Functor& functor)
{ return ThroughFiltration<Functor, Filtration>(filtration, functor); }


template<class Map, class Filtration>
class DimensionFunctor
{
    public:
                                DimensionFunctor(const Map& map, const Filtration& filtration):
                                    map_(map), filtration_(filtration)
                                {}

        template<class key_type>
        Dimension               operator()(key_type i) const                    { return filtration_.simplex(map_[i]).dimension(); }

    private:
        const Map&              map_;
        const Filtration&       filtration_;
};

template<class Map, class Filtration>
DimensionFunctor<Map, Filtration>
make_dimension_functor(const Map& map, const Filtration& filtration)
{ return DimensionFunctor<Map, Filtration>(map, filtration); }


#include "filtration.hpp"

#endif // __FILTRATION_H__
