/**
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2005 -- 2010
 */

#ifndef __LSVINEYARD_H__
#define __LSVINEYARD_H__

#include <iostream>

#include "topology/simplex.h"
#include "topology/dynamic-persistence.h"
#include "topology/lowerstarfiltration.h"
#include "topology/vineyard.h"

#include <utilities/indirect.h>

#include <geometry/simulator.h>
#include <geometry/kinetic-sort.h>
#include <geometry/linear-kernel.h>

#include <boost/tuple/tuple.hpp>
namespace b  = boost;


template<class Vertex_, class VertexEvaluator_, class Simplex_ = Simplex<Vertex_>, class Filtration_ = Filtration<Simplex_> >
class LSVineyard
{
    public:
        typedef                     LSVineyard                                          Self;

        typedef                     Vertex_                                             Vertex;
        typedef                     VertexEvaluator_                                    VertexEvaluator;
        typedef                     typename VertexEvaluator::result_type               VertexValue;

        typedef                     Simplex_                                            Simplex;
        typedef                     Filtration_                                         LSFiltration;
        typedef                     typename LSFiltration::Index                        LSFIndex;

        typedef                     LinearKernel<VertexValue>                           KineticKernel;
        typedef                     Simulator<KineticKernel>                            KineticSimulator;
        class                       KineticVertexType;
        class                       KineticVertexComparison;
        class                       TrajectoryExtractor;
        typedef                     typename OrderContainer<KineticVertexType>::Container
                                                                                        VertexContainer;
        typedef                     typename VertexContainer::iterator                  VertexIndex;

        struct                      AttachmentData: public VineData
        {
            void                    set_attachment(VertexIndex v)                       { attachment = v; }
            VertexIndex             attachment;
        };
        typedef                     DynamicPersistenceTrails<AttachmentData>            Persistence;
        typedef                     typename Persistence::OrderIndex                    Index;
        typedef                     typename Persistence::iterator                      iterator;

        typedef                     typename Persistence::template SimplexMap<LSFiltration>
                                                                                        PFMap;

        class                       Evaluator;
        class                       StaticEvaluator;
        class                       KineticEvaluator;
        class                       DimensionFromIterator;

        typedef                     std::map<Vertex, LSFIndex>                          VertexLSFIndexMap;
        typedef                     ThroughEvaluatorComparison<VertexEvaluator>         VertexComparison;
        class                       VertexAttachmentComparison;
        typedef                     MaxVertexComparison<Simplex, VertexComparison>      SimplexComparison;

        class                       TranspositionVisitor;
        friend class                TranspositionVisitor;

        typedef                     Vineyard<Index, iterator, Evaluator>                Vnrd;

    public:
        template<class VertexIterator>
                                    LSVineyard(VertexIterator begin, VertexIterator end,
                                               LSFiltration& filtration,
                                               const VertexEvaluator& veval = VertexEvaluator());
                                    ~LSVineyard();

        void                        compute_vineyard(const VertexEvaluator& veval);
        bool                        transpose_vertices(VertexIndex vi);

        const LSFiltration&         filtration() const                                  { return filtration_; }
        const Vnrd&                 vineyard() const                                    { return vineyard_; }
        const Persistence&          persistence() const                                 { return persistence_; }
        const VertexComparison&     vertex_comparison() const                           { return vcmp_; }
        const VertexEvaluator&      vertex_evaluator() const                            { return veval_; }
        const SimplexComparison&    simplex_comparison() const                          { return scmp_; }

        VertexValue                 vertex_value(const Vertex& v) const                 { return veval_(v); }
        VertexValue                 simplex_value(const Simplex& s) const               { return vertex_value(*std::max_element(s.vertices().begin(), s.vertices().end(), vcmp_)); }
        const Simplex&              pfmap(iterator i) const                             { return pfmap_[i]; }
        const Simplex&              pfmap(Index i) const                                { return pfmap_[i]; }
        VertexIndex                 filtration_attachment(LSFIndex i) const             { return (persistence().begin() + (i - filtration().begin()))->attachment; }

        Index                       index(iterator i) const                             { return persistence_.index(i); }

    public:
        // For Kinetic Sort
        void                        swap(VertexIndex a, KineticSimulator* simulator);

    private:
        void                        change_evaluator(Evaluator* eval);
        void                        set_attachment(iterator i, VertexIndex vi)          { persistence_.modifier()(i, boost::bind(&AttachmentData::set_attachment, bl::_1, vi)); }
        void                        transpose_filtration(iterator i)                    { filtration_.transpose(filtration_.begin() + (i - persistence_.begin())); }

        bool                        verify_pairing() const;

        typedef                     b::tuple< b::reference_wrapper<const Simplex>,
                                              b::reference_wrapper<const typename Persistence::Element> >
                                                                                        SimplexPersistenceElementTuple;
        struct                      AttachmentCmp;

    private:
        VertexContainer             vertices_;

        VertexEvaluator             veval_;
        VertexComparison            vcmp_;
        SimplexComparison           scmp_;

        LSFiltration&               filtration_;
        Persistence                 persistence_;
        PFMap                       pfmap_;

        Vnrd                        vineyard_;
        Evaluator*                  evaluator_;
        unsigned                    time_count_;

#if 0
    private:
        // Serialization
        friend class boost::serialization::access;

        LSVineyard()                                                                    {}

        template<class Archive>
        void serialize(Archive& ar, version_type )
        {
            ar & BOOST_SERIALIZATION_NVP(grid_stack_);
            ar & BOOST_SERIALIZATION_NVP(vertices_);
            ar & BOOST_SERIALIZATION_NVP(filtration_);
        };
#endif
};

//BOOST_CLASS_EXPORT(LSVineyard)

template<class V, class VE, class S, class C>
std::ostream&
operator<<(std::ostream& out, const typename LSVineyard<V,VE,S,C>::VertexIndex& vi)
{ return out << vi->vertex(); }

template<class V, class VE, class S, class C>
std::ostream&
operator<<(std::ostream& out, const typename LSVineyard<V,VE,S,C>::KineticVertexType& v)
{ return out << v.vertex(); }

template<class V, class VE, class S, class C>
class LSVineyard<V,VE,S,C>::KineticVertexType
{
    public:
                                KineticVertexType(const Vertex& v):
                                    vertex_(v)                                              {}

        Vertex                  vertex() const                                              { return vertex_; }
        void                    set_vertex(Vertex v)                                        { vertex_ = v; }

        LSFIndex                simplex_index() const                                       { return simplex_index_; }
        void                    set_simplex_index(LSFIndex i)                               { simplex_index_ = i; }

    private:
        Vertex                  vertex_;
        LSFIndex                simplex_index_;
};

template<class V, class VE, class S, class C>
class LSVineyard<V,VE,S,C>::TrajectoryExtractor: public std::unary_function<VertexIndex, typename KineticSimulator::Function>
{
    public:
        typedef                 typename KineticSimulator::Function                         Function;

                                TrajectoryExtractor(const VertexEvaluator& veval0,
                                                    const VertexEvaluator& veval1):
                                    veval0_(veval0), veval1_(veval1)                        {}


        Function                operator()(VertexIndex i) const                             { VertexValue v0 = veval0_(i->vertex()), v1 = veval1_(i->vertex()); return Function(v0, v1 - v0); }

    private:
        const VertexEvaluator&  veval0_, veval1_;
};

template<class V, class VE, class S, class C>
class LSVineyard<V,VE,S,C>::KineticVertexComparison: public std::binary_function<const KineticVertexType&, const KineticVertexType&, bool>
{
    public:
                                KineticVertexComparison(const VertexComparison& vcmp):
                                    vcmp_(vcmp)                                             {}

        bool                    operator()(const KineticVertexType& v1, const KineticVertexType& v2) const
        { return vcmp_(v1.vertex(), v2.vertex()); }

    private:
        VertexComparison            vcmp_;
};

template<class V, class VE, class S, class C>
class LSVineyard<V,VE,S,C>::TranspositionVisitor: public Persistence::TranspositionVisitor
{
    public:
        typedef                 typename Persistence::TranspositionVisitor                  Parent;
        typedef                 typename LSVineyard<V,VE,S,C>::iterator                     iterator;
        typedef                 typename LSVineyard<V,VE,S,C>::Index                        Index;

                                TranspositionVisitor(LSVineyard& v): lsvineyard_(v)         {}

        void                    transpose(iterator i)                                       { lsvineyard_.transpose_filtration(i); }
        void                    switched(iterator i, SwitchType type)                       { lsvineyard_.vineyard_.switched(index(i), index(boost::next(i))); }

    private:
        Index                   index(iterator i)                                           { return lsvineyard_.index(i); }

        LSVineyard&             lsvineyard_;
};

template<class V, class VE, class S, class C>
class LSVineyard<V,VE,S,C>::Evaluator: public std::unary_function<Index, RealType>
{
    public:
        virtual ~Evaluator() {}
        virtual RealType        time() const                                                =0;
        virtual RealType        operator()(Index i) const                                   =0;
        virtual Dimension       dimension(Index i) const                                    =0;
        virtual RealType        operator()(iterator i) const                                { return operator()(&*i); }
        virtual Dimension       dimension(iterator i) const                                 { return dimension(&*i); }
};

template<class V, class VE, class S, class C>
class LSVineyard<V,VE,S,C>::DimensionFromIterator: std::unary_function<iterator, Dimension>
{
    public:
                                DimensionFromIterator(const PFMap& pfmap): pfmap_(pfmap)    {}

        Dimension               operator()(iterator i) const                                { return pfmap_[i].dimension(); }

    private:
        const PFMap&            pfmap_;
};

#include "lsvineyard.hpp"

#endif // __LSVINEYARD_H__
