#ifndef __PERSISTENCE_DIAGRAM_H__
#define __PERSISTENCE_DIAGRAM_H__

#include <utilities/types.h>

#include <vector>
#include <iostream>
#include <cmath>

#include <boost/compressed_pair.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>

/**
 * Class: PDPoint
 *
 * Stores birth-death pair plus any additional information provided by `Data` template parameter.
 */
template<class Data_ = Empty<> >
class PDPoint
{
    public:
        typedef                 Data_                                       Data;

                                PDPoint(const PDPoint& other):
                                    point_(other.point_)                    {}
                                PDPoint(RealType x = 0, RealType y = 0, const Data& data = Data());

        RealType                x() const                                   { return point_.first().first; }
        RealType                y() const                                   { return point_.first().second; }
        const Data&             data() const                                { return point_.second(); }
        Data&                   data()                                      { return point_.second(); }

        std::ostream&           operator<<(std::ostream& out) const         { return (out << x() << " " << y()); } // << " " << data()); }

        struct Visitor
        {
            template<class Iterator>
            void                point(Iterator i, PDPoint& p) const         {}
        };

    private:
        RealType&               x()                                         { return point_.first().first; }
        RealType&               y()                                         { return point_.first().second; }

    private:
        boost::compressed_pair<std::pair<RealType, RealType>, Data>       point_;

    private:
        /* Serialization */
        friend class boost::serialization::access;

        template<class Archive>
        void                    serialize(Archive& ar, version_type );
};

template<class Data>
std::ostream&                   operator<<(std::ostream& out, const PDPoint<Data>& point)
{ return (point.operator<<(out)); }

template<class Point, class Iterator, class Evaluator, class Visitor>
boost::optional<Point>
make_point(Iterator i, const Evaluator& evaluator, const Visitor& visitor);

template<class Point, class Iterator, class Evaluator>
boost::optional<Point>
make_point(Iterator i, const Evaluator& evaluator)
{ return make_point<Point>(i, evaluator, Point::Visitor()); }


/**
 * Class: PersistenceDiagram
 *
 * Stores birth-death pairs, i.e. points in the extended plane. Each point can also store
 * additional information described by `Data_` template parameter.
 */
template<class Data_ = Empty<> >
class PersistenceDiagram
{
    public:
        typedef                 Data_                                       Data;
        typedef                 PDPoint<Data>                               Point;
        typedef                 std::vector<Point>                          PointVector;
        typedef                 typename PointVector::const_iterator        const_iterator;

                                PersistenceDiagram()                        {}

                                PersistenceDiagram( Dimension dimension ):
                                    dimension_( dimension )                  {}

        template<class OtherData>
                                PersistenceDiagram(const PersistenceDiagram<OtherData>& other);

        template<class Iterator, class Evaluator>
                                PersistenceDiagram(Iterator bg, Iterator end,
                                                   const Evaluator& eval = Evaluator());

        template<class Iterator, class Evaluator, class Visitor>
                                PersistenceDiagram(Iterator bg, Iterator end,
                                                   const Evaluator& eval = Evaluator(),
                                                   const Visitor& visitor = Visitor());

        template<class Iterator, class Evaluator, class Visitor>
        void                    init(Iterator bg, Iterator end,
                                     const Evaluator& eval = Evaluator(),
                                     const Visitor& visitor = Visitor());

        const_iterator          begin() const                               { return points_.begin(); }
        const_iterator          end() const                                 { return points_.end(); }
        size_t                  size() const                                { return points_.size(); }

        void                    push_back(const Point& point)               { points_.push_back(point); }

        std::ostream&           operator<<(std::ostream& out) const;

        Dimension               dimension() const                           { return dimension_; }

    private:
        PointVector             points_;
        Dimension               dimension_;

    private:
        /* Serialization */
        friend class boost::serialization::access;

        template<class Archive>
        void                    serialize(Archive& ar, version_type );
};

template<class Data>
std::ostream&                   operator<<(std::ostream& out, const PersistenceDiagram<Data>& pd)
{ return (pd.operator<<(out)); }

// Function: init_diagram_vector()
template<class Diagrams, class Iterator, class Evaluator, class DimensionExtractor>
void                    init_diagrams(Diagrams& diagrams,
                                      Iterator bg, Iterator end,
                                      const Evaluator& evaluator = Evaluator(),
                                      const DimensionExtractor& dimension = DimensionExtractor());

template<class Diagrams, class Iterator, class Evaluator, class DimensionExtractor, class Visitor>
void                    init_diagrams(Diagrams& diagrams,
                                      Iterator bg, Iterator end,
                                      const Evaluator& evaluator = Evaluator(),
                                      const DimensionExtractor& dimension = DimensionExtractor(),
                                      const Visitor& visitor = Visitor());

// Class: Linfty
// Functor that computes L infinity norm between two points
template<class Point1, class Point2>
struct Linfty
{
    RealType            operator()(const Point1& p1, const Point2& p2) const        { return std::max(std::abs(p1.x() - p2.x()),
                                                                                                      std::abs(p1.y() - p2.y())); }

    template<class Point>
    RealType            diagonal(const Point& p) const                              { return std::abs(p.y() - p.x())/2; }
};

// Function: bottleneck_distance(dgm1, dgm2)
// Computes bottleneck distance between the two diagrams.
template<class Diagram1,
         class Diagram2,
         class Norm>
RealType                bottleneck_distance(const Diagram1& dgm1, const Diagram2& dgm2, const Norm& norm = Norm());

template<class Diagram1,
         class Diagram2>
RealType                bottleneck_distance(const Diagram1& dgm1, const Diagram2& dgm2)
{ return bottleneck_distance(dgm1, dgm2, Linfty<typename Diagram1::Point, typename Diagram2::Point>()); }

template<class Diagram>
RealType                wasserstein_distance(const Diagram& dgm1, const Diagram& dgm2, unsigned p);


#include "persistence-diagram.hpp"

#endif // __PERSISTENCE_DIAGRAM_H__
