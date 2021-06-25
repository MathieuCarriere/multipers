#include <boost/serialization/vector.hpp>
#include <boost/serialization/nvp.hpp>

#include "utilities/munkres/munkres.h"

using boost::serialization::make_nvp;

template<class D>
PDPoint<D>::
PDPoint(RealType x, RealType y, const Data& data)
{
    point_.first().first = x;
    point_.first().second = y;
    point_.second() = data;
}


template<class D>
template<class OtherData>
PersistenceDiagram<D>::
PersistenceDiagram(const PersistenceDiagram<OtherData>& other)
{
    points_.reserve(other.size());
    for (typename PersistenceDiagram<OtherData>::PointVector::const_iterator cur = points_.begin();
                                                                             cur != points_.end(); ++cur)
        push_back(Point(cur->x(), cur->y()));
}

template<class D>
template<class Iterator, class Evaluator>
PersistenceDiagram<D>::
PersistenceDiagram(Iterator bg, Iterator end, const Evaluator& eval)
{
    init(bg, end, eval, Point::Visitor());
}

template<class D>
template<class Iterator, class Evaluator, class Visitor>
PersistenceDiagram<D>::
PersistenceDiagram(Iterator bg, Iterator end, const Evaluator& eval, const Visitor& visitor)
{
    init(bg, end, eval, visitor);
}

template<class D>
template<class Iterator, class Evaluator, class Visitor>
void
PersistenceDiagram<D>::
init(Iterator bg, Iterator end, const Evaluator& evaluator, const Visitor& visitor)
{
    for (Iterator cur = bg; cur != end; ++cur)
        if (cur->sign())
        {
            boost::optional<Point> p = make_point(cur, evaluator, visitor);
            if (p)  push_back(*p);
        }
}

template<class Point, class Iterator, class Evaluator, class Visitor>
boost::optional<Point>
make_point(Iterator i, const Evaluator& evaluator, const Visitor& visitor)
{
    RealType x = evaluator(&*i);
    RealType y = Infinity;
    if (&*(i->pair) != &*i)
        y = evaluator(&*(i->pair));

    Point p(x,y);
    visitor.point(i, p);

    if (x == y) return boost::optional<Point>();

    return p;
}

template<class Diagrams, class Iterator, class Evaluator, class DimensionExtractor>
void    init_diagrams(Diagrams& diagrams,
                      Iterator bg, Iterator end,
                      const Evaluator& evaluator,
                      const DimensionExtractor& dimension)
{
    // FIXME: this is specialized for Diagrams that is std::map
    typedef             typename Diagrams::mapped_type              PDiagram;

    init_diagrams(diagrams, bg, end, evaluator, dimension, typename PDiagram::Point::Visitor());
}

template<class Diagrams, class Iterator, class Evaluator, class DimensionExtractor, class Visitor>
void    init_diagrams(Diagrams& diagrams,
                      Iterator bg, Iterator end,
                      const Evaluator& evaluator,
                      const DimensionExtractor& dimension,
                      const Visitor& visitor)
{
    // FIXME: this is specialized for Diagrams that is std::map
    typedef             typename Diagrams::mapped_type              PDiagram;

    for (Iterator cur = bg; cur != end; ++cur)
        if (cur->sign())
        {
            boost::optional<typename PDiagram::Point> p = make_point<typename PDiagram::Point>(cur, evaluator, visitor);
            if (p)
                diagrams[dimension(&*cur)].push_back(*p);
        }
}

template<class D>
std::ostream&
PersistenceDiagram<D>::
operator<<(std::ostream& out) const
{
    for (const_iterator cur = begin(); cur != end(); ++cur)
        out << *cur << std::endl;
    return out;
}

template<class D>
template<class Archive>
void
PDPoint<D>::
serialize(Archive& ar, version_type )
{
    ar & make_nvp("x", x());
    ar & make_nvp("y", y());
    ar & make_nvp("data", data());
}

template<class D>
template<class Archive>
void
PersistenceDiagram<D>::
serialize(Archive& ar, version_type )
{
    ar & make_nvp("points", points_);
}


/**
 * Some structures to compute bottleneck distance between two persistence diagrams (in bottleneck_distance() function below)
 * by setting up bipartite graphs, and finding maximum cardinality matchings in them using Boost Graph Library.
 */
#include <boost/iterator/counting_iterator.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

struct Edge: public std::pair<unsigned, unsigned>
{
    typedef         std::pair<unsigned, unsigned>                       Parent;

                    Edge(unsigned v1, unsigned v2, RealType d):
                        Parent(v1, v2), distance(d)                     {}

    bool            operator<(const Edge& other) const                  { return distance < other.distance; }

    RealType        distance;
};
typedef std::vector<Edge>               EdgeVector;
typedef EdgeVector::const_iterator      EV_const_iterator;

struct CardinaliyComparison
{
    typedef         boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>         Graph;
    typedef         std::vector<boost::graph_traits<Graph>::vertex_descriptor>                  MatchingVector;

                    CardinaliyComparison(unsigned size, EV_const_iterator begin):
                        max_size(size), bg(begin), last(bg), g(2*max_size), mates(2*max_size)
                    { boost::add_edge(bg->first, bg->second, g); }

    bool            operator()(EV_const_iterator i1, EV_const_iterator i2)
    {
        //std::cout << "Max size: " << max_size << std::endl;
        //std::cout << "Comparing: (" << i1->first << ", " << i1->second << ") and "
        //          <<            "(" << i2->first << ", " << i2->second << ")" << std::endl;

        // FIXME: the matching is being recomputed from scratch every time, this should be fixed
        if (i2 > last)
            do
            {
                ++last;
                boost::add_edge(last->first, last->second, g);
            } while (last != i2);
        else
            do
            {
                boost::remove_edge(last->first, last->second, g);
            } while (--last != i2);

        edmonds_maximum_cardinality_matching(g, &mates[0]);
        //std::cout << "Found matching of size: " << matching_size(g, &mates[0]) << std::endl;
        return matching_size(g, &mates[0]) == max_size;
    }

    unsigned                max_size;
    EV_const_iterator       bg;
    EV_const_iterator       last;
    Graph                   g;
    MatchingVector          mates;
};

// Bottleneck distance
template<class Diagram1, class Diagram2, class Norm>
RealType                bottleneck_distance(const Diagram1& dgm1, const Diagram2& dgm2, const Norm& norm)
{
    typedef         typename Diagram1::const_iterator                   Citer1;
    typedef         typename Diagram2::const_iterator                   Citer2;

    const unsigned  max_size = dgm1.size() + dgm2.size();

    // Compute all the edges and sort them by distance
    EdgeVector   edges;

    // Connect all diagonal points to each other
    for (unsigned i = dgm1.size(); i < max_size; ++i)
        for (unsigned j = max_size + dgm2.size(); j < 2*max_size; ++j)
            edges.push_back(Edge(i, j, 0));

    // Edges between real points
    unsigned i = 0;
    for (Citer1 cur1 = dgm1.begin(); cur1 != dgm1.end(); ++cur1)
    {
        unsigned j = max_size;
        for (Citer2 cur2 = dgm2.begin(); cur2 != dgm2.end(); ++cur2)
            edges.push_back(Edge(i,j++, norm(*cur1, *cur2)));

        ++i;
    }

    // Edges between real points and their corresponding diagonal points
    i = 0;
    for (Citer1 cur1 = dgm1.begin(); cur1 != dgm1.end(); ++cur1, ++i)
        edges.push_back(Edge(i, max_size + dgm2.size() + i, norm.diagonal(*cur1)));
    i = max_size;
    for (Citer2 cur2 = dgm2.begin(); cur2 != dgm2.end(); ++cur2, ++i)
        edges.push_back(Edge(dgm1.size() + (i - max_size), i, norm.diagonal(*cur2)));


    std::sort(edges.begin(), edges.end());
    //for (i = 0; i < edges.size(); ++i)
    //    std::cout << "Edge: " << edges[i].first << " " << edges[i].second << " " << edges[i].distance << std::endl;

    // Perform cardinality based binary search
    typedef boost::counting_iterator<EV_const_iterator>         EV_counting_const_iterator;
    EV_counting_const_iterator bdistance = std::upper_bound(EV_counting_const_iterator(edges.begin()),
                                                            EV_counting_const_iterator(edges.end()),
                                                            edges.begin(),
                                                            CardinaliyComparison(max_size, edges.begin()));

    return (*bdistance)->distance;
}

// Wasserstein distance
template<class Diagram>
RealType
wasserstein_distance(const Diagram& dgm1, const Diagram& dgm2, unsigned p)
{
    typedef         RealType                    Distance;
    typedef         typename Diagram::Point     Point;
    typedef         Linfty<Point, Point>        Norm;

    unsigned size = dgm1.size() + dgm2.size();
    Norm norm;

    // Setup the matrix
    Matrix<Distance>        m(size,size);
    for (unsigned i = 0; i < dgm1.size(); ++i)
        for (unsigned j = 0; j < dgm2.size(); ++j)
        {
            const Point& p1 = *(dgm1.begin() + i);
            const Point& p2 = *(dgm2.begin() + j);
            m(i,j) = pow(norm(p1, p2),  p);
            m(j + dgm1.size(), i + dgm2.size()) = 0;
        }

    for (unsigned i = 0; i < dgm1.size(); ++i)
        for (unsigned j = dgm2.size(); j < size; ++j)
        {
            const Point& p1 = *(dgm1.begin() + i);
            m(i,j) = pow(norm.diagonal(p1), p);
        }

    for (unsigned j = 0; j < dgm2.size(); ++j)
        for (unsigned i = dgm1.size(); i < size; ++i)
        {
            const Point& p2 = *(dgm2.begin() + j);
            m(i,j) = pow(norm.diagonal(p2), p);
        }

    // Compute weighted matching
    Munkres munkres;
    munkres.solve(m);

    // Assume everything is assigned (i.e., that we have a perfect matching)
    Distance sum = 0;
    for (unsigned i = 0; i < size; i++)
        for (unsigned j = 0; j < size; j++)
            if (m(i,j) == 0)
            {
                //std::cout << i << ": " << j << '\n';
                //sum += m[i][j];
                if (i >= dgm1.size())
                {
                    if (j >= dgm2.size())
                        sum += 0;
                    else
                    {
                        const Point& p2 = *(dgm2.begin() + j);
                        sum += pow(norm.diagonal(p2), p);
                    }
                } else
                {
                    if (j >= dgm2.size())
                    {
                        const Point& p1 = *(dgm1.begin() + i);
                        sum += pow(norm.diagonal(p1), p);
                    } else
                    {
                        const Point& p1 = *(dgm1.begin() + i);
                        const Point& p2 = *(dgm2.begin() + j);
                        sum += pow(norm(p1, p2),  p);
                    }
                }
                break;
            }

    return sum;
}
