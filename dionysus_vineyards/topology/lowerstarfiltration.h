/*
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2005 -- 2008
 */

#ifndef __LOWERSTARFILTRATION_H__
#define __LOWERSTARFILTRATION_H__

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/nvp.hpp>

/**
 * Struct: MaxVertexComparison
 *
 * Functor that determines which simplex has a higher vertex with respect to VertexComparison_, breaking ties by dimension
 */
template<class Simplex_, class VertexComparison_>
struct MaxVertexComparison
{
    typedef                     VertexComparison_                                   VertexComparison;
    typedef                     Simplex_                                            Simplex;
    typedef                     typename Simplex::Vertex                            Vertex;

                                MaxVertexComparison(const VertexComparison& vcmp):
                                    vcmp_(vcmp)                                     {}

    bool                        operator()(const Simplex& s1, const Simplex& s2) const
    {
        const Vertex& max1 = *std::max_element(s1.vertices().begin(), s1.vertices().end(), vcmp_);
        const Vertex& max2 = *std::max_element(s2.vertices().begin(), s2.vertices().end(), vcmp_);
        
        bool less = vcmp_(max1, max2), 
             more = vcmp_(max2, max1);
        
        if (!less && !more)     // equal
            return s1.dimension() < s2.dimension();

        return less;
    }

    VertexComparison            vcmp_;
};


/**
 * Map from i-th vertex to its index in the filtration.
 */
template<class Index_, class Filtration_>
class VertexSimplexMap
{
    public:
        typedef                 Index_                                              Index;
        typedef                 Filtration_                                         Filtration;
        typedef                 std::vector<Index>                                  VertexVector;
                                
                                VertexSimplexMap(Index begin, Index end, const Filtration& f)
        {
            for (Index cur = begin; cur != end; ++cur)
                if (f.simplex(cur).dimension() == 0)
                    vertices_.push_back(cur);
        }

    private:
        VertexVector            vertices_;
};


#endif // __LOWERSTARFILTRATION_H__
