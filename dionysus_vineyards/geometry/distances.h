#ifndef __DISTANCES_H__
#define __DISTANCES_H__

#include <vector>

/**
 * Class: ExplicitDistances 
 * Stores the pairwise distances of Distances_ instance passed at construction. 
 * It's a protypical Distances template argument for the Rips complex.
 */
template<class Distances_>
class ExplicitDistances
{
    public:
        typedef             Distances_                                      Distances;
        typedef             size_t                                          IndexType;
        typedef             typename Distances::DistanceType                DistanceType;

                            ExplicitDistances(IndexType size):
                                size_(size), 
                                distances_(size*(size + 1)/2 + size)        {}
                            ExplicitDistances(const Distances& distances);

        DistanceType        operator()(IndexType a, IndexType b) const;
        DistanceType&       operator()(IndexType a, IndexType b);

        size_t              size() const                                    { return size_; }
        IndexType           begin() const                                   { return 0; }
        IndexType           end() const                                     { return size(); }

    private:
        std::vector<DistanceType>                   distances_;
        size_t                                      size_;
};


/**
 * Class: PairwiseDistances
 * Given a Container_ of points and a Distance_, it computes distances between elements 
 * in the container (given as instances of Index_ defaulted to unsigned) using the Distance_ functor.
 *
 * Container_ is assumed to be an std::vector. That simplifies a number of things.
 */
template<class Container_, class Distance_, typename Index_ = unsigned>
class PairwiseDistances
{
    public:
        typedef             Container_                                      Container;
        typedef             Distance_                                       Distance;
        typedef             Index_                                          IndexType;
        typedef             typename Distance::result_type                  DistanceType;


                            PairwiseDistances(const Container& container, 
                                              const Distance& distance = Distance()):
                                container_(container), distance_(distance)  {}

        DistanceType        operator()(IndexType a, IndexType b) const      { return distance_(container_[a], container_[b]); }

        size_t              size() const                                    { return container_.size(); }
        IndexType           begin() const                                   { return 0; }
        IndexType           end() const                                     { return size(); }

    private:
        const Container&    container_;
        Distance            distance_;
};

#include "distances.hpp"

#endif // __DISTANCES_H__
