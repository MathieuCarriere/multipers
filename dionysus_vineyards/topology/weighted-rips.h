#ifndef __WEIGHTED_RIPS_H__
#define __WEIGHTED_RIPS_H__

#include <vector>
#include <string>
#include "simplex.h"
#include "rips.h"
#include <boost/iterator/counting_iterator.hpp>

/**
 * WeightedRipsSimplex class
 * 
 * This class sits as an invisible layer between the Simplex datatype passed
 * to WeightedRips and the class itself. The need for this layer is the need
 * to store the ``value'' (max inter-vertex distance) of each simplex in the
 * Weighted Rips complex--something that the user of the class does not need
 * to be aware of.
 */

template<class Simplex_, class DistanceType_>
class WeightedRipsSimplex : public Simplex_
{
    public:
        typedef          typename Simplex_::Vertex                 Vertex;
        typedef          typename Simplex_::VertexContainer        VertexContainer;
        typedef          DistanceType_                             DistanceType;

        WeightedRipsSimplex(Simplex_ s) : Simplex_(s)  { }

        void             setSimplexValue(const DistanceType &sv) { simplexValue = sv; }
        DistanceType     getSimplexValue() const       { return    simplexValue;      }

    protected:
        DistanceType     simplexValue;
};

/**
 * WeightedRips class
 *
 * Class providing basic operations to work with Rips complexes. It implements Bron-Kerbosch algorithm, 
 * and provides simple wrappers for various functions.
 *
 * Distances_ is expected to define types IndexType and DistanceType as well as 
 *               provide operator()(...) which given two IndexTypes should return 
 *               the distance between them. There should be methods begin() and end() 
 *               for iterating over IndexTypes as well as a method size().
 */
template<class Distances_, class Simplex_ = Simplex<typename Distances_::IndexType> >
class WeightedRips : public Rips<Distances_, Simplex_>
{
    public:

        /* redeclaring the typedefs because they cannot be inherited at compile-time */
        typedef             Distances_                                      Distances; 
        typedef             typename Distances::IndexType                   IndexType;
        typedef             typename Distances::DistanceType                DistanceType;

        typedef             WeightedRipsSimplex<Simplex_, DistanceType>     Simplex;
        typedef             typename Simplex::Vertex                        Vertex;             // should be the same as IndexType
        typedef             typename Simplex::VertexContainer               VertexContainer;

        class               Evaluator;
        class               Comparison;

    public:
                            WeightedRips(const Distances& distances):
                                Rips<Distances_, Simplex_>(distances)                             {}

                            template<class Functor>
                            void generate(Dimension k, DistanceType max, const Functor& f) const;

};

/**
 * DistanceDataStackingFunctor class
 * 
 * Class providing a functor that is to be called by WeightedRips::generate(). This functor
 * takes as an argument (to its constructor) the original functor passed by the user to
 * generate(), and a new ``double'' functor is created. Assuming that the functor acts on
 * simplices, first the value of the simplex is computed (the radius at which the simplex
 * appears in the weighted Rips complex), the data field of the simplex is populated with
 * this value, and then the original functor is called (it has no idea that it was
 * intercepted).
 */

template<class Rips_, class Functor_>
class DistanceDataStackingFunctor
{
	public:
		typedef     typename Rips_::Simplex         Simplex_;

		DistanceDataStackingFunctor(const Rips_ &r, const Functor_ &f):
			rips(r), original_functor(f) { }

		void operator()(const Simplex_ &s) const
		{
			Simplex_         s_new(s);
			s_new.setSimplexValue (rips.distance(s_new, s_new));
			original_functor      (s_new);
		}

	private:
		const Rips_       &rips;
		const Functor_    &original_functor;
};

template<class Distances_, class Simplex_>
template<class Functor>
void WeightedRips<Distances_, Simplex_>::generate(Dimension k, DistanceType max, const Functor &f) const
{
	Rips<Distances_,Simplex_>::generate(k, max, DistanceDataStackingFunctor<WeightedRips<Distances_, Simplex_>,Functor>(*this, f));
}

template<class Distances_, class Simplex_>
class WeightedRips<Distances_, Simplex_>::Evaluator: public Rips<Distances_,Simplex_>::Evaluator
{
    public:
                            Evaluator(const Distances& distances): 
                                Rips<Distances_, Simplex_>::Evaluator(distances)                       {}

        DistanceType       operator()(const Simplex& s) const { return s.getSimplexValue(); }
};

template<class Distances_, class Simplex_>
class WeightedRips<Distances_, Simplex_>::Comparison: public Rips<Distances_,Simplex_>::Comparison
{
    public:
                            Comparison(const Distances& distances):
                                Rips<Distances_, Simplex_>::Comparison(distances)                            {}

        bool                operator()(const Simplex& s1, const Simplex& s2) const    
        { 
                            if (s1.dimension() != s2.dimension())
                                return s1.dimension() < s2.dimension();
                            return s1.getSimplexValue() < s2.getSimplexValue();
        }
};

#endif // __WEIGHTED_RIPS_H__
