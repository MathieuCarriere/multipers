#ifndef __COMPLEX_TRAITS_H__
#define __COMPLEX_TRAITS_H__

#include <vector>
#include <map>
#include "utilities/property-maps.h"

// Class: ComplexTraits
// Class that extracts various properties of a complex.
//
// Parameter:
//   Complex -      the type whose traits we are extracting
template<class Complex_>
struct ComplexTraits
{
    /* 
     * Typedefs: Associated types
     * Types extractable by the traits class.
     *
     * Complex -            the type of the complex
     * Simplex -            the type of the simplex stored in the complex
     * Index -              the type used to refer to the elements of the complex
     * IndexComparison -    the type used to compare indices
     */
    typedef         Complex_                                                Complex;
    typedef         typename Complex::Simplex                               Simplex;
    typedef         typename Complex::Index                                 Index;
    typedef         std::less<Index>                                        IndexComparison;

    /* Typedefs: Maps
     * Maps to go from Indicies to Simplices and back.
     *
     * SimplexIndexMap -    the map from Simplex to the Index (could be implemented via a binary search, or a map)
     */
    typedef         AssociativeMap<std::map<Simplex, Index, 
                                   typename Simplex::VertexComparison> >    SimplexIndexMap;

    // Function: simplex_index_map() 
    // Initializes an SimplexIndexMap given a Complex.
    static SimplexIndexMap
                    simplex_index_map(const Complex& complex)               { return SimplexIndexMap(complex.begin(),
                                                                                                     complex.end()); } 

    static Index    begin(const Complex& complex)                           { return complex.begin(); }
    static Index    end(const Complex& complex)                             { return complex.end(); }

    static const Simplex&
                    simplex(Index i)                                        { return *i; }
};


// Specializations
template<class Simplex_>
struct ComplexTraits<std::vector<Simplex_> >
{
    typedef         std::vector<Simplex_>                                   Complex;
    typedef         Simplex_                                                Simplex;
    typedef         typename Complex::const_iterator                        Index;
    typedef         std::less<Index>                                        IndexComparison;

    typedef         BinarySearchMap<Simplex, Index,
                                    typename Simplex::VertexComparison>     SimplexIndexMap;

    static SimplexIndexMap
                    simplex_index_map(const Complex& complex)               { return SimplexIndexMap(complex.begin(),
                                                                                                     complex.end()); }
    static SimplexIndexMap
                    simplex_index_map(Index bg, Index end)                  { return SimplexIndexMap(bg, end); }

    static Index    begin(const Complex& complex)                           { return complex.begin(); }
    static Index    end(const Complex& complex)                             { return complex.end(); }

    static const Simplex&
                    simplex(Index i)                                        { return *i; }
};

#endif // __COMPLEX_TRAITS_H__
