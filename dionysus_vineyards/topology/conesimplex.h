/**
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2007
 */

#ifndef __CONESIMPLEX_H__
#define __CONESIMPLEX_H__

#include <list>
#include <iostream>

#include "utilities/types.h"
#include "simplex.h"


template<class V, class T = Empty>
class ConeSimplex: public Simplex<V,T>
{
    public:
        typedef     Simplex<V,T>                                        Parent;
        typedef     ConeSimplex<V,T>                                    Self;

    public:
                                ConeSimplex(const Self& s): 
                                    Parent(s), coned_(s.coned_)         {}
                                ConeSimplex(const Parent& parent, 
                                            bool coned = false):
                                    Parent(parent), coned_(coned)       {}
        
        Cycle                   boundary() const;
        bool                    coned() const                           { return coned_; }
        Dimension               dimension() const                       { return coned_ ? (Parent::dimension() + 1) : Parent::dimension(); }
        
        bool                    operator==(const Self& other) const     { return !(coned_ ^ other.coned_) && Parent::operator==(other); }

        std::ostream&           operator<<(std::ostream& out) const;

        struct ConedVertexComparison;
        
    private:
        bool                    coned_;
};

template<class V, class T>
struct ConeSimplex<V,T>::ConedVertexComparison: public typename Simplex<V,T>::VertexComparison
{
        typedef     typename Simplex<V,T>::VertexComparison         Parent; 
    
        bool                    operator()(const Self& a, const Self& b) const       
        { 
            if (a.coned() ^ b.coned())
                return b.coned();                   // coned simplices shall come after non-coned ones
            else
                return Parent::operator()(a,b);     // within coned/non-coned order by vertices
        }
};

template<class V, class T>
std::ostream&       operator<<(std::ostream& out, const ConeSimplex<V,T>& s)  { return s.operator<<(out); }

#include "conesimplex.hpp"

#endif // __CONESIMPLEX_H__
