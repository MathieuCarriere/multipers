/*
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2005 -- 2009
 */

#ifndef __VINEYARD_H__
#define __VINEYARD_H__

#include "utilities/types.h"
#include <list>
#include <string>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/list.hpp>
    
#include <boost/iterator/iterator_traits.hpp>


class Knee;
class Vine;

/**
 * Vineyard class. Keeps track of vines and knees. switched() is the key function called
 * when pairing switches.
 *
 * \ingroup topology
 */
template<class Index_, class Iterator_, class Evaluator_>
class Vineyard
{
    public:
        typedef                         Index_                                          Index;
        typedef                         Iterator_                                       Iterator;
        typedef                         Evaluator_                                      Evaluator;

        typedef                         std::list<Vine>                                 VineList;
        typedef                         std::list<VineList>                             VineListList;
        typedef                         std::vector<VineListList::iterator>             VineListVector;
                                        
    public:
                                        Vineyard(Evaluator* eval = 0): 
                                            evaluator(eval)                             {}

        void                            start_vines(Iterator bg, Iterator end);
        void                            switched(Index i, Index j);
        template<class Iter>
        void                            record_knee(Iter i);
        void                            record_diagram(Iterator bg, Iterator end);

        void                            set_evaluator(Evaluator* eval)                  { evaluator = eval; }

        void                            save_edges(const std::string& filename, bool skip_infinite = false) const;
        void                            save_vines(const std::string& filename, bool skip_infinite = false) const;
        std::vector<std::vector<std::vector<double>>>             get_vines(const int& discard) const;
        std::vector<std::vector<std::vector<double>>>             get_dgms(const int& discard, const int& num) const;

    private:
        template<class Iter>
        void                            start_vine(Iter i);

    private:
        VineListList                    vines;            // stores vine lists
        VineListVector                  vines_vector;     // stores pointers (iterators) to vine lists
        Evaluator*                      evaluator;
};

/**
 * Knee class stores the knee in R^3.
 *
 * \ingroup topology
 */
class Knee
{
    public:
        RealType                birth;
        RealType                death;
        RealType                time;
            
                                // Default parameters for serialization
                                Knee(RealType b = 0, RealType d = 0, RealType t = 0):
                                    birth(b), death(d), time(t)
                                {}
                                Knee(const Knee& other): 
                                    birth(other.birth), death(other.death), time(other.time)
                                {}

        bool                    is_diagonal() const                             { return birth == death; }
        bool                    is_infinite() const                             { return (death == Infinity) || (birth == Infinity); }

        std::ostream&           operator<<(std::ostream& out) const             { return out << "(" << birth << ", " 
                                                                                                    << death << ", " 
                                                                                                    << time  << ")"; }
    
    private:
        friend class boost::serialization::access;

        template<class Archive>
        void                    serialize(Archive& ar, version_type );
};

std::ostream& operator<<(std::ostream& out, const Knee& k)                      { return k.operator<<(out); }

/**
 * Vine is a list of Knees
 */
class Vine: public std::list<Knee>
{   
    public:
        typedef                 std::list<Knee>                                 VineRepresentation;
        typedef                 VineRepresentation::const_iterator              const_knee_iterator;
        
                                Vine()                                          {}
                                Vine(const Vine& other): 
                                    VineRepresentation(other)                   {}
                                Vine(const VineRepresentation& other): 
                                    VineRepresentation(other)                   {}
                                Vine(const Knee& k)                             { add(k); }
        
        void                    add(RealType b, RealType d, RealType t)         { push_back(Knee(b,d,t)); }
        void                    add(const Knee& k)                              { push_back(k); }

        std::ostream&           operator<<(std::ostream& out) const             { std::copy(begin(), end(), std::ostream_iterator<Knee>(out, " ")); return out; }

        using VineRepresentation::begin;
        using VineRepresentation::end;
        using VineRepresentation::front;
        using VineRepresentation::back;
        using VineRepresentation::size;
        using VineRepresentation::empty;

    protected:
        using VineRepresentation::push_back;

    private:
        friend class boost::serialization::access;

        template<class Archive>
        void                    serialize(Archive& ar, version_type );
};

std::ostream& operator<<(std::ostream& out, const Vine& v)                      { return v.operator<<(out); }


class VineData
{
    public:
        void        set_vine(Vine* vine) const                                          { vine_ = vine; }
        Vine*       vine() const                                                        { return vine_; }

    private:
        mutable Vine*       vine_;      // cheap trick to work around MultiIndex's constness
};


#include "vineyard.hpp"

#endif // __VINEYARD_H__
