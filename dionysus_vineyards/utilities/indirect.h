#ifndef __INDIRECT_H__
#define __INDIRECT_H__

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <iterator>

// TODO: write documentation

/**************
 * Comparison *
 **************/

template<class Comparison_>
struct IndirectComparison
{
    typedef         Comparison_                             Comparison;

                    IndirectComparison(const Comparison& cmp):
                        cmp_(cmp)
                    {}

    template<class Iterator>
    bool            operator()(Iterator a, Iterator b) const
    { return cmp_(*a, *b); }

    const Comparison&               cmp_;
};
template<class Comparison>
IndirectComparison<Comparison> make_indirect_comparison(const Comparison& cmp)
{ return IndirectComparison<Comparison>(cmp); }


template<class Comparison_>
class FirstComparison
{
    public:
        typedef     Comparison_                             Comparison;

                    FirstComparison(Comparison cmp): 
                        cmp_(cmp)                           {}

        template<class Pair>
        bool        operator()(const Pair& a, const Pair& b) const
        { return cmp_(a.first, b.first); }

    private:
        Comparison  cmp_;
};
template<class Comparison>
FirstComparison<Comparison> make_first_comparison(const Comparison& cmp)
{ return FirstComparison<Comparison>(cmp); }


template<class Comparison>
struct ThreeOutcomeCompare: public Comparison
{
    typedef         typename Comparison::first_argument_type                            first_argument_type;
    typedef         typename Comparison::second_argument_type                           second_argument_type;

                    ThreeOutcomeCompare(const Comparison& cmp = Comparison()): Comparison(cmp)              {}

    int             compare(first_argument_type a, second_argument_type b) const
    {   if (this->operator()(a,b))          return -1;
        else if (this->operator()(b,a))     return 1;
        else                                return 0;
    }
};

template<class Evaluator_>
class ThroughEvaluatorComparison
{
    public:
        typedef                 Evaluator_                                                  Evaluator;

                                ThroughEvaluatorComparison(const Evaluator& eval): 
                                    eval_(eval)                                             {}

        template<class T>                            
        bool                    operator()(T a, T b) const                                  { return (eval_(a) < eval_(b)); }

    private:
        const Evaluator&        eval_;
};


/*************
 * Iterators *
 *************/

// Iterates over the difference of the two sorted sequences, dereferencing into the first sequence
template<class Iterator1, class Iterator2, class StrictWeakOrdering>
class difference_iterator: public boost::iterator_facade<difference_iterator<Iterator1, Iterator2, StrictWeakOrdering>,
                                                         typename std::iterator_traits<Iterator1>::value_type,
                                                         boost::forward_traversal_tag>
{
    public:
        typedef     typename std::iterator_traits<Iterator1>::reference                 reference;

                    difference_iterator(Iterator1 cur1, Iterator1 end1,
                                        Iterator2 cur2, Iterator2 end2,
                                        const StrictWeakOrdering& cmp = StrictWeakOrdering()):
                        cur1_(cur1), end1_(end1), cur2_(cur2), end2_(end2), cmp_(cmp)   { catchup(); }

    private:
        friend class boost::iterator_core_access;

        void        increment()                                                         { ++cur1_; catchup(); }
        bool        equal(const difference_iterator& other) const                       { return (cur1_ == other.cur1_) && (cur2_ == other.cur2_); }
        reference   dereference() const                                                 { return *cur1_; }

    private:
        Iterator1   cur1_, end1_;
        Iterator2   cur2_, end2_;
        const StrictWeakOrdering&   cmp_;

    private:
        void        catchup()
        {
            while ((cur1_ != end1_) && (cur2_ != end2_))
            {
                if      (cmp_(*cur1_, *cur2_))      break;
                else if (cmp_(*cur2_, *cur1_))      ++cur2_;
                else                                { ++cur1_; ++cur2_; }
            }

            if (cur1_ == end1_)                     cur2_ = end2_;
        }
};

template<class Iterator1, class Iterator2, class StrictWeakOrdering>
difference_iterator<Iterator1, Iterator2, StrictWeakOrdering>
make_difference_iterator(Iterator1 cur1, Iterator1 end1, Iterator2 cur2, Iterator2 end2, const StrictWeakOrdering& cmp)
{ return difference_iterator<Iterator1, Iterator2, StrictWeakOrdering>(cur1, end1, cur2, end2, cmp); }

// Iterates over the intersection of the two sorted sequences, dereferencing into the first sequence
template<class Iterator1, class Iterator2, class StrictWeakOrdering>
class intersection_iterator: public boost::iterator_facade<intersection_iterator<Iterator1, Iterator2, StrictWeakOrdering>,
                                                           typename std::iterator_traits<Iterator1>::value_type,
                                                           boost::forward_traversal_tag>
{
    public:
        typedef     typename std::iterator_traits<Iterator1>::reference                 reference;

                    intersection_iterator(Iterator1 cur1, Iterator1 end1,
                                          Iterator2 cur2, Iterator2 end2,
                                          const StrictWeakOrdering& cmp = StrictWeakOrdering()):
                        cur1_(cur1), end1_(end1), cur2_(cur2), end2_(end2), cmp_(cmp)   { catchup(); }

    private:
        friend class boost::iterator_core_access;

        void        increment()                                                         { ++cur1_; ++cur2_; catchup(); }
        bool        equal(const intersection_iterator& other) const                     { return (cur1_ == other.cur1_) && (cur2_ == other.cur2_); }
        reference   dereference() const                                                 { return *cur1_; }

    private:
        Iterator1   cur1_, end1_;
        Iterator2   cur2_, end2_;
        const StrictWeakOrdering&   cmp_;

    private:
        void        catchup()
        {
            while ((cur1_ != end1_) && (cur2_ != end2_))
            {
                if      (cmp_(*cur1_, *cur2_))      ++cur1_;
                else if (cmp_(*cur2_, *cur1_))      ++cur2_;
                else                                break;
            }

            if ((cur1_ == end1_) || (cur2_ == end2_))
            {
                cur1_ = end1_;
                cur2_ = end2_;
            }
        }
};

template<class Iterator1, class Iterator2, class StrictWeakOrdering>
intersection_iterator<Iterator1, Iterator2, StrictWeakOrdering>
make_intersection_iterator(Iterator1 cur1, Iterator1 end1, Iterator2 cur2, Iterator2 end2, const StrictWeakOrdering& cmp)
{ return intersection_iterator<Iterator1, Iterator2, StrictWeakOrdering>(cur1, end1, cur2, end2, cmp); }

#endif // __INDIRECT_H__
