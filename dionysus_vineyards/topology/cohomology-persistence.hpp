#include <utilities/boost.h>
#include <queue>
#include <vector>
#include <limits>

#include <utilities/log.h>
#include <utilities/indirect.h>
#include <utilities/counter.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/foreach.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
namespace bl = boost::lambda;

#include <boost/make_shared.hpp>

#ifdef LOGGING
static rlog::RLogChannel* rlCohomology =                DEF_CHANNEL("topology/cohomology",        rlog::Log_Debug);
#endif

#ifdef COUNTERS
static Counter*  cCohomologyAddBasic =                  GetCounter("cohomology/add/basic");
static Counter*  cCohomologyAddComparison =             GetCounter("cohomology/add/comparison");
static Counter*  cCohomologyElementCount =              GetCounter("cohomology/elements");
static Counter*  cCohomologyCocycleCount =              GetCounter("cohomology/cocycles");
static Counter*  cCohomologyCandidatesCount =           GetCounter("cohomology/candidates");
#endif // COUNTERS

template<class BirthInfo, class SimplexData, class Field>
class CohomologyPersistence<BirthInfo, SimplexData, Field>::CompareSNode
{
    public:
        bool        operator()(const SNode& s1, const SNode& s2) const                  { return s1.si->order < s2.si->order; }
};
    

struct Alternator
{
    typedef         int         result_type;
    int             operator()(int x) const                                             { return (x % 2)*2 - 1; }
};
    

template<class BirthInfo, class SimplexData, class Field>
template<class BI>
typename CohomologyPersistence<BirthInfo, SimplexData, Field>::IndexDeathCocycle
CohomologyPersistence<BirthInfo, SimplexData, Field>::
add(BI begin, BI end, BirthInfo birth, bool store, const SimplexData& sd, bool image)
{
    // Set coefficient to be an iterator over (-1)^i
    
    return add(boost::make_transform_iterator(boost::make_counting_iterator(1), Alternator()),
               begin, end, 
               birth, store, sd, image);
}


template<class BirthInfo, class SimplexData, class Field>
template<class BI, class CI>
typename CohomologyPersistence<BirthInfo, SimplexData, Field>::IndexDeathCocycle
CohomologyPersistence<BirthInfo, SimplexData, Field>::
add(CI coefficient_iter, BI begin, BI end, BirthInfo birth, bool store, const SimplexData& sd, bool image)
{
    // Create simplex representation
    simplices_.push_back(SHead(sd, simplices_.empty() ? 0 : (simplices_.back().order + 1)));
    SimplexIndex    si = boost::prior(simplices_.end());

    // Find out if there are cocycles that evaluate to non-zero on the new simplex
    typedef         std::list<CocycleCoefficientPair>           Candidates;
    Candidates      candidates, candidates_bulk;
    rLog(rlCohomology, "Boundary");

    for (BI cur = begin; cur != end; ++cur)
    {
        FieldElement coefficient = field_.init(*coefficient_iter++);
        SimplexIndex cursi = *cur;

        rLog(rlCohomology, "  %d %d", cursi->order, coefficient);
        BOOST_FOREACH(const SNode& zcur, std::make_pair(cursi->row.begin(), cursi->row.end()))
            candidates_bulk.push_back(std::make_pair(zcur.ci, field_.mul(coefficient, zcur.coefficient)));
    }

    candidates_bulk.sort(make_first_comparison(make_indirect_comparison(std::less<Cocycle>())));
    CountBy(cCohomologyCandidatesCount, candidates_bulk.size());
    
#if LOGGING    
    rLog(rlCohomology,  "  Candidates bulk");
    for (typename Candidates::iterator cur  = candidates_bulk.begin(); 
                                       cur != candidates_bulk.end(); ++cur)
        rLog(rlCohomology, "    %d %d", cur->first->order, cur->second);
#endif

    // Remove duplicates
    {
        typename Candidates::const_iterator cur = candidates_bulk.begin();
        while (cur != candidates_bulk.end())
        {
            typename Candidates::const_iterator next = cur;
            FieldElement sum = field_.zero();
            while (next != candidates_bulk.end() && next->first == cur->first) 
            { 
                sum = field_.add(sum, next->second);
                ++next; 
            }
    
            rLog(rlCohomology, "  Bulk candidate %d sum %d", cur->first->order, sum);
            if (!field_.is_zero(sum))
                candidates.push_back(CocycleCoefficientPair(cur->first, sum));
    
            cur = next;
        }
    }

    // Birth
    if (candidates.empty())
    {
        // rLog(rlCohomology, "  Birth occurred");
        if (!store)
        {
            simplices_.pop_back();
            boost::shared_ptr<ZColumn> p = boost::make_shared<ZColumn>();
            return boost::make_tuple(simplices_.end(), Death(), p);
        }
        
        signed order = 0;
        if (image)
            if (image_begin_ == cocycles_.end())
                order = std::numeric_limits<signed>::min();
            else
                order = image_begin_->order + 1;
        else
            if (!cocycles_.empty() && cocycles_.front().order >= 0)     // we have something outside the image
                order = cocycles_.front().order + 1;

        CocycleIndex nw;
        if (image)
        {
            image_begin_ = cocycles_.insert(image_begin_, Cocycle(birth, order));
            nw = image_begin_;
        } else
        {
            cocycles_.push_front(Cocycle(birth, order));
            nw = cocycles_.begin();
        }
        
        rLog(rlCohomology,  "Birth: %d", nw->order);

        // set up the cocycle
        ZColumn& cocycle = nw->zcolumn;
        cocycle.push_back(SNode(si, field_.id(), nw));
        si->row.push_back(cocycle.front());
        rLog(rlCohomology,  "  Cocyle: %d", si->order);
        
        Count(cCohomologyElementCount);
        Count(cCohomologyCocycleCount);

        boost::shared_ptr<ZColumn> p = boost::make_shared<ZColumn>();
        return boost::make_tuple(si, Death(), p);
    }

    // Death
    rLog(rlCohomology,  "Death");

#if LOGGING
    // Debug only, output candidates
    rLog(rlCohomology,  "  Candidates");
    for (typename Candidates::iterator cur  = candidates.begin(); cur != candidates.end(); ++cur)
        rLog(rlCohomology, "    %d %d", cur->first->order, cur->second);
#endif

    CocycleCoefficientPair& z   = candidates.front();
    Death d                     = z.first->birth;
    rLog(rlCohomology, "  Order: %d", z.first->order);
    if (z.first->order >= 0)    // if death outside image
        d = Death();            // no death occurs outside the image
    else
        if (z.first == image_begin_)
            ++image_begin_;

    // add z to everything else in candidates
    for (typename Candidates::iterator cur  = boost::next(candidates.begin()); 
                                       cur != candidates.end(); ++cur)
    {
        CountBy(cCohomologyElementCount, -cur->first->zcolumn.size());
        add_cocycle(*cur, z);
        CountBy(cCohomologyElementCount, cur->first->zcolumn.size());
    }

    for (typename ZColumn::iterator cur = z.first->zcolumn.begin(); cur != z.first->zcolumn.end(); ++cur)
        cur->unlink();
    
    CountBy(cCohomologyElementCount, -z.first->zcolumn.size());
    CountBy(cCohomologyCocycleCount, -1);
    boost::shared_ptr<ZColumn> p = boost::make_shared<ZColumn>();
    p->swap(z.first->zcolumn);
    cocycles_.erase(z.first);

    return boost::make_tuple(si, d, p);
}
        
template<class BirthInfo, class SimplexData, class Field>
void
CohomologyPersistence<BirthInfo, SimplexData, Field>::
show_cocycles() const
{
    std::cout << "Cocycles: " << cocycles_.size() << std::endl;
    for (typename Cocycles::const_iterator cur = cocycles_.begin(); cur != cocycles_.end(); ++cur)
    {
        // std::cout << cur->order << " (" << cur->birth << "): ";
        for (typename ZColumn::const_iterator zcur = cur->zcolumn.begin(); zcur != cur->zcolumn.end(); ++zcur)
            std::cout << zcur->coefficient << " * " << zcur->si->order << ", ";
        std::cout << std::endl;
    }
}

template<class BirthInfo, class SimplexData, class Field>
void
CohomologyPersistence<BirthInfo, SimplexData, Field>::
add_cocycle(CocycleCoefficientPair& to, CocycleCoefficientPair& from)
{
    rLog(rlCohomology,  "Adding cocycle %d to %d", from.first->order, to.first->order);

    FieldElement    multiplier = field_.neg(field_.div(to.second, from.second));
    CocycleIndex    ci = to.first;
    CompareSNode    cmp;

    // Insert at the end optimization
    if (cmp(to.first->zcolumn.back(), from.first->zcolumn.front()) && 
        to.first->zcolumn.capacity() >= (to.first->zcolumn.size() + from.first->zcolumn.size()))
    {
        BOOST_FOREACH(const SNode& fs, from.first->zcolumn)
        {
            to.first->zcolumn.push_back(SNode(fs.si, field_.mul(multiplier, fs.coefficient), ci));
            fs.si->row.push_back(to.first->zcolumn.back());
            Count(cCohomologyAddBasic);
        }

        return;
    }

    ZColumn         nw;
    typename ZColumn::iterator tcur = to.first->zcolumn.begin();
    typename ZColumn::iterator fcur = from.first->zcolumn.begin();
    while (tcur != to.first->zcolumn.end() && fcur != from.first->zcolumn.end())
    {
        rLog(rlCohomology, "  %d %d", tcur->si->order, fcur->si->order);
        Count(cCohomologyAddComparison);
        Count(cCohomologyAddBasic);
        if (cmp(*tcur, *fcur))
        {
            nw.push_back(*tcur);
            ++tcur;
        }
        else if (cmp(*fcur, *tcur))
        {
            nw.push_back(SNode(fcur->si, field_.mul(multiplier, fcur->coefficient), ci));
            ++fcur;
        }
        else        // equality
        {
            FieldElement res = field_.mul(multiplier, tcur->coefficient);
            res = field_.add(fcur->coefficient, res);
            if (!field_.is_zero(res))
                nw.push_back(SNode(fcur->si, res, ci));
            ++tcur; ++fcur;
        }
    }
    for (; tcur != to.first->zcolumn.end(); ++tcur)
    {
        rLog(rlCohomology, "  %d", tcur->si->order);
        Count(cCohomologyAddBasic);
        nw.push_back(SNode(*tcur));
    }
    for (; fcur != from.first->zcolumn.end(); ++fcur)
    {
        rLog(rlCohomology, "  %d", fcur->si->order);
        Count(cCohomologyAddBasic);
        nw.push_back(SNode(fcur->si, field_.mul(multiplier, fcur->coefficient), ci));
    }


    for (typename ZColumn::iterator cur = to.first->zcolumn.begin(); cur != to.first->zcolumn.end(); ++cur)
        cur->unlink();

    to.first->zcolumn.swap(nw);

    for (typename ZColumn::iterator cur = to.first->zcolumn.begin(); cur != to.first->zcolumn.end(); ++cur)
        cur->si->row.push_back(*cur);
}
