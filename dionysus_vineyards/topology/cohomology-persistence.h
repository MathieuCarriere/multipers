#ifndef __COHOMOLOGY_PERSISTENCE_H__
#define __COHOMOLOGY_PERSISTENCE_H__

#if DEBUG_CONTAINERS
    #include <debug/list>
    #include <debug/vector>
    namespace s = std::__debug;
    #warning "Using debug/list and debug/vector in CohomologyPersistence"
#else
    #include <list>
    #include <vector>
    namespace s = std;
#endif

#include <vector>
#include <list>
#include <utility>

#include <topology/field-arithmetic.h>
#include "utilities/types.h"

#include <boost/optional.hpp>
#include <boost/intrusive/list.hpp>
namespace bi = boost::intrusive;

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>


template<class BirthInfo_, class SimplexData_ = Empty<>, class Field_ = ZpField>
class CohomologyPersistence
{
    public:
        typedef             BirthInfo_                                                  BirthInfo;
        typedef             SimplexData_                                                SimplexData;
        typedef             Field_                                                      Field;

        typedef             typename Field::Element                                     FieldElement;


                            CohomologyPersistence(const Field& field = Field()):
                                field_(field), image_begin_(cocycles_.end())            {}


        // An entry in a cocycle column
        struct  SNode;      // members: si, coefficient, ci
        typedef             s::vector<SNode>                                            ZColumn;
        typedef             bi::list<SNode, bi::constant_time_size<false> >             ZRow;
        class   CompareSNode;

        struct  SHead;      // members: row, order
        typedef             s::list<SHead>                                              Simplices;
        typedef             typename Simplices::iterator                                SimplexIndex;

        struct  Cocycle;    // members: zcolumn, birth, order
        typedef             s::list<Cocycle>                                            Cocycles;
        typedef             typename Cocycles::iterator                                 CocycleIndex;
        typedef             std::pair<CocycleIndex, FieldElement>                       CocycleCoefficientPair;

        typedef             boost::optional<BirthInfo>                                  Death;
        typedef             boost::shared_ptr<ZColumn>                                  CocyclePtr;
        typedef             boost::tuple<SimplexIndex, Death, CocyclePtr>               IndexDeathCocycle;

        // return either a SimplexIndex or a Death
        // BI = BoundaryIterator; it should dereference to a SimplexIndex
        template<class BI>
        IndexDeathCocycle   add(BI begin, BI end, BirthInfo b, bool store = true, const SimplexData& sd = SimplexData(), bool image = true);
        
        // if sign needs to be specified explicitly, provide (parallel) coefficient_iter
        template<class BI, class CI>
        IndexDeathCocycle   add(CI coefficient_iter, BI begin, BI end, BirthInfo b, bool store = true, const SimplexData& sd = SimplexData(), bool image = true);

        void                show_cocycles() const;
        CocycleIndex        begin()                                                     { return image_begin_; }
        CocycleIndex        end()                                                       { return cocycles_.end(); }

    private:
        void                add_cocycle(CocycleCoefficientPair& z1, CocycleCoefficientPair& z2);

    private:
        Simplices           simplices_;
        Cocycles            cocycles_;
        CocycleIndex        image_begin_;
        Field               field_;
};
        
// Simplex representation
template<class BirthInfo_, class SimplexData_, class Field_>
struct CohomologyPersistence<BirthInfo_, SimplexData_, Field_>::SHead: public SimplexData
{
                    SHead(const SHead& other):
                        SimplexData(other), order(other.order)                  {}  // don't copy row since we can't
                    SHead(const SimplexData& sd, unsigned o): 
                        SimplexData(sd), order(o)                               {}

    // intrusive list corresponding to row of s in Z^*, not ordered in any particular order
    ZRow            row;
    unsigned        order;
};

// An entry in a cocycle column; it's also an element in an intrusive list, hence the list_base_hook<>
typedef             bi::list_base_hook<bi::link_mode<bi::auto_unlink> >         auto_unlink_hook;
template<class BirthInfo_, class SimplexData_, class Field_>
struct CohomologyPersistence<BirthInfo_, SimplexData_, Field_>::SNode: public auto_unlink_hook
{
                    SNode(const SNode& other):
                        si(other.si), coefficient(other.coefficient), 
                        ci(other.ci)                                            {}

                    SNode(SimplexIndex sidx, FieldElement coef, CocycleIndex cidx): 
                        si(sidx), coefficient(coef), ci(cidx)                   {}

    SimplexIndex    si;
    FieldElement    coefficient;

    CocycleIndex    ci;                         // TODO: is there no way to get rid of this overhead?

    void            unlink()                    { auto_unlink_hook::unlink(); }
};

template<class BirthInfo_, class SimplexData_, class Field_>
struct CohomologyPersistence<BirthInfo_, SimplexData_, Field_>::Cocycle
{
                    Cocycle(const BirthInfo& b, unsigned o):
                        birth(b), order(o)                                      {}

    ZColumn         zcolumn;
    BirthInfo       birth;
    signed          order;

    bool            operator<(const Cocycle& other) const                       { return order > other.order; }
    bool            operator==(const Cocycle& other) const                      { return order == other.order; }
};


#include "cohomology-persistence.hpp"

#endif // __COHOMOLOGY_PERSISTENCE_H__
