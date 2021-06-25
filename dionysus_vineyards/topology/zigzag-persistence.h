#ifndef __ZIGZAG_PERSISTENCE_H__
#define __ZIGZAG_PERSISTENCE_H__

#include "cycles.h"
#include "utilities/types.h"
#include <sstream>

#if DEBUG_CONTAINERS
    #include <debug/list>
    using std::__debug::list;
    #warning "Using debug/list in ZigzagPersistence"
#else
    #include <list>
    using std::list;
#endif


/**
 * Class: ZigzagPersistence
 * TODO: this should probably be parametrized by Chain or Field
 */
template<class BirthID_ = Empty<>, class SimplexData_ = Empty<> >
class ZigzagPersistence
{
    public:
        typedef                         BirthID_                                BirthID;
        typedef                         SimplexData_                            SimplexData;

        struct ZNode;
        struct BNode;
        struct SimplexNode;

        typedef                         list<ZNode>                                     ZList;
        typedef                         typename ZList::iterator                        ZIndex;
        typedef                         list<BNode>                                     BList;
        typedef                         typename BList::iterator                        BIndex;
        typedef                         list<SimplexNode>                               SimplexList;
        typedef                         typename SimplexList::iterator                  SimplexIndex;

        // TODO: should all chains be DequeChains? probably not
        typedef                         typename DequeChains<ZIndex>::Chain             ZRow;
        typedef                         typename DequeChains<ZIndex>::Chain             BColumn;
        typedef                         typename VectorChains<BIndex>::Chain            BRow;
        typedef                         typename VectorChains<BIndex>::Chain            CRow;
        typedef                         typename VectorChains<SimplexIndex>::Chain      ZColumn;
        typedef                         typename VectorChains<SimplexIndex>::Chain      CColumn;

        typedef                         boost::optional<BirthID>                        Death;
        typedef                         std::pair<SimplexIndex, Death>                  IndexDeathPair;

        struct ZNode
        {
                                        ZNode(int o, BIndex l):
                                            order(o), low(l)                            {}

            int                         order;
            ZColumn                     z_column;
            BRow                        b_row;
            BIndex                      low;            // which BColumn has this ZIndex as low

            BirthID                     birth;          // TODO: do we need to do empty-member optimization?
                                                        //       i.e., does it ever make sense for birth to be empty?
        };

        struct BNode
        {
                                        BNode(unsigned o): order(o)                     {}

            unsigned                    order;
            BColumn                     b_column;
            CColumn                     c_column;
        };

        struct SimplexNode: public SimplexData
        {
                                        SimplexNode(unsigned o, ZIndex l):
                                            order(o), low(l)                            {}

            unsigned                    order;
            ZRow                        z_row;
            CRow                        c_row;
            ZIndex                      low;            // which ZColumn has this SimplexNode as low
#if ZIGZAG_CONSISTENCY
            ZColumn                     boundary;       // NB: debug only
#endif
        };

        // Constructor: ZigzagPersistence()
                                        ZigzagPersistence()                             {}

        // Function: add(bdry, birth)
        IndexDeathPair                  add(ZColumn bdry, const BirthID& birth = BirthID())         { ZigzagVisitor zzv; return add<ZigzagVisitor>(bdry, birth, zzv); }

        // Function: remove(s, birth)
        Death                           remove(SimplexIndex s, const BirthID& birth = BirthID())    { ZigzagVisitor zzv; return remove<ZigzagVisitor>(s, birth, zzv); }


        ZIndex                          begin()                                                     { return z_list.begin(); }
        ZIndex                          end()                                                       { return z_list.end(); }

        bool                            is_alive(ZIndex i) const                                    { return i->low == b_list.end(); }
        bool                            is_alive(ZNode zn) const                                    { return zn.low == b_list.end(); }

    protected:
        // Function: add(s)
        template<class Visitor>
        IndexDeathPair                  add(ZColumn bdry, const BirthID& birth, Visitor& visitor);

        // Function: remove(s)
        template<class Visitor>
        Death                           remove(SimplexIndex s, const BirthID& birth, Visitor& visitor);

        // Struct: ZigzagVisitor
        // Various methods of an instance of this class are called at different stages of addition and removal algorithm.
        // NB: currently the places where it's called are catered for image zigzags, in the future this could be expanded
        //     to provide simple support for other algorithms
        // TODO: not obvious that the methods should be const (and therefore the reference passed to add() and remove())
        //       revisit when working on ImageZigzag
        struct ZigzagVisitor
        {
            SimplexIndex                new_simplex(ZigzagPersistence& zz);

            // Function: new_z_in_add(zz, z, u)
            // Called when a new cycle is born after adding a simplex. The method is expected to add an element to z_list, and return its ZIndex.
            ZIndex                      new_z_in_add(ZigzagPersistence& zz, const ZColumn& z, const BRow& u);

            BIndex                      select_j_in_remove(ZigzagPersistence& zz, const CRow& c_row);

            ZIndex                      new_z_in_remove(ZigzagPersistence& zz);

            void                        erasing_z(ZigzagPersistence&, ZIndex)           {}

            Death                       death(ZigzagPersistence& zz, ZIndex dying_z);
        };

    public:
        // Debug; non-const because Indices are iterators, and not const_iterators
        void                            show_all();
        bool                            check_consistency(SimplexIndex s_skip, ZIndex z_skip, BIndex b_skip);
        bool                            check_consistency()                             { return check_consistency(s_list.end(), z_list.end(), b_list.end()); }

    protected:
        ZList                           z_list;
        BList                           b_list;
        SimplexList                     s_list;

        /* Helper functors */
        template<class Member, class Element>                                           struct Appender;
        template<class Member, class Element>                                           struct Remover;
        template<class Member, class Chain>                                             struct Adder;

        template<class Member, class Element>
        Appender<Member, Element>       make_appender(Member m, Element e) const        { return Appender<Member, Element>(m,e); }
        template<class Member, class Element>
        Remover<Member, Element>        make_remover(Member m, Element e) const         { return Remover<Member, Element>(m,e); }
        template<class Member, class Chain>
        Adder<Member, Chain>            make_adder(Member m, Chain& c) const            { return Adder<Member, Chain>(m, c); }

        template<class Index, class IndexFrom, class PrimaryMember, class SecondaryMember>
        void                            add_chains(Index bg, Index end, IndexFrom j, PrimaryMember pm, SecondaryMember sm);
        template<class IndexTo, class IndexFrom, class PrimaryMemberTo, class SecondaryMemberTo, class PrimaryMemberFrom>
        void                            add_chains(IndexTo bg, IndexTo end, IndexFrom j,
                                                   PrimaryMemberTo   pmt, SecondaryMemberTo smt,
                                                   PrimaryMemberFrom pmf);
        template<class Index, class PrimaryMember, class SecondaryMember>
        void                            add_chain(Index to, Index from,
                                                  PrimaryMember   pmt, SecondaryMember smt);
        template<class IndexTo, class IndexFrom, class PrimaryMemberTo, class SecondaryMemberTo, class PrimaryMemberFrom>
        void                            add_chain(IndexTo to, IndexFrom from,
                                                  PrimaryMemberTo   pmt, SecondaryMemberTo smt,
                                                  PrimaryMemberFrom pmf);
        template<class IndexTo, class IndexFrom, class PrimaryMember, class SecondaryMember, class DualPrimaryMember, class DualSecondaryMember>
        void                            change_basis(IndexTo bg, IndexTo end, IndexFrom j,
                                                     PrimaryMember pm, SecondaryMember sm,
                                                     DualPrimaryMember dpm, DualSecondaryMember dsm);

    public:
        struct OrderComparison
        {
            template<class T>
            bool                        operator()(T a, T b) const                      { return a->order < b->order; }
        }                               cmp;

        struct OrderOutput
        {
            template<class T>
            std::string                 operator()(T a) const                           { std::stringstream s; s << a->order; return s.str(); }
        }                               out;
};

#include "zigzag-persistence.hpp"

#endif // __ZIGZAG_PERSISTENCE_H__
