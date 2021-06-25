#include <utilities/log.h>
#include <utilities/boost.h>
#include <boost/iterator/filter_iterator.hpp>
#include <algorithm>
#include <utilities/indirect.h>
#include <functional>

#ifdef LOGGING
static rlog::RLogChannel* rlZigzagAdd =                   DEF_CHANNEL("topology/persistence/zigzag/add",        rlog::Log_Debug);
static rlog::RLogChannel* rlZigzagRemove =                DEF_CHANNEL("topology/persistence/zigzag/remove",     rlog::Log_Debug);
static rlog::RLogChannel* rlZigzagAddChain =              DEF_CHANNEL("topology/persistence/zigzag/addchain",   rlog::Log_Debug);
static rlog::RLogChannel* rlZigzagCheckConsistency=       DEF_CHANNEL("topology/persistence/zigzag/check",      rlog::Log_Debug);
#endif // LOGGING

#ifdef COUNTERS
static Counter*  cZigzagAdd =                             GetCounter("zigzag/add");
static Counter*  cZigzagRemove =                          GetCounter("zigzag/remove");
static Counter*  cZigzagConsistency =                     GetCounter("zigzag/consistency");
#endif // COUNTERS

template<class BID, class SD>
template<class Visitor>
typename ZigzagPersistence<BID,SD>::IndexDeathPair
ZigzagPersistence<BID,SD>::
add(ZColumn bdry, const BirthID& birth, Visitor& visitor)
{
    Count(cZigzagAdd);

    rLog(rlZigzagAdd,       "Entered ZigzagPersistence::add()");
    rLog(rlZigzagAdd,       "  Boundary: %s", bdry.tostring(out).c_str());
    rLog(rlZigzagAdd,       "  Boundary size: %d", bdry.size());
    AssertMsg(check_consistency(), "Must be consistent before addition");

    SimplexIndex last_s     = visitor.new_simplex(*this);
    last_s->low             = z_list.end();
#if ZIGZAG_CONSISTENCY
    last_s->boundary        = bdry;     // NB: debug only
#endif

    rLog(rlZigzagAdd,   "  Reducing among cycles");
    // Reduce bdry among the cycles
    rLog(rlZigzagAdd,       "    Boundary: %s", bdry.tostring(out).c_str());
    BColumn v;                // representation of the boundary in the cycle basis
    while (!bdry.empty())
    {
        SimplexIndex l      = bdry.back();
        ZIndex k            = l->low;
        v.append(k, cmp);
        bdry.add(k->z_column, cmp);
        rLog(rlZigzagAdd,       "    Boundary: %s", bdry.tostring(out).c_str());
    }
    rLog(rlZigzagAdd,   "  Reduced among cycles");

    // Reduce v among boundaries
    BRow u;
    while (!(v.empty()))
    {
        ZIndex l = v.back();
        BIndex k = l->low;

        if (k == b_list.end())
            break;

        u.append(k, cmp);
        v.add(k->b_column, cmp);
    }
    rLog(rlZigzagAdd,   "  Reduced among boundaries");

    if (v.empty())
    {
        rLog(rlZigzagAdd,       "  Birth case in add");

        // Figure out the new cycle z
        ZColumn z;
        std::for_each(u.begin(), u.end(), make_adder(&BNode::c_column, z));
        z.append(last_s, cmp);

        // Birth
        ZIndex new_z                = visitor.new_z_in_add(*this, z, u);
        new_z->birth                = birth;

        // Set s_row
        std::for_each(z.begin(), z.end(), make_appender(&SimplexNode::z_row, new_z));

        // Set z_column
        new_z->z_column.swap(z);

        // Set low
        new_z->low                  = b_list.end();
        last_s->low                 = new_z;

        return std::make_pair(last_s, Death());
    } else
    {
        rLog(rlZigzagAdd,       "  Death case in add");

        // Death
        unsigned order              = b_list.empty() ? 0 : boost::prior(b_list.end())->order + 1;
        b_list.push_back(BNode(order));
        BIndex last_b               = boost::prior(b_list.end());

        // Set b_column and low
        last_b->b_column.swap(v);
        last_b->b_column.back()->low = last_b;

        // Set b_row
        std::for_each(last_b->b_column.begin(), last_b->b_column.end(), make_appender(&ZNode::b_row, last_b));

        // Set c_column
        CColumn& c                  = last_b->c_column;
        std::for_each(u.begin(), u.end(), make_adder(&BNode::c_column, c));
        c.append(last_s, cmp);

        // Set c_row
        std::for_each(c.begin(), c.end(), make_appender(&SimplexNode::c_row, last_b));

        return std::make_pair(last_s, visitor.death(*this, last_b->b_column.back()));
    }
}


template<class BID, class SD>
template<class Visitor>
typename ZigzagPersistence<BID,SD>::Death
ZigzagPersistence<BID,SD>::
remove(SimplexIndex s, const BirthID& birth, Visitor& visitor)
{
    Count(cZigzagRemove);

    rLog(rlZigzagRemove,        "Entered ZigzagPersistence::remove(%d)", s->order);
    AssertMsg(check_consistency(), "Must be consistent before removal");

    if (s->z_row.empty())
    {
        AssertMsg(!(s->c_row.empty()),  "Birth after removal, row in C must be non-empty");

        // Birth
        //show_all();
        rLog(rlZigzagRemove,        "Birth case in remove");

        // Prepend DC[j] = ZB[j] to Z
        rLog(rlZigzagRemove,        "Computing the column DC[j] = ZB[j] to prepend to Z");
        BIndex j                    = visitor.select_j_in_remove(*this, s->c_row);
        rLog(rlZigzagRemove,        "  j = %d", j->order);

        ZColumn z;
        std::for_each(j->b_column.begin(), j->b_column.end(), make_adder(&ZNode::z_column, z));

        ZIndex first_z              = visitor.new_z_in_remove(*this);
        first_z->birth              = birth;
        std::for_each(z.begin(),           z.end(),           make_appender(&SimplexNode::z_row, first_z));
        first_z->z_column.swap(z);
        first_z->low                = b_list.end();

        rLog(rlZigzagRemove,        "  Prepended %d [%s]", first_z->order, z.tostring(out).c_str());
        //AssertMsg(check_consistency(),  "Must be consistent after prepending DC[j] = ZB[j] to Z");

        // Prepend row of s in C to B
        rLog(rlZigzagRemove,        "Prepending the row of s to B");
        first_z->b_row = s->c_row;      // copying instead of swapping is inefficient,
                                        // but it simplifies logic when subtracting chains later
        std::for_each(first_z->b_row.begin(), first_z->b_row.end(), make_appender(&BNode::b_column, first_z));
        //AssertMsg(check_consistency(),  "Must be consistent after prepending row of s to B");

#if ZIGZAG_CONSISTENCY
        {
            ZColumn zz;
            std::for_each(j->b_column.begin(), j->b_column.end(), make_adder(&ZNode::z_column, zz));
            AssertMsg(zz.empty(),       "ZB[j] must be 0 after we prepended the row of s in C to B");
        }
#endif

        typedef         std::not_equal_to<BIndex>       NotEqualBIndex;

        // Subtract C[j] from every column of C that contains s
        AssertMsg(s->c_row == first_z->b_row,   "s->c_row == first_z->b_row before subtracting C[j]");
        rLog(rlZigzagRemove,        "Subtracting C[j]=[%s] from every column of C that contains s=%d with row [%s]",
                                    j->c_column.tostring(out).c_str(),
                                    s->order, s->c_row.tostring(out).c_str());
        add_chains(boost::make_filter_iterator(std::bind2nd(NotEqualBIndex(), j), first_z->b_row.begin(), first_z->b_row.end()),
                   boost::make_filter_iterator(std::bind2nd(NotEqualBIndex(), j), first_z->b_row.end(),   first_z->b_row.end()),
                   j, &BNode::c_column, &SimplexNode::c_row);
        add_chain(j, j, &BNode::c_column, &SimplexNode::c_row);
        // TODO: that's how it was done before, now it can be removed
        //       add_chains(first_z->b_row.rbegin(), first_z->b_row.rend(), j, &BNode::c_column, &SimplexNode::c_row);
        //AssertMsg(check_consistency(s_list.end(), z_list.begin(), b_list.end()),  "Must be consistent after subtracting C[j] in remove::birth");

        // Subtract B[j] from every other column of B that has l
        ZIndex l                    = j->b_column.back();
        BRow   l_row                = l->b_row;
        rLog(rlZigzagRemove,    "Subtracting B[j], j is %d, l is %d, l_row: [%s]",
                                j->order, l->order, l_row.tostring(out).c_str());
        add_chains(boost::make_filter_iterator(std::bind2nd(NotEqualBIndex(), j), l_row.rbegin(), l_row.rend()),
                   boost::make_filter_iterator(std::bind2nd(NotEqualBIndex(), j), l_row.rend(),   l_row.rend()),
                   j, &BNode::b_column, &ZNode::b_row);
        j->b_column.back()->low = b_list.end();     // redundant since l will be deleted (here for check_consistency only)
        add_chain(j, j, &BNode::b_column, &ZNode::b_row);
        // TODO: investigate why this works for ordinary zigzag, but fails for the image zigzag
        //AssertMsg(check_consistency(s_list.end(), first_z, b_list.end()),  "Must be consistent after subtracting B[j] in remove::birth");


        // Drop j, l, and s
        //
        // l->z_column is the only non-empty thing, but we drop it,
        // the basis is preserved because we added first_z
        l->z_column.back()->low     = z_list.end();
        std::for_each(l->z_column.begin(), l->z_column.end(), make_remover(&SimplexNode::z_row, l));

        //show_all();
        rLog(rlZigzagRemove,        "l=%d has z_column: [%s]", l->order, l->z_column.tostring(out).c_str());

        AssertMsg(l->b_row.empty(),     "b_row of l must be empty before erasing in remove::birth");
        AssertMsg(s->z_row.empty(),     "z_row of s must be empty before erasing in remove::birth");
        rLog(rlZigzagRemove,            "s->c_row: [%s]", s->c_row.tostring(out).c_str());
        if (!s->c_row.empty())
        {
            rLog(rlZigzagRemove,        "s->c_row[0]: [%s]", s->c_row.front()->c_column.tostring(out).c_str());
            rLog(rlZigzagRemove,        "   b_column: [%s]", s->c_row.front()->b_column.tostring(out).c_str());
        }
        AssertMsg(s->c_row.empty(),     "c_row of s must be empty before erasing in remove::birth");
        visitor.erasing_z(*this, l);
        b_list.erase(j);
        z_list.erase(l);
        s_list.erase(s);
        AssertMsg(check_consistency(s_list.end(), first_z, b_list.end()),  "Must be consistent before reducing Z in remove::birth");

        // Reduce Z
        rLog(rlZigzagRemove,        "Reducing Z");
        SimplexIndex ls = first_z->z_column.back();
        while(ls->low != first_z)
        {
            if (ls->low == z_list.end())    { ls->low = first_z; break; }

            // if ls->low precedes first_z, swap them
            if (cmp(ls->low, first_z))      std::swap(ls->low, first_z);

            add_chain(first_z, ls->low, &ZNode::b_row, &BNode::b_column);
            add_chain(ls->low, first_z, &ZNode::z_column, &SimplexNode::z_row);
            std::swap(ls->low, first_z);

            ls = first_z->z_column.back();
        }
        AssertMsg(check_consistency(),  "Must be consistent at the end of birth case in remove");

        return Death();
    } else
    {
        // Death
        rLog(rlZigzagRemove,        "Death case in remove");

        ZIndex j                    = s->z_row.front();
        CRow c_row                  = s->c_row;

        rLog(rlZigzagRemove,        "j=%d, j->b_row=[%s]", j->order, j->b_row.tostring(out).c_str());
        rLog(rlZigzagRemove,        "s=%d, s->c_row=[%s]", s->order, s->c_row.tostring(out).c_str());
        rLog(rlZigzagRemove,        "s=%d, s->z_row=[%s]", s->order, s->z_row.tostring(out).c_str());

        // Subtract Z[j] from every chain in C that contains s
        // (it's Ok to go in the forward order since we are subtracting a column in Z from C)
        add_chains(c_row.begin(), c_row.end(), j, &BNode::c_column, &SimplexNode::c_row, &ZNode::z_column);
        AssertMsg(check_consistency(),  "Must be consistent after subtracting Z[j] from C");

        // Change basis to remove s from Z
        // Compute reducers --- columns that we will be adding to other columns
        ZRow z_row                  = s->z_row;
        typedef typename ZRow::reverse_iterator             ZRowReverseIterator;
        typedef std::list<ZRowReverseIterator>              ReducersContainer;
        ReducersContainer  reducers;                        // list of ZColumns that we will be adding to other columns
        reducers.push_back(boost::prior(z_row.rend()));     // j is the first reducer
        AssertMsg(*(reducers.back()) == j,  "The first element of reducers should be j");
        SimplexIndex low            = j->z_column.back();
        rLog(rlZigzagRemove,        "   Added reducer %d [%s] with low=%d",
                                    j->order, j->z_column.tostring(out).c_str(),
                                    low->order);
        for (typename ZRow::iterator cur = z_row.begin(); cur != z_row.end(); ++cur)
            if (cmp((*cur)->z_column.back(), low))
            {
                reducers.push_back(ZRowReverseIterator(boost::next(cur)));
                low = (*(reducers.back()))->z_column.back();
                rLog(rlZigzagRemove,        "   Added reducer %d [%s] with low=%d",
                                            (*cur)->order, (*cur)->z_column.tostring(out).c_str(),
                                            low->order);
                rLog(rlZigzagRemove,        "   reducers.back(): %d [%s] with low=%d",
                                            (*(reducers.back()))->order,
                                            (*(reducers.back()))->z_column.tostring(out).c_str(),
                                            (*(reducers.back()))->z_column.back()->order);
            }
        rLog(rlZigzagRemove,        " Reducers size: %d, s is %d",
                                    reducers.size(), s->order);

        //show_all();

        // Add each reducer to the columns that follow them until the next reducer
        // NB: processing reducers in the reverse order fixes a bug in the paper,
        //     in step Remove.Death.1.4, where the matrix B is updated incorrectly.
        //     I can't find a mention of this bug in my notes, but the fact
        //     that it's fixed in the code suggests that I knew about it. (Or
        //     most likely I didn't recognize that what the paper said is not
        //     exactly what I meant.)
        typename ReducersContainer::reverse_iterator    cur     = reducers.rbegin();
        ZRowReverseIterator                             zcur    = z_row.rbegin();

        while (cur != reducers.rend())
        {
            rLog(rlZigzagRemove,        " Cur reducer: %d [%s]", (**cur)->order,
                                                                 (**cur)->z_column.tostring(out).c_str());
            change_basis(zcur, *cur, **cur,
                         &ZNode::z_column, &SimplexNode::z_row,
                         &ZNode::b_row,    &BNode::b_column);
            if (cur != reducers.rbegin())
            {
                AssertMsg((*zcur)->z_column.back() == (**cur)->z_column.back(),
                          "The back of the z_columns must be the same.");
                (*zcur)->z_column.back()->low = *zcur;
            }
            else
                (**cur)->z_column.back()->low = z_list.end();

            zcur = *cur++;
            // This makes it inconsistent until the next iteration of this update loop
        }

        // Drop j and s
        Death d                     = visitor.death(*this, j);

        if (j->z_column.back()->low == j)
            j->z_column.back()->low = z_list.end();
        std::for_each(j->z_column.begin(), j->z_column.end(), make_remover(&SimplexNode::z_row, j));
        rLog(rlZigzagRemove,            "j->b_row: [%s]", j->b_row.tostring(out).c_str());
        if (!j->b_row.empty())
        {
            rLog(rlZigzagRemove,        "j->b_row[0]: [%s]", j->b_row.front()->b_column.tostring(out).c_str());
            rLog(rlZigzagRemove,        "   c_column: [%s]", j->b_row.front()->c_column.tostring(out).c_str());
        }
        AssertMsg(j->b_row.empty(),     "b_row of j must be empty before erasing in remove(). Most likely you are trying to remove a simplex whose coface is still in the complex.");
        AssertMsg(s->z_row.empty(),     "z_row of s must be empty before erasing in remove()");
        AssertMsg(s->c_row.empty(),     "c_row of s must be empty before erasing in remove()");
        visitor.erasing_z(*this, j);
        z_list.erase(j);
        s_list.erase(s);

        //show_all();

        AssertMsg(check_consistency(),  "Must be consistent when done in remove()");

        return d;
    }
}

template<class BID, class SD>
void
ZigzagPersistence<BID,SD>::
show_all()
{
    std::cout << "s_list:" << std::endl;
    for (SimplexIndex cur = s_list.begin(); cur != s_list.end(); ++cur)
    {
        std::cout << "  " << cur->order << ":" << std::endl;

        std::cout << "    z_row: ";
        for (typename ZRow::const_iterator zcur = cur->z_row.begin(); zcur != cur->z_row.end(); ++zcur)
            std::cout << (*zcur)->order << " ";
        std::cout << std::endl;

        std::cout << "    c_row: ";
        for (typename CRow::const_iterator ccur = cur->c_row.begin(); ccur != cur->c_row.end(); ++ccur)
            std::cout << (*ccur)->order << " ";
        std::cout << std::endl;

        std::cout << "    low: ";
        if (cur->low != z_list.end())
            std::cout << cur->low->order;
        else
            std::cout << "none";
        std::cout << std::endl;
    }

    std::cout << "z_list:" << std::endl;
    for (ZIndex cur = z_list.begin(); cur != z_list.end(); ++cur)
    {
        std::cout << "  " << cur->order << ":" << std::endl;

        std::cout << "    birth: " << cur->birth << std::endl;

        std::cout << "    z_column: ";
        for (typename ZColumn::const_iterator zcur = cur->z_column.begin(); zcur != cur->z_column.end(); ++zcur)
            std::cout << (*zcur)->order << " ";
        std::cout << std::endl;

        std::cout << "    b_row: ";
        for (typename BRow::const_iterator bcur = cur->b_row.begin(); bcur != cur->b_row.end(); ++bcur)
            std::cout << (*bcur)->order << " ";
        std::cout << std::endl;

        std::cout << "    low: ";
        if (cur->low != b_list.end())
            std::cout << cur->low->order;
        else
            std::cout << "none";
        std::cout << std::endl;
    }

    std::cout << "b_list:" << std::endl;
    for (BIndex cur = b_list.begin(); cur != b_list.end(); ++cur)
    {
        std::cout << "  " << cur->order << ":" << std::endl;

        std::cout << "    b_column: ";
        for (typename BColumn::const_iterator bcur = cur->b_column.begin(); bcur != cur->b_column.end(); ++bcur)
            std::cout << (*bcur)->order << " ";
        std::cout << std::endl;

        std::cout << "    c_column: ";
        for (typename CColumn::const_iterator ccur = cur->c_column.begin(); ccur != cur->c_column.end(); ++ccur)
            std::cout << (*ccur)->order << " ";
        std::cout << std::endl;
    }
}

template<class BID, class SD>
bool
ZigzagPersistence<BID,SD>::
check_consistency(SimplexIndex, ZIndex, BIndex)
{
#ifdef ZIGZAG_CONSISTENCY
    #warning "Checking consistency in ZigzagPersistence"

    Count(cZigzagConsistency);
    for (SimplexIndex cur = s_list.begin(); cur != s_list.end(); ++cur)
    {
        if (cur == s_skip) continue;
        //rLog(rlZigzagCheckConsistency,      "SimplexIndex cur: %d", cur->order);
        for (typename ZRow::const_iterator zcur = cur->z_row.begin(); zcur != cur->z_row.end(); ++zcur)
            if (std::find((*zcur)->z_column.begin(), (*zcur)->z_column.end(), cur) == (*zcur)->z_column.end())
            {
                rError("In check_consistency(): SimplexNode %d not found in z_column of %d", cur->order, (*zcur)->order);
                return false;
            }
        for (typename CRow::const_iterator ccur = cur->c_row.begin(); ccur != cur->c_row.end(); ++ccur)
            if (std::find((*ccur)->c_column.begin(), (*ccur)->c_column.end(), cur) == (*ccur)->c_column.end())
            {
                rError("In check_consistency(): SimplexNode %d not found in c_column of %d", cur->order, (*ccur)->order);
                return false;
            }
        if (cur->low != z_list.end())
            AssertMsg(!(cur->low->z_column.empty()),        "z_column must not be empty");
        if (cur->low != z_list.end() && cur->low->z_column.back() != cur)
        {
            rError("low of SimplexNode %d is incorrect", cur->order);
            return false;
        }
    }

    for (ZIndex cur = z_list.begin(); cur != z_list.end(); ++cur)
    {
        if (cur == z_skip) continue;

        //rLog(rlZigzagCheckConsistency,      "ZIndex cur: %d", cur->order);
        for (typename ZColumn::const_iterator scur = cur->z_column.begin(); scur != cur->z_column.end(); ++scur)
            if (std::find((*scur)->z_row.begin(), (*scur)->z_row.end(), cur) == (*scur)->z_row.end())
            {
                rError("In check_consistency(): ZNode %d not found in z_row of %d", cur->order, (*scur)->order);
                return false;
            }
        for (typename BRow::const_iterator bcur = cur->b_row.begin(); bcur != cur->b_row.end(); ++bcur)
            if (std::find((*bcur)->b_column.begin(), (*bcur)->b_column.end(), cur) == (*bcur)->b_column.end())
            {
                rError("In check_consistency(): ZNode %d not found in b_column of %d", cur->order, (*bcur)->order);
                return false;
            }
        if (cur->low != b_list.end() && cur->low->b_column.back() != cur)
        {
            rError("low of ZNode %d is incorrect", cur->order);
            return false;
        }
        if (cur->z_column.back()->low != cur)
        {
            rError("The low of the back of the z_column must be set correctly");
            rError("  %d [%s], its back %d with low=%d", cur->order,
                                                         cur->z_column.tostring(out).c_str(),
                                                         cur->z_column.back()->order,
                                                         (cur->z_column.back()->low == z_list.end()) ? 0 : cur->z_column.back()->low->order);
            return false;
        }
    }

    for (BIndex cur = b_list.begin(); cur != b_list.end(); ++cur)
    {
        if (cur == b_skip) continue;

        //rLog(rlZigzagCheckConsistency,      "BIndex cur: %d", cur->order);
        for (typename BColumn::const_iterator zcur = cur->b_column.begin(); zcur != cur->b_column.end(); ++zcur)
            if (std::find((*zcur)->b_row.begin(), (*zcur)->b_row.end(), cur) == (*zcur)->b_row.end())
            {
                rError("In check_consistency(): BNode %d not found in b_row of %d", cur->order, (*zcur)->order);
                return false;
            }
        for (typename CColumn::const_iterator scur = cur->c_column.begin(); scur != cur->c_column.end(); ++scur)
            if (std::find((*scur)->c_row.begin(), (*scur)->c_row.end(), cur) == (*scur)->c_row.end())
            {
                rError("In check_consistency(): BNode %d not found in c_row of %d", cur->order, (*scur)->order);
                return false;
            }
        if (!(cur->b_column.empty() || cur->b_column.back()->low == cur))
        {
            rError("The low of the back of the b_column must be set correctly");
            return false;
        }

        // ZB == DC
        ZColumn zb, dc;
        std::for_each(cur->b_column.begin(), cur->b_column.end(), make_adder(&ZNode::z_column, zb));
        std::for_each(cur->c_column.begin(), cur->c_column.end(), make_adder(&SimplexNode::boundary, dc));
        zb.add(dc, cmp);
        if (!zb.empty())
        {
            rError("   b_column: [%s]",    cur->b_column.tostring(out).c_str());
            rError("   c_column: [%s]",    cur->c_column.tostring(out).c_str());
            rError("   zb - dc:  [%s]",    zb.tostring(out).c_str());
            rError("ZB = DC");
            return false;
        }
    }
#endif

    return true;
}

/* Private */

// Class: Appender
//
// Functor that appends given element to the given member of whatever parameter it is invoked with
template<class BID, class SD>
template<class Member, class Element>
struct ZigzagPersistence<BID,SD>::Appender
{
                Appender(Member mm, Element ee):
                    m(mm), e(ee)                        {}

    template<class T>
    void        operator()(T& a)                        { ((*a).*m).append(e, cmp); }

    Member          m;
    Element         e;
    OrderComparison cmp;
};

// Class: Remover
//
// Functor that removes given element from the given member of whatever parameter it is invoked with
template<class BID, class SD>
template<class Member, class Element>
struct ZigzagPersistence<BID,SD>::Remover
{
                Remover(Member mm, Element ee):
                    m(mm), e(ee)                        {}

    template<class T>
    void        operator()(T& a)                        { ((*a).*m).remove(e); }

    Member  m;
    Element e;
};

// Class: Adder
//
// Functor that adds the given member of whatever it is invoked with to the given chain
template<class BID, class SD>
template<class Member, class Chain>
struct ZigzagPersistence<BID,SD>::Adder
{
                Adder(Member mm, Chain& cc):
                    m(mm), c(cc)                        {}

    template<class T>
    void        operator()(T& a)                        { c.add((*a).*m, cmp); }

    Member          m;
    Chain&          c;
    OrderComparison cmp;
};


// Function: add_chains()
//
// Special case of add_chains where all Indexes are the same, and
// therefore PrimaryMemberFrom and PrimaryMemberTo are the same
template<class BID, class SD>
template<class Index, class IndexFrom, class PrimaryMember, class SecondaryMember>
void
ZigzagPersistence<BID,SD>::
add_chains(Index bg, Index end, IndexFrom j, PrimaryMember pm, SecondaryMember sm)
{
    add_chains(bg, end, j, pm, sm, pm);
}

// Function: add_chains()
//
// Adds PrimaryMember pm of j to pm of every element in the range [bg,end)
// Fixes SecondaryMembers by adding and removing the corresponding elements.
// For example, if we add a column to a number of other columns, then PrimaryMember is that
// column member, and SecondaryMember is the corresponding row member.
template<class BID, class SD>
template<class IndexTo, class IndexFrom, class PrimaryMemberTo, class SecondaryMemberTo, class PrimaryMemberFrom>
void
ZigzagPersistence<BID,SD>::
add_chains(IndexTo bg, IndexTo end, IndexFrom j, PrimaryMemberTo pmt, SecondaryMemberTo smt, PrimaryMemberFrom pmf)
{
    for (IndexTo cur = bg; cur != end; ++cur)
        add_chain(*cur, j, pmt, smt, pmf);
}

// Function: change_basis()
//
// Changes basis by adding PrimaryMember pm of j to pm of every element in range [bg, end).
// In parallel it performs the reverse (complementary) update on the dual members, i.e.
// column and row operations are performed in sync, so that the product of the two matrices doesn't change
template<class BID, class SD>
template<class IndexTo, class IndexFrom, class PrimaryMember, class SecondaryMember, class DualPrimaryMember, class DualSecondaryMember>
void
ZigzagPersistence<BID,SD>::
change_basis(IndexTo bg, IndexTo end, IndexFrom j, PrimaryMember pm, SecondaryMember sm, DualPrimaryMember dpm, DualSecondaryMember dsm)
{
    for (IndexTo cur = bg; cur != end; ++cur)
    {
        add_chain(*cur, j,  pm,  sm,  pm);
        add_chain(j, *cur, dpm, dsm, dpm);
    }
}

template<class BID, class SD>
template<class Index, class PrimaryMember, class SecondaryMember>
void
ZigzagPersistence<BID,SD>::
add_chain(Index to, Index from, PrimaryMember pm, SecondaryMember sm)
{
    add_chain(to, from, pm, sm, pm);
}

// Function: add_chain()
//
// Adds PrimaryMemberFrom pmf of `from` to PrimaryMemberTo pmt of `to`.
// Fixes SecondaryMemberTos. See add_chains().
template<class BID, class SD>
template<class IndexTo, class IndexFrom, class PrimaryMemberTo, class SecondaryMemberTo, class PrimaryMemberFrom>
void
ZigzagPersistence<BID,SD>::
add_chain(IndexTo to, IndexFrom from, PrimaryMemberTo pmt, SecondaryMemberTo smt, PrimaryMemberFrom pmf)
{
    rLog(rlZigzagAddChain,  "Adding %d [%s] to %d [%s]",
                            (*from).order,
                            ((*from).*pmf).tostring(out).c_str(),
                            (*to).order,
                            ((*to).*pmt).tostring(out).c_str());

    // Fix secondaries
    std::for_each(make_intersection_iterator(((*from).*pmf).begin(),  ((*from).*pmf).end(),
                                               ((*to).*pmt).begin(),    ((*to).*pmt).end(),
                                             cmp),
                  make_intersection_iterator(((*from).*pmf).end(),    ((*from).*pmf).end(),
                                               ((*to).*pmt).end(),      ((*to).*pmt).end(),
                                             cmp),
                  make_remover(smt, to));
    std::for_each(make_difference_iterator(((*from).*pmf).begin(),    ((*from).*pmf).end(),
                                             ((*to).*pmt).begin(),      ((*to).*pmt).end(),
                                           cmp),
                  make_difference_iterator(((*from).*pmf).end(),      ((*from).*pmf).end(),
                                             ((*to).*pmt).end(),        ((*to).*pmt).end(),
                                           cmp),
                  make_appender(smt, to));

    // Add primaries
    ((*to).*pmt).add((*from).*pmf, cmp);
    rLog(rlZigzagAddChain,  "Got %s", ((*to).*pmt).tostring(out).c_str());
}


/* ZigzagVisitor */
template<class BID, class SD>
typename ZigzagPersistence<BID,SD>::SimplexIndex
ZigzagPersistence<BID,SD>::ZigzagVisitor::
new_simplex(ZigzagPersistence& zz)
{
    unsigned order      = zz.s_list.empty() ? 0 : boost::prior(zz.s_list.end())->order + 1;
    zz.s_list.push_back(SimplexNode(order, zz.z_list.end()));
    return boost::prior(zz.s_list.end());
}

template<class BID, class SD>
typename ZigzagPersistence<BID,SD>::ZIndex
ZigzagPersistence<BID,SD>::ZigzagVisitor::
new_z_in_add(ZigzagPersistence& zz, const ZColumn&, const BRow&)
{
    int order                   = zz.z_list.empty() ? 0 : boost::prior(zz.z_list.end())->order + 1;
    zz.z_list.push_back(ZNode(order, zz.b_list.end()));
    return boost::prior(zz.z_list.end());
}

template<class BID, class SD>
typename ZigzagPersistence<BID,SD>::BIndex
ZigzagPersistence<BID,SD>::ZigzagVisitor::
select_j_in_remove(ZigzagPersistence&, const CRow& c_row)
{
    return c_row.front();
}

template<class BID, class SD>
typename ZigzagPersistence<BID,SD>::ZIndex
ZigzagPersistence<BID,SD>::ZigzagVisitor::
new_z_in_remove(ZigzagPersistence& zz)
{
    int order                   = zz.z_list.empty() ? 0 : zz.z_list.begin()->order - 1;
    zz.z_list.push_front(ZNode(order, zz.b_list.end()));
    return zz.z_list.begin();
}

template<class BID, class SD>
typename ZigzagPersistence<BID,SD>::Death
ZigzagPersistence<BID,SD>::ZigzagVisitor::
death(ZigzagPersistence&, ZIndex dying_z)
{
    return Death(dying_z->birth);
}
