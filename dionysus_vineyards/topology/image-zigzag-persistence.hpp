#ifdef LOGGING
static rlog::RLogChannel* rlImageZigzag =                   DEF_CHANNEL("topology/persistence/zigzag/image",        rlog::Log_Debug);
#endif // LOGGING


template<class BID>
typename ImageZigzagPersistence<BID>::SimplexIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
new_simplex(Parent& zz)
{
    SimplexIndex si = Parent::ZigzagVisitor::new_simplex(zz);
    si->subcomplex = subcomplex;
    rLog(rlImageZigzag,     "New simplex %d (inL=%d)", si->order, si->subcomplex);
    return si;
}

template<class BID>
typename ImageZigzagPersistence<BID>::ZIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
new_z_in_add(Parent& zz, const ZColumn& z, const BRow& u)
{ 
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 

    // Check if z is entirely in the subcomplex
    if (in_subcomplex(z))
        return append_in_image(zz);
    else
    {
        append_in_cokernel(zz);
        
        // Simplex we are adding is in the subcomplex
        if (subcomplex)
        {
            rLog(rlImageZigzag,     "Modifying boundaries");

            SimplexIndex s      = boost::prior(izz.s_list.end());
            AssertMsg(s->subcomplex,        "The new simplex must be in the subcomplex");
            rLog(rlImageZigzag,     "  s=%d", s->order);
            rLog(rlImageZigzag,     "  u=[%s]", u.tostring(izz.out).c_str());

            BIndex max = u.front();
            for (typename BRow::const_iterator cur = u.begin(); cur != u.end(); ++cur)
                if ((*cur)->b_column.back()->order > max->b_column.back()->order)     
                    max = *cur;
            
            // Replace B[max] with sum of b_columns from u
            BColumn sum;
            std::for_each(u.begin(), u.end(), izz.make_adder(&BNode::b_column, sum));
#if ZIGZAG_CONSISTENCY
            AssertMsg(in_subcomplex(s->boundary), "Boundary of s must be in the subcomplex");
#endif
            std::for_each(max->b_column.begin(), max->b_column.end(), izz.make_remover(&ZNode::b_row, max));
            rLog(rlImageZigzag,     "max->b_column=[%s]", max->b_column.tostring(izz.out).c_str());
            rLog(rlImageZigzag,     "          sum=[%s]", sum.tostring(izz.out).c_str());
            AssertMsg(sum.back() == max->b_column.back(), "Must be replacing by a column with the same low");
            max->b_column.swap(sum);
            std::for_each(max->b_column.begin(), max->b_column.end(), izz.make_appender(&ZNode::b_row, max));
            // NB: low doesn't need to be adjusted (see AssertMsg above)

            // Wipe out C[max], and replace it with s
            std::for_each(max->c_column.begin(), max->c_column.end(), izz.make_remover(&SimplexNode::c_row, max));
            max->c_column.clear();
            max->c_column.append(s, izz.cmp);
            AssertMsg(s->c_row.empty(),     "s was just added, so it cannot appear in any bounding chain");
            s->c_row.append(max, izz.cmp);
        }
        
        return boost::prior(izz.z_list.end());
    }
}


template<class BID>
typename ImageZigzagPersistence<BID>::BIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
select_j_in_remove(Parent& zz, const CRow& c_row)
{
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 
    for (typename CRow::const_iterator cur = c_row.begin(); cur != c_row.end(); ++cur)
        if ((*cur)->b_column.back()->order <= izz.im_last->order)
        {
            birth_in_image = true;
            return *cur;
        }

    return Parent::ZigzagVisitor::select_j_in_remove(zz, c_row);
}

template<class BID>
typename ImageZigzagPersistence<BID>::ZIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
new_z_in_remove(Parent& zz)
{ 
    if (birth_in_image)
        return prepend_in_image(zz);
    else
        return prepend_in_cokernel(zz);
}

template<class BID>
void
ImageZigzagPersistence<BID>::ImageZZVisitor::
erasing_z(Parent& zz, ZIndex j)
{ 
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 
    if          (j == izz.im_last)          --(izz.im_last);
    else if     (j == izz.cok_begin)        ++(izz.cok_begin);
}

template<class BID>
typename ImageZigzagPersistence<BID>::Death
ImageZigzagPersistence<BID>::ImageZZVisitor::
death(Parent& zz, ZIndex dying_z)
{
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 
    if (izz.im_last == izz.z_list.end() || dying_z->order > izz.im_last->order)
        return Death();
    else
        return Death(dying_z->birth);
}

template<class BID>
typename ImageZigzagPersistence<BID>::ZIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
append_in_image(Parent& zz)
{
    rLog(rlImageZigzag,     "Appending in image");
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 

    // if no cycles in the image
    if (izz.im_last == izz.z_list.end())
    {
        izz.z_list.push_front(ZNode(izz.im_order_begin, izz.b_list.end()));
        return (izz.im_last = izz.z_list.begin());
    } else
    {
        izz.z_list.insert(boost::next(izz.im_last), ZNode(izz.im_last->order + 1, izz.b_list.end()));
        return ++(izz.im_last);
    }
}

template<class BID>
typename ImageZigzagPersistence<BID>::ZIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
append_in_cokernel(Parent& zz)
{
    rLog(rlImageZigzag,     "Appending in cokernel");
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 

    // if no cycles in the cokernel
    if (izz.cok_begin == izz.z_list.end())
    {
        izz.z_list.push_back(ZNode(izz.cok_order_begin, izz.b_list.end()));
        izz.cok_begin = boost::prior(izz.z_list.end());
    } else
    {
        izz.z_list.push_back(ZNode(boost::prior(izz.z_list.end())->order + 1, izz.b_list.end()));
    }

    return boost::prior(izz.z_list.end());
}

template<class BID>
typename ImageZigzagPersistence<BID>::ZIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
prepend_in_image(Parent& zz)
{
    rLog(rlImageZigzag,     "Prepending in image");
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 

    // if no cycles in the image
    if (izz.im_last == izz.z_list.end())
    {
        izz.z_list.push_front(ZNode(izz.im_order_begin, izz.b_list.end()));
        return (izz.im_last = izz.z_list.begin());
    } else
    {
        izz.z_list.push_front(ZNode(izz.z_list.begin()->order - 1, izz.b_list.end()));
        return izz.z_list.begin();
    }
}

template<class BID>
typename ImageZigzagPersistence<BID>::ZIndex
ImageZigzagPersistence<BID>::ImageZZVisitor::
prepend_in_cokernel(Parent& zz)
{
    rLog(rlImageZigzag,     "Prepending in cokernel");
    ImageZigzagPersistence& izz = static_cast<ImageZigzagPersistence&>(zz); 

    // if no cycles in the cokernel
    if (izz.cok_begin == izz.z_list.end())
    {
        izz.z_list.push_back(ZNode(izz.cok_order_begin, izz.b_list.end()));
        izz.cok_begin = boost::prior(izz.z_list.end());
    } else
    {
        izz.z_list.insert(izz.cok_begin, ZNode(izz.cok_begin->order - 1, izz.b_list.end()));
        --(izz.cok_begin);
    }
    return izz.cok_begin;
}

template<class BID>
bool
ImageZigzagPersistence<BID>::ImageZZVisitor::
in_subcomplex(const ZColumn& z)
{
    for (typename ZColumn::const_iterator cur = z.begin(); cur != z.end(); ++cur)
        if (!((*cur)->subcomplex))
            return false;
    return true;
}
