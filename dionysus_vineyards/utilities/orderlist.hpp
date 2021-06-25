/* Implementations */

#ifdef LOGGING
static rlog::RLogChannel* rlOrderList = 					DEF_CHANNEL("utilities/orderlist", rlog::Log_Debug);
#endif // LOGGING

// This cannot be in the header file to conform to the C++ standard
template<class T>
const float OrderList<T>::density_threshold = 1.2;

template<class T>
void 
OrderList<T>::
swap(iterator i, iterator j)
{
	typename Parent::iterator i_base = i.get_base();
	typename Parent::iterator j_base = j.get_base();
	std::swap(i_base->tag, j_base->tag);

	// Exchange the actual elements in the list --- so that iterators behave as expected
	typename Parent::iterator after_j = boost::next(j_base);	
	Parent::splice(i_base, *this, j_base);
	Parent::splice(after_j, *this, i_base);
}

template<class T>
template<class BinaryPredicate>
void 
OrderList<T>::
sort(BinaryPredicate cmp)
{
	Parent::sort(OrderListNodeComparison<T, BinaryPredicate>(cmp));
	OrderType cur_order = 0;
	for (typename Parent::iterator cur = Parent::begin(); cur != Parent::end(); ++cur)
		cur->tag = cur_order++;
}

template<class T>
typename OrderList<T>::iterator 
OrderList<T>::
push_back(const_reference x)
{
	if (empty()) 
		Parent::push_back(NodeType(x, 0));
	else
		Parent::push_back(NodeType(x, last().get_base()->tag + 1));
	
	return last();
}

template<class T>
typename OrderList<T>::iterator 
OrderList<T>::
insert(iterator p, const_reference x)
{
	typename Parent::iterator p_base = p.get_base();
	OrderType tag = (p_base++)->tag + 1;
	typename Parent::iterator new_base = Parent::insert(p_base, NodeType(x, tag));

	if (p_base->tag != tag)
		return iterator(new_base);

	// Find non-overflowing region
	unsigned int num_elements = 1, maximum = 1, lower = tag, upper = tag, level = 0;
	float inv_density = 1;
	typename Parent::iterator prev = p_base, next = p_base;
	--(--prev); ++next; 		// move prev back twice to skip over the newly inserted element

	do
	{
		lower &= ~(1 << level);
		upper |= (1 << level);
		maximum <<= 1; inv_density *= density_threshold;
		++level;

		while (prev != Parent::end() && prev->tag >= lower) { --prev; ++num_elements; }
		while (next != Parent::end() && next->tag <= upper) { ++next; ++num_elements; }
	} while (inv_density * num_elements >= maximum);
	++num_elements;			// for the extra element inserted

	rLog(rlOrderList, "%i, %i, %i", num_elements, lower, upper);
	rLog(rlOrderList, "prev is at the end: %i", (prev == Parent::end()));
	rLog(rlOrderList, "next is at the end: %i", (next == Parent::end()));
	
	// Reorder
	AssertMsg((upper - lower + 1)/num_elements > 0, "Spacing between new tags must be non-zero");
	for (unsigned int i = 0; i < num_elements; ++i)
	{
		(++prev)->tag = lower + i*((upper - lower + 1)/num_elements);
		rLog(rlOrderList, "%i", prev->tag);
		AssertMsg(prev->tag != 0 || prev == Parent::begin(), "Cannot assign 0 tag except at the beginning of OrderList");
	}

	AssertMsg(++prev == next, "prev + num_elements != next in OrderList::insert()");

	return iterator(new_base);
}

template<class T>
void
OrderList<T>::
show_elements() const
{
	for (const_iterator cur = begin(); cur != end(); ++cur)
		std::cout << *(cur.get_base()) << std::endl;
	std::cout << std::endl;
}

/* OrderComparison */
template<class T>
int 
OrderList<T>::OrderComparison::
compare(ComparableType a, ComparableType b) const
{
	if (a.get_base()->tag == b.get_base()->tag)			return 0;
	if (a.get_base()->tag < b.get_base()->tag)			return -1;
	return 1;
}

template<class T>
bool
OrderList<T>::OrderComparison::
operator()(ComparableType a, ComparableType b) const
{
	return (compare(a,b) < 0);
}
