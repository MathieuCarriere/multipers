/* Implementations */

template<class T>
void 
ConsistencyList<T>::
swap(iterator i, iterator j)
{
	Parent::swap(i.get_base(),j.get_base());
}

template<class T>
template<class BinaryPredicate>
void 
ConsistencyList<T>::
sort(BinaryPredicate cmp)
{
	Parent::sort(cmp);
	OrderType cur_order = 0;
	for (typename Parent::iterator cur = begin(); cur != end(); ++cur)
		cur->consistency = cur_order++;
}


template<class T>
typename ConsistencyList<T>::iterator 
ConsistencyList<T>::
push_back(const_reference x)
{
	++sz;
	return Parent::push_back(NodeType(x, sz));
}

template<class T>
typename ConsistencyList<T>::iterator 
ConsistencyList<T>::
insert(iterator p, const_reference x)
{
	++sz;
	return Parent::insert(p.get_base(), NodeType(x, sz));
}


/* OrderComparison */
template<class T>
int 
ConsistencyList<T>::
OrderComparison::compare(ComparableType a, ComparableType b) const
{
	typename Parent::OrderComparison cmp;
	return cmp.compare(a.get_base(), b.get_base());
}


/* LessThanComparison */
template<class T>
int ConsistencyList<T>::LessThanComparison::compare(ComparableType a, ComparableType b) const
{ return Parent::compare(a,b); }

template<class T>
bool ConsistencyList<T>::LessThanComparison::operator()(ComparableType a, ComparableType b) const
{ return compare(a,b) == -1; }


/* GreaterThanComparison */
template<class T>
int ConsistencyList<T>::GreaterThanComparison::compare(ComparableType a, ComparableType b) const
{ return -Parent::compare(a,b); }

template<class T>
bool ConsistencyList<T>::GreaterThanComparison::operator()(ComparableType a, ComparableType b) const
{ return compare(a,b) == -1; }


/* ConsistencyComparison */
template<class T>
int ConsistencyList<T>::ConsistencyComparison::compare(ComparableType a, ComparableType b) const
{ 
	if (a.get_base()->consistency < b.get_base()->consistency) 			return -1;
	else if (a.get_base()->consistency == b.get_base()->consistency)	return 0;
	else																return 1;
}

template<class T>
bool ConsistencyList<T>::ConsistencyComparison::operator()(ComparableType a, ComparableType b) const
{ return compare(a,b) == -1; }

