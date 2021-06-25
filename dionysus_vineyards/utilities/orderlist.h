/*
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2006
 *
 * Implements the simplified order list data strcutre given in ``Two Simplified
 * Algorithms for Maintaining Order in a List'' by Bender et al.
 *
 * Indirection is not implemented, so the insertion cost is amortized O(log n),
 * while comparison and deletion are O(1).
 */

#ifndef __ORDERLIST_H__
#define __ORDERLIST_H__

#include "log.h"

#include <iterator>
#include <iostream>
#include <list>

#include "types.h"

#include <utilities/boost.h>
#include <boost/iterator/iterator_adaptor.hpp>
//#include <boost/type_traits/is_convertible.hpp>
//#include <boost/utility/enable_if.hpp>


typedef					unsigned int					OrderType;

template<class T> 		struct 	OrderListNode;
template<class T>		class   OrderListIterator;
template<class T>		class   const_OrderListIterator;

/**
 * OrderList stores a list of objects while maintaining their total order under
 * the operations of insert(), swap(), and delete(). push_back() is provided for
 * uniformity, OrderComparison member class carries out comparison of the
 * elements.
 */
template<class T>
class OrderList: public std::list<OrderListNode<T> >
{
	public:
		class 			OrderComparison;

		/// OrderComparison type
		//typedef			OrderComparison								OrderComparison;

		typedef			OrderListNode<T>							NodeType;
		typedef			OrderList<T>								Self;
		typedef			std::list<NodeType >						Parent;

		typedef			T											value_type;
		typedef			T&											reference;
		typedef			const T&									const_reference;
		typedef			OrderListIterator<T>						iterator;
		typedef			const_OrderListIterator<T>					const_iterator;

						OrderList()									{}
						~OrderList() 								{ clear(); }

		/// \name Order operations
		void			swap(iterator i, iterator j);				///< Exchanges the order of simplices pointed to by i and j
		template<class BinaryPredicate>
		void			sort(BinaryPredicate cmp);					///< Sorts the elements in accordance with cmp

		/// \name Container operations
		/// @{
		// Explicit calls instead of using declarations for Doxygen
		iterator		push_back(const_reference x);
		iterator		insert(iterator predecessor, const_reference x);	///< Inserts x immediately after predecessor (has to be a valid iterator)
		void			erase(iterator x)							{ Parent::erase(x.get_base()); }

		void			clear()										{ return Parent::clear(); }
		bool			empty() const								{ return Parent::empty(); }
		SizeType		size()	const								{ return Parent::size(); }
		iterator		begin()										{ return iterator(Parent::begin()); }
		const_iterator	begin() const								{ return const_iterator(Parent::begin()); }
		iterator		end()										{ return iterator(Parent::end()); }
		const_iterator	end() const									{ return const_iterator(Parent::end()); }
		reference		back()										{ return Parent::back(); }
		const_reference	back() const								{ return Parent::back(); }
		void			pop_back()									{ return Parent::pop_back(); }

		iterator		last()										{ return iterator(boost::prior(end())); }
		const_iterator	last() const								{ return const_iterator(boost::prior(end())); }
		/// @}

		/// \name Debugging operations
		/// @{
		void			show_elements() const;
		/// @}

	private:
		static const float density_threshold;
};

/// Basic comparison that LessThan and GreaterThan derive from
template<class T>
class OrderList<T>::OrderComparison
{
	public:
		typedef			typename OrderList<T>::const_iterator		ComparableType;
		int 			compare(ComparableType a, ComparableType b) const;				/// (-1,0,1) = a (precedes, ==, succeeds) b
		bool			operator()(ComparableType a, ComparableType b) const;
};

/// Structure storing auxilliary information requred for each node of OrderList
template<class T>
struct OrderListNode
{
	OrderListNode(const T& d, unsigned int t):
		data(d), tag(t)
	{}

	T 				data;
	OrderType		tag;

	std::ostream& 	operator<<(std::ostream& out) const				{ return out << data << ": " << tag; }
};

template<class T, class BinaryPredicate>
class OrderListNodeComparison
{
	public:
		typedef 		OrderListNode<T>								Node;
		OrderListNodeComparison(BinaryPredicate cmp): cmp_(cmp) 		{}
		bool operator()(const Node& a, const Node& b) const				{ return cmp_(a.data, b.data); }

	private:
		BinaryPredicate cmp_;
};

template<class T>
std::ostream&			operator<<(std::ostream& out, const OrderListNode<T>& n)	{ return n.operator<<(out); }

template<class T>
class OrderListIterator: public boost::iterator_adaptor<OrderListIterator<T>,
						 								typename OrderList<T>::Parent::iterator,
						 								T>
{
	private:
		struct			enabler										{};

	public:
		typedef			typename OrderList<T>::Parent				OrderListParent;
		typedef 		boost::iterator_adaptor<OrderListIterator<T>,
												typename OrderListParent::iterator,
												T>					Parent;
		typedef			typename Parent::reference					reference;
		typedef			typename Parent::base_type					base_type;

						OrderListIterator()							{}
						OrderListIterator(const typename OrderListParent::iterator& iter):
    						OrderListIterator::iterator_adaptor_(iter)
						{}
						OrderListIterator(const OrderListIterator<T>& other):
							OrderListIterator::iterator_adaptor_(other.base())
						{}

	private:
		friend class	boost::iterator_core_access;
		reference		dereference() const							{ return Parent::base_reference()->data; }
		base_type&		get_base()									{ return Parent::base_reference(); }

		friend class 	OrderList<T>;
};

template<class T>
class const_OrderListIterator: public boost::iterator_adaptor<const_OrderListIterator<T>,
						 									  typename OrderList<T>::Parent::const_iterator,
						 									  const T>
{
	private:
		struct			enabler										{};

	public:
		typedef			typename OrderList<T>::Parent				OrderListParent;
		typedef 		boost::iterator_adaptor<const_OrderListIterator<T>,
												typename OrderListParent::const_iterator,
												const T>			Parent;
		typedef			typename Parent::reference					reference;
		typedef			typename Parent::base_type					base_type;

						const_OrderListIterator()					{}
						const_OrderListIterator(const typename OrderListParent::const_iterator& iter):
    						const_OrderListIterator::iterator_adaptor_(iter)
						{}
						const_OrderListIterator(const const_OrderListIterator<T>& other):
							const_OrderListIterator::iterator_adaptor_(other.base())
						{}
						const_OrderListIterator(const OrderListIterator<T>& other):
							const_OrderListIterator::iterator_adaptor_(other.base())
						{}

	private:
		friend class	boost::iterator_core_access;
		reference		dereference() const							{ return Parent::base_reference()->data; }
		const base_type&
						get_base()									{ return Parent::base_reference(); }

		friend class 	OrderList<T>;
};


#include "orderlist.hpp"

#endif // __ORDERLIST_H__
