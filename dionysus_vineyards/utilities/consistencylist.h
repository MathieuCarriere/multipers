/*
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2006
 */

#ifndef __CONSISTENCYLIST_H__
#define __CONSISTENCYLIST_H__

#include "orderlist.h"

#include <iterator>
#include <iostream>
#include <list>
#include "types.h"

#include <utilities/boost.h>
#include <boost/iterator/iterator_adaptor.hpp>
//#include <boost/type_traits/is_convertible.hpp>
//#include <boost/utility/enable_if.hpp>


template<class T> 		struct 	ConsistencyListNode;
template<class T>		class   ConsistencyListIterator;
template<class T>		class   const_ConsistencyListIterator;

/**
 * ConsistencyList adds consistency order to OrderList which remains unchanged
 * through the lifetime of the list.
 */
template<class T>
class ConsistencyList: public OrderList<ConsistencyListNode<T> >
{
	public:
		class 			OrderComparison;
		class 			LessThanComparison;
		class 			GreaterThanComparison;
		class 			ConsistencyComparison;

		/// \name Comparison types
		/// @{
		// Explicit typedefs for Doxygen
		//typedef			LessThanComparison							LessThanComparison;
		//typedef			GreaterThanComparison						GreaterThanComparison;
		//typedef			ConsistencyComparison						ConsistencyComparison;
		/// @}

		typedef			ConsistencyListNode<T>						NodeType;
		typedef			ConsistencyList<T>							Self;
		typedef			OrderList<NodeType >						Parent;

		typedef			T											value_type;
		typedef			T&											reference;
		typedef			const T&									const_reference;
		typedef			ConsistencyListIterator<T>					iterator;
		typedef			const_ConsistencyListIterator<T>			const_iterator;

						ConsistencyList(): sz(0)					{}
						~ConsistencyList() 							{ clear(); }

		/// \name Order operations
		void			swap(iterator i, iterator j);				///< Exchanges the order of simplices pointed to by i and j
		template<class BinaryPredicate>
		void			sort(BinaryPredicate cmp);					///< Sorts the elements in accordance with cmp

		/// \name Container operations
		/// @{
		// Explicit calls instead of using declarations for Doxygen
		iterator		push_back(const_reference x);
		iterator		insert(iterator predecessor, const_reference x);	///< Inserts x immediately after predecessor (has to be a valid iterator)
		void			erase(iterator x)							{ Parent::erase(x.get_base()); --sz; }

		void			clear()										{ return Parent::clear(); }
		bool			empty() const								{ return Parent::empty(); }
		SizeType		size()	const								{ return sz; }
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

	private:
		unsigned int	sz;
};

/// Basic comparison that LessThan and GreaterThan derive from
template<class T>
class ConsistencyList<T>::OrderComparison
{
	public:
		typedef			typename ConsistencyList<T>::const_iterator		ComparableType;

		int 			compare(ComparableType a, ComparableType b) const;				/// (-1,0,1) = a (precedes, ==, succeeds) b
};

/// Determines if the first element is less than the second one
template<class T>
class ConsistencyList<T>::LessThanComparison: public OrderComparison
{
	public:
		typedef			OrderComparison								Parent;
		typedef			typename Parent::ComparableType				ComparableType;

		int 			compare(ComparableType a, ComparableType b) const;
		bool 			operator()(ComparableType a, ComparableType b) const;
};

/// Determines if the first element is greater than the second one
template<class T>
class ConsistencyList<T>::GreaterThanComparison: public OrderComparison
{
	public:
		typedef			OrderComparison								Parent;
		typedef			typename Parent::ComparableType				ComparableType;

		int 			compare(ComparableType a, ComparableType b) const;
		bool 			operator()(ComparableType a, ComparableType b) const;
};

/// Determines the order of the two elements in the consistency order (that doesn't change during the execution)
template<class T>
class ConsistencyList<T>::ConsistencyComparison
{
	public:
		typedef			typename ConsistencyList<T>::const_iterator		ComparableType;

		int 			compare(ComparableType a, ComparableType b) const;				///< (-1,0,1) = a (precedes, ==, succeeds) b
		bool 			operator()(ComparableType a, ComparableType b) const;
};

/// Structure storing auxilliary information requred for each node of ConsistencyList
template<class T>
struct ConsistencyListNode
{
	ConsistencyListNode(const T& d, unsigned int c):
		data(d), consistency(c)
	{}

	T 				data;
	OrderType		consistency;

	std::ostream& 		operator<<(std::ostream& out) const					{ return out << consistency << ": " << data; }
};

template<class T>
std::ostream&			operator<<(std::ostream& out, const ConsistencyListNode<T>& n)	{ return n.operator<<(out); }

template<class T>
class ConsistencyListIterator: public boost::iterator_adaptor<ConsistencyListIterator<T>,
						 								typename ConsistencyList<T>::Parent::iterator,
						 								T>
{
	private:
		struct			enabler										{};

	public:
		typedef			typename ConsistencyList<T>::Parent			ConsistencyListParent;
		typedef 		boost::iterator_adaptor<ConsistencyListIterator<T>,
												typename ConsistencyListParent::iterator,
												T>					Parent;
		typedef			typename Parent::reference					reference;
		typedef			typename Parent::base_type					base_type;

						ConsistencyListIterator()					{}
						ConsistencyListIterator(const typename ConsistencyListParent::iterator& iter):
    						ConsistencyListIterator::iterator_adaptor_(iter)
						{}
						ConsistencyListIterator(const ConsistencyListIterator<T>& other):
							ConsistencyListIterator::iterator_adaptor_(other.base())
						{}

	private:
		friend class	boost::iterator_core_access;
		reference		dereference() const							{ return Parent::base_reference()->data; }
		base_type&		get_base()									{ return Parent::base_reference(); }

		friend class 	ConsistencyList<T>;
};

template<class T>
class const_ConsistencyListIterator: public boost::iterator_adaptor<const_ConsistencyListIterator<T>,
						 									  typename ConsistencyList<T>::Parent::const_iterator,
						 									  const T>
{
	private:
		struct			enabler										{};

	public:
		typedef			typename ConsistencyList<T>::Parent				ConsistencyListParent;
		typedef 		boost::iterator_adaptor<const_ConsistencyListIterator<T>,
												typename ConsistencyListParent::const_iterator,
												const T>			Parent;
		typedef			typename Parent::reference					reference;
		typedef			typename Parent::base_type					base_type;

						const_ConsistencyListIterator()					{}
						const_ConsistencyListIterator(const typename ConsistencyListParent::const_iterator& iter):
    						const_ConsistencyListIterator::iterator_adaptor_(iter)
						{}
						const_ConsistencyListIterator(const const_ConsistencyListIterator<T>& other):
							const_ConsistencyListIterator::iterator_adaptor_(other.base())
						{}
						const_ConsistencyListIterator(const ConsistencyListIterator<T>& other):
							const_ConsistencyListIterator::iterator_adaptor_(other.base())
						{}

	private:
		friend class	boost::iterator_core_access;
		reference		dereference() const							{ return Parent::base_reference()->data; }
		const base_type&
						get_base()									{ return Parent::base_reference(); }

		friend class 	ConsistencyList<T>::OrderComparison;
};


#include "consistencylist.hpp"

#endif // __CONSISTENCYLIST_H__
