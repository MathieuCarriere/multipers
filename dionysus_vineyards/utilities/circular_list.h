/**
 * Singly-linked circular list implementation by Donovan Rebbechi <elflord@pegasus.rutgers.edu>.
 * Tiny modifications by DM to make ISO C++ compliant (i.e., so that the compiler stops complaining), 
 * added size(), and boost serialization
 */

#ifndef __CIRCULAR_LIST_H__
#define __CIRCULAR_LIST_H__

#include <iterator>

#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/binary_object.hpp>

#include "types.h"


template <class T>
class List
{
	struct Node
	{
		Node(const T& x,Node* y = 0):m_data(x),m_next(y){}
		T m_data;
		Node* m_next;
	
		private:
			// Serialization
			friend class boost::serialization::access;

			template<class Archive>
			void serialize(Archive& ar, version_type )
			{
				ar & make_nvp("data", m_data);
				ar & make_nvp("next", m_next);
			}
	};

	Node* m_node;
	size_t sz;
	
	private:
		// Serialization
		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive& ar, version_type )
		{ ar & make_nvp("root", m_node); }

public:
	class const_iterator;
	
	class iterator
		: public std::iterator<std::forward_iterator_tag, T>
	{
		Node* m_rep;
	public:
		friend class List;
		friend class const_iterator;
		typedef T value_type;
		typedef T& reference;
		typedef T* pointer;
		typedef int difference_type;
		typedef std::forward_iterator_tag iterator_category;

		inline iterator(Node* x=0):m_rep(x){}
		inline iterator(const iterator& x):m_rep(x.m_rep) {}
		inline iterator(const const_iterator& x):m_rep(const_cast<Node*>(x.m_rep)) {}			// not very safe
		inline iterator& operator=(const iterator& x)
		{ 
			m_rep=x.m_rep; return *this; 
		}
		inline iterator& operator++()
		{ 
			m_rep = m_rep->m_next; return *this; 
		}
		inline iterator operator++(int)
		{ 
			iterator tmp(*this); m_rep = m_rep->m_next; return tmp; 
		}
		inline reference operator*() const { return m_rep->m_next->m_data; }
		inline pointer operator->() const { return m_rep->m_next; }
		inline bool operator==(const iterator& x) const
		{
			return m_rep == x.m_rep; 
		}	
		inline bool operator!=(const iterator& x) const
		{
			return m_rep != x.m_rep; 
		}	

	};

	class const_iterator
		: public std::iterator<std::forward_iterator_tag, const T>
	{
		const Node* m_rep;
	public:
		friend class List;
		friend class iterator;
		typedef T value_type;
		typedef T& reference;
		typedef T* pointer;
		typedef int difference_type;
		typedef std::forward_iterator_tag iterator_category;
		inline const_iterator(const Node* x=0):m_rep(x){}
		inline const_iterator(const const_iterator& x):m_rep(x.m_rep) {}
		inline const_iterator(const iterator& x):m_rep(x.m_rep){}
		inline const_iterator& operator=(const const_iterator& x)
		{ 
			m_rep=x.m_rep; return *this; 
		}
		inline const_iterator& operator=(const iterator& x)
		{ 
			m_rep=x.m_rep; return *this; 
		}		
		inline const_iterator& operator++()
		{ 
			m_rep = m_rep->m_next; return *this; 
		}
		inline const_iterator operator++(int)
		{ 
			const_iterator tmp(*this); m_rep = m_rep->m_next; return tmp; 
		}
		inline reference operator*() const { return m_rep->m_next->m_data; }
		inline pointer operator->() const { return m_rep->m_next; }
		inline bool operator==(const const_iterator& x) const
		{
			return m_rep == x.m_rep; 
		}
		inline bool operator!=(const const_iterator& x) const
		{
			return m_rep != x.m_rep; 
		}
	};

	typedef			T&				reference;
	typedef			const T&		const_reference;
	typedef			T*				pointer;
	typedef			const T*		const_pointer;

	List() : m_node(new Node(T())), sz(0) { m_node->m_next = m_node; }

	List(const List& L) : m_node(new Node(T())), sz(L.sz)
	{
		m_node->m_next=m_node;
		for ( const_iterator i = L.begin(); i!= L.end(); ++i )
			push_front(*i);
		reverse();
	}

	void reverse()
	{
		if (empty())
			return;
		Node* new_m_node = m_node->m_next->m_next;
		Node* p = m_node; Node* i = m_node->m_next; Node* n;
		do  
		{
			n = i->m_next;
			i->m_next = p;
			p = i; i = n;
		}
		while (p!=m_node);
		m_node = new_m_node;
	}

	void swap(List& x)
	{
		std::swap(m_node, x.m_node);
		std::swap(sz, x.sz);
	}

	List& operator=(const List& x)
	{
		List tmp(x);
		swap(tmp);
		return *this;
	}

	~List() { clear(); delete m_node; }
	void clear() { while (!empty()) pop_front(); sz = 0; }


	inline size_t size() const					{ return sz; }

	inline void push_front(const T&x)
	{
		insert (begin(),x);
	}
	inline void push_back(const T& x)
	{
		insert (end(),x);
	}
	inline void pop_front()
	{
		erase(begin());
	}
	inline bool empty() const { return m_node == m_node->m_next; }

	inline T& front() { return *begin(); }
	inline const T& front() const { return *begin(); }
	inline T& back() { return m_node->m_data; }
	inline const T& back() const { return m_node->m_data; }

	inline iterator begin() { return iterator(m_node->m_next); }
	inline iterator end() { return iterator(m_node); }
	inline const_iterator begin() const 
	{ 
		return const_iterator(m_node->m_next); 
	}
	inline const_iterator end() const 
	{ 
		return const_iterator(m_node); 
	}

	iterator erase (iterator x)
	{
		if (x == end())
			return x;
		
		if (x.m_rep->m_next == m_node)
			m_node = x.m_rep;

		Node* tmp = x.m_rep->m_next;
		x.m_rep->m_next = x.m_rep->m_next->m_next;
		delete tmp;
		--sz;

		return x;
	}

	iterator insert (iterator x, const T& y)
	{
		Node* tmp = new Node(y,x.m_rep->m_next);
		if (x.m_rep == m_node)		// if inserting before end
			m_node = tmp;
		x.m_rep->m_next = tmp;
		++sz;
		
		return x++;					// x points to the new data, return it, and slide it over
	}

	// rotate x to beginning
	void rotate (iterator x)
	{
		if (x == end())
			return;
		Node* sentinel = m_node->m_next;
		m_node->m_next = m_node->m_next->m_next;
		sentinel->m_next = x.m_rep->m_next;
		x.m_rep->m_next = sentinel;
		m_node = x.m_rep; 
	}
};

#endif // __CIRCULAR_LIST_H__
