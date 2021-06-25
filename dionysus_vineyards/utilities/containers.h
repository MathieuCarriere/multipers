#ifndef __CONTAINERS_H__
#define __CONTAINERS_H__

#include "circular_list.h"

#if DEBUG_CONTAINERS
    #include <debug/vector>
    #include <debug/deque>
    using std::__debug::vector;
    using std::__debug::deque;
#else
    #include <vector>
    #include <deque>
    using std::vector;
    using std::deque;
#endif

// TODO: write documentation

template<class Container_, class Comparison_ = std::less<typename Container_::value_type> >
struct ContainerTraits
{
    typedef     Container_                                                                  Container;
    typedef     typename Container::value_type                                              Item;
    typedef     typename Container::const_reference                                         const_reference;
    typedef     Comparison_                                                                 Comparison;

    static void reserve(Container& c, size_t sz)                                            {}
    static void sort(Container& c, const Comparison& cmp = Comparison())                    { c.sort(cmp); }
    static void push_front(Container& c, const_reference x)                                 { c.push_front(x); }
};

/**
 * Class: SizeStorage
 *
 * This class expresses how size should be stored for various containers to be able 
 * to deduce it in constant time. By default it is stored explicitly, so that if the 
 * Container is std::list everything works. However, specialization is available for 
 * std::vector which uses its builtin size() function.
 */
template<class Container_>
class SizeStorage
{
    public:
        typedef         Container_                                                          Container;
        typedef         SizeStorage<Container>                                              Self;

                        SizeStorage(size_t size = 0): size_(size)                           {}

        Self&           operator+=(size_t inc)                                              { size_ += inc; return *this; }
        Self&           operator-=(size_t dec)                                              { size_ -= dec; return *this; }
        Self&           operator++()                                                        { ++size_; return *this; }
        Self            operator++(int)                                                     { Self tmp = *this; size_++; return tmp; }
        Self&           operator--()                                                        { --size_; return *this; }
        Self            operator--(int)                                                     { Self tmp = *this; size_--; return tmp; }
        size_t          size(const Container& c) const                                      { return size_; }
        void            swap(SizeStorage& other)                                            { std::swap(size_, other.size_); }

        void            clear()                                                             { size_ = 0; }

    private:
        size_t          size_;
};

/**
 * Class: CountingBackInserter<C>
 *
 * Derives from std::back_insert_iterator<C> and SizeStorage<C>, 
 * and keeps track of the number of inserted elements.
 */
template<class C>
struct CountingBackInserter: public std::back_insert_iterator<C>, 
                             public SizeStorage<C>
{
    typedef                     CountingBackInserter                            Self;
    typedef                     std::back_insert_iterator<C>                    ParentIterator;
    typedef                     SizeStorage<C>                                  ParentSize;

                                CountingBackInserter(C& c):
                                    ParentIterator(c)                           {}

    Self&                       operator++()                                    { ParentSize::operator++(); ParentIterator::operator++(); return *this; }
    Self                        operator++(int i)                               { Self tmp = *this; ParentSize::operator++(i); ParentIterator::operator++(i); return tmp; }
};

/**
 * Class: PushBackFunctor<Container>
 *
 * Performs the same task as std::back_insert_iterator<Container>, but as a functor.
 */
template<class Container_>
class PushBackFunctor
{
    public:
        typedef                 Container_                                      Container;
        typedef                 typename Container::value_type                  value_type;

                                PushBackFunctor(Container& container):
                                    container_(container)                       {}

        void                    operator()(const value_type& v) const           { container_.push_back(v); }

    private:
        Container&              container_;
};

template<class Container>
PushBackFunctor<Container>
make_push_back_functor(Container& container)                                    { return PushBackFunctor<Container>(container); }

/**
 * Class: InsertFunctor<Container>
 *
 * Performs insertions of its arguments into the given container.
 */
template<class Container_>
class InsertFunctor
{
    public:
        typedef                 Container_                                      Container;
        typedef                 typename Container::value_type                  value_type;

                                InsertFunctor(Container& container):
                                    container_(container)                       {}

        void                    operator()(const value_type& v) const           { container_.insert(v); }

    private:
        Container&              container_;
};

template<class Container>
InsertFunctor<Container>
make_insert_functor(Container& container)                                       { return InsertFunctor<Container>(container); }

/* Specializations */

template<class T, class Comparison_>
struct ContainerTraits<vector<T>, Comparison_>
{
    typedef     T                                                                           Item;
    typedef     vector<T>                                                                   Container;
    typedef     typename Container::const_reference                                         const_reference;
    typedef     Comparison_                                                                 Comparison;

    static void reserve(Container& c, size_t sz)                                            { c.reserve(sz); }
    static void sort(Container& c, const Comparison& cmp = Comparison())                    { std::sort(c.begin(), c.end(), cmp); }
    static void push_front(Container& c, const_reference x)                                 { c.insert(c.begin(), x); }
};

template<class T, class Comparison_>
struct ContainerTraits<deque<T>, Comparison_>
{
    typedef     T                                                                           Item;
    typedef     deque<T>                                                                    Container;
    typedef     typename Container::const_reference                                         const_reference;
    typedef     Comparison_                                                                 Comparison;

    static void reserve(Container& c, size_t sz)                                            { }
    static void sort(Container& c, const Comparison& cmp = Comparison())                    { std::sort(c.begin(), c.end(), cmp); }
    static void push_front(Container& c, const_reference x)                                 { c.push_front(x); }
};

template<class T, class Comparison_>
struct ContainerTraits<List<T>, Comparison_>
{
    typedef     T                                                                           Item;
    typedef     List<T>                                                                     Container;
    typedef     typename Container::const_reference                                         const_reference;
    typedef     Comparison_                                                                 Comparison;

    static void reserve(Container& c, size_t sz)                                            { }
    static void sort(Container& c, const Comparison& cmp = Comparison())                    
    { 
        vector<Item> tmp(c.begin(), c.end());
        std::sort(tmp.begin(), tmp.end(), cmp);
        std::copy(tmp.begin(), tmp.end(), c.begin());
    }
    static void push_front(Container& c, const_reference x)                                 { c.push_front(x); }
};

// TODO: specialize for List (singly-linked list)

template<class T>
class SizeStorage<vector<T> >
{
    public:
        typedef         vector<T>                                                           Container;
        typedef         SizeStorage<Container>                                              Self;

                        SizeStorage()                                                       {}
                        SizeStorage(size_t)                                                 {}

        Self&           operator+=(size_t)                                                  { return *this; }
        Self&           operator-=(size_t)                                                  { return *this; }
        Self&           operator++()                                                        { return *this; }
        Self            operator++(int)                                                     { return *this; }
        Self&           operator--()                                                        { return *this; }
        Self            operator--(int)                                                     { return *this; }
        size_t          size(const Container& c) const                                      { return c.size(); }
        void            swap(SizeStorage&)                                                  {}

        void            clear()                                                             {}
};

template<class T>
class SizeStorage<deque<T> >
{
    public:
        typedef         deque<T>                                                            Container;
        typedef         SizeStorage<Container>                                              Self;

                        SizeStorage()                                                       {}
                        SizeStorage(size_t)                                                 {}

        Self&           operator+=(size_t)                                                  { return *this; }
        Self&           operator-=(size_t)                                                  { return *this; }
        Self&           operator++()                                                        { return *this; }
        Self            operator++(int)                                                     { return *this; }
        Self&           operator--()                                                        { return *this; }
        Self            operator--(int)                                                     { return *this; }
        size_t          size(const Container& c) const                                      { return c.size(); }
        void            swap(SizeStorage&)                                                  {}

        void            clear()                                                             {}
};

#endif // __CONTAINERS_H__
