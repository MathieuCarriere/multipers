#include <algorithm>
#include <vector>
#include "utilities/containers.h"

#include <boost/serialization/split_member.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/binary_object.hpp>
#include <utilities/boost.h>

#include "utilities/log.h"
#include "utilities/counter.h"

using boost::serialization::make_nvp;
using boost::serialization::make_binary_object;

#ifdef LOGGING
static rlog::RLogChannel* rlChain =                 DEF_CHANNEL( "topology/chain", rlog::Log_Debug);
#endif // LOGGING

#ifdef COUNTERS
static Counter*  cChainAddBasic =                   GetCounter("chain/add/basic");
static Counter*  cChainAddComparison =              GetCounter("chain/add/comparison");
#endif // COUNTERS

template<class C>
ChainWrapper<C>::
ChainWrapper()
{}

template<class C>
ChainWrapper<C>::
ChainWrapper(const ChainWrapper& c): 
    ChainRepresentation(c), Size(c)
{}

template<class C>
template<class Iterator>
ChainWrapper<C>::
ChainWrapper(Iterator bg, Iterator end): 
    ChainRepresentation(bg, end), Size(ChainRepresentation::size())
{}

template<class C>
template<class ConsistencyCmp>
void
ChainWrapper<C>::
append(const_reference x, const ConsistencyCmp& cmp)                        
{ 
    Size::operator++();

    // First try the special cases that x goes at the end
    if (empty() || cmp(back(), x))
        push_back(x); 
    // Then try the special case that x goes at the front
    else if (cmp(x, front()))
        ContainerTraits<C,ConsistencyCmp>::push_front(*this, x);
    else
        insert(std::upper_bound(begin(), end(), x, cmp), x);
}
        
template<class C>
template<class OrderComparison>
typename ChainWrapper<C>::const_reference               
ChainWrapper<C>::
top(const OrderComparison& cmp) const
{ 
    AssertMsg(!empty(), "Chain must not be empty for low()");
    return *std::min_element(begin(), end(), cmp);
}

template<class C>
void 
ChainWrapper<C>::
swap(ChainWrapper& c)
{
    ChainRepresentation::swap(c);
    Size::swap(c);
}

template<class C>
void 
ChainWrapper<C>::
clear()
{
    ChainRepresentation::clear();
    Size::clear();
}

template<class C>
template<class ConsistencyComparison>
void 
ChainWrapper<C>::
sort(const ConsistencyComparison& cmp)
{ 
    ContainerTraits<C,ConsistencyComparison>::sort(*this, cmp);
}

template<class C>
boost::optional<typename ChainWrapper<C>::const_iterator>
ChainWrapper<C>::
contains(const_reference x) const
{
    const_iterator res = std::find(begin(), end(), x);
    return boost::make_optional(res != end(), res);
}

template<class C>
boost::optional<typename ChainWrapper<C>::iterator>
ChainWrapper<C>::
contains(const_reference x)
{
    iterator res = std::find(begin(), end(), x);
    return boost::make_optional(res != end(), res);
}

template<class C>
bool
ChainWrapper<C>::
remove_if_contains(const_reference x)
{
    boost::optional<iterator> i = contains(x);
    if (i)
    {
        remove(*i);
        return true;
    } else
        return false;
}

template<class C>
template<class OutputMap>
std::string
ChainWrapper<C>::
tostring(const OutputMap& outmap) const
{
    std::string str;
    for (const_iterator cur = begin(); cur != end(); ++cur)
    {
        if (cur != begin()) str += ", ";
        str += outmap(*cur);
    }
    return str;
}

template<class C>
template<class ConsistencyCmp>
typename ChainWrapper<C>::Self& 
ChainWrapper<C>::
add(const Self& c, const ConsistencyCmp& cmp)
{
    // TODO: tmp-based addition is necessary and useful for Containers that are vectors, 
    //       however, I believe it creates costly overhead for Containers that are lists.
    //       Need to put some thought into this.
    ChainRepresentation     tmp;

    CountingBackInserter<ChainRepresentation> bi(tmp);
    std::set_symmetric_difference(begin(), end(), c.begin(), c.end(), bi, cmp);

    CountBy(cChainAddBasic, size() + c.size() - (size() + c.size() - tmp.size())/2);
    
    static_cast<ChainRepresentation*>(this)->swap(tmp);
    static_cast<Size*>(this)->swap(bi);   

    return *this;
}

template<class C>
template<class OrderComparison>
typename ChainWrapper<C>::const_reference 
ChainWrapper<C>::
get_first(const OrderComparison& cmp) const
{ return top(cmp); }

        
template<class C>
template<class Archive> 
void                        
ChainWrapper<C>::
serialize(Archive& ar, boost::serialization::version_type )
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Parent);
}
