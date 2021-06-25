#ifndef __TYPES_H__
#define __TYPES_H__

#include <limits>
#include <iostream>

/* Types */
typedef 	bool					Sign;
typedef		short int		        Dimension;
const 		Sign	 				POS = true;
const 		Sign					NEG = false;
typedef		double					RealType;
typedef		unsigned int			SizeType;

static RealType Infinity = std::numeric_limits<RealType>::infinity();

typedef 	const unsigned int&		version_type;

// Empty is made a template so that we don't have to compile and deal with a library
// solely for its operator<<(out, e) function
template<typename T = void>
struct      Empty                   {};


struct      use_default             {};

template<class T, class Default>
struct      if_default
{ typedef   T           type; };

template<class Default>
struct      if_default<use_default, Default>
{ typedef   Default     type; };


template<typename T>
std::ostream& operator<<(std::ostream& out, Empty<T> e) { return out; }

enum        SwitchType
{
            DiffDim     = 0,
            Case1       = 0x4,
            Case12      = 0x5,
            Case112     = 0x6,
            Case2       = 0x8,
            Case212     = 0x9,
            Case3       = 0x10,
            Case31      = 0x11,
            Case4       = 0x20,
};

// Nothing to do for serializing Empty, but still need to provide this function
namespace boost {
namespace serialization {

template<class Archive, class T>
void serialize(Archive & ar, Empty<T>&, const unsigned int )
{}

} // namespace serialization
} // namespace boost

#endif // __TYPES_H__
