#ifndef __NUMBER_TRAITS_H__
#define __NUMBER_TRAITS_H__

template<class NumberType_>
class number_traits
{
	public:
		typedef							NumberType_									NumberType;

		static NumberType&				normalize(NumberType& n)					{ return n; }
};

#endif
