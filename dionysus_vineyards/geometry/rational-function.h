#ifndef __RATIONAL_FUNCTION_H__
#define __RATIONAL_FUNCTION_H__

#include <iostream>
#include "number-traits.h"

template <class Polynomial_>
class RationalFunction
{
	public:
		typedef 				Polynomial_													Polynomial;
		typedef					typename Polynomial::coeff_t								CoefficientType;
		typedef					typename Polynomial::value_type								ValueType;

		/// \name Constructors
		/// @{
								RationalFunction():		
									numerator_(CoefficientType(0)), 
									denominator_(CoefficientType(1)) 						{}
								RationalFunction(const Polynomial& p):
									numerator_(p), 
									denominator_(CoefficientType(1))						{}
								RationalFunction(const Polynomial& num, const Polynomial& denom):
									numerator_(num), denominator_(denom)					{ normalize(); }
								RationalFunction(const RationalFunction& other):
									numerator_(other.numerator_), 
									denominator_(other.denominator_)						{ normalize(); }
		/// @}

		/// \name Operators
		/// @{
		RationalFunction		operator-()	const;
		RationalFunction		operator+(const RationalFunction& o) const;
		RationalFunction		operator-(const RationalFunction& o) const;
		RationalFunction		operator*(const RationalFunction& o) const;
		RationalFunction		operator/(const RationalFunction& o) const;
		RationalFunction		operator+(const CoefficientType& a) const;
		RationalFunction		operator-(const CoefficientType& a) const;
		RationalFunction		operator*(const CoefficientType& a) const;
		RationalFunction		operator/(const CoefficientType& a) const;
		/// @}
		
		/// \name Modifiers
		/// @{
		RationalFunction& 		operator+=(const RationalFunction& o);
		RationalFunction&		operator-=(const RationalFunction& o);
		RationalFunction&		operator*=(const RationalFunction& o);
		RationalFunction&		operator/=(const RationalFunction& o);
		/// @}
		
		/// \name Assignment
		/// @{
		RationalFunction&		operator=(const RationalFunction& o);
		//RationalFunction&		operator=(const Polynomial& o);
		/// @}
		
		/// \name Evaluation
		/// @{
		ValueType				operator()(const ValueType& t) const;
		bool 					operator==(const RationalFunction& o) const;
		bool 					operator!=(const RationalFunction& o) const					{ return !operator==(o); }
		/// @}
								
		/// \name Accessors
		/// @{
		const Polynomial&		numerator() const											{ return numerator_; }
		const Polynomial&		denominator() const											{ return denominator_; }
		/// @}
		
		RationalFunction&		normalize();
						
	private:
		Polynomial				numerator_, denominator_;
};

template<class P> 
struct number_traits<RationalFunction<P> > 
{
	typedef 			RationalFunction<P>				NumberType;
	static NumberType&	normalize(NumberType& n)		{ return n.normalize(); }
};


template<class Polynomial_>
std::ostream& 
operator<<(std::ostream& out, const RationalFunction<Polynomial_>& r)
{ return out << r.numerator() << " / " << r.denominator(); }//   << ", gcd: " << gcd(r.numerator(), r.denominator()); }

template<class Polynomial_>
inline RationalFunction<Polynomial_>
operator*(const typename RationalFunction<Polynomial_>::CoefficientType& a, const RationalFunction<Polynomial_>& r)
{ return (r * a); }

template<class Polynomial_>
inline RationalFunction<Polynomial_>
operator+(const typename RationalFunction<Polynomial_>::CoefficientType& a, const RationalFunction<Polynomial_>& r)
{ return (r + a); }

template<class Polynomial_>
inline RationalFunction<Polynomial_> 
operator-(const typename RationalFunction<Polynomial_>::NT& a, const RationalFunction<Polynomial_>& r)
{ return -(r - a); }


#include "rational-function.hpp"

#endif // __RATIONAL_FUNCTION_H__
