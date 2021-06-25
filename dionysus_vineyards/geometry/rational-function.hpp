/* Operators */
template<class P>
RationalFunction<P>
RationalFunction<P>::
operator-()	const											
{ return RationalFunction(-numerator_,denominator_); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator+(const RationalFunction& o) const	
{ return RationalFunction(numerator_*o.denominator_ + o.numerator_*denominator_, denominator_*o.denominator_); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator-(const RationalFunction& o) const
{ RationalFunction tmp(*this); tmp -= o; return tmp; }
//{ return RationalFunction(numerator_*o.denominator_ - o.numerator_*denominator_, denominator_*o.denominator_); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator*(const RationalFunction& o) const
{ return RationalFunction(numerator_*o.numerator_, denominator_*o.denominator_); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator/(const RationalFunction& o) const
{ return RationalFunction(numerator_*o.denominator_, denominator_*o.numerator_); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator+(const typename RationalFunction<P>::CoefficientType& a) const
{ return RationalFunction(numerator_ + a*denominator_, denominator_); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator-(const typename RationalFunction<P>::CoefficientType& a) const
{ return operator+(-a); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator*(const typename RationalFunction<P>::CoefficientType& a) const
{ return RationalFunction(a*numerator_, denominator_); }

template<class P>
RationalFunction<P>
RationalFunction<P>::
operator/(const typename RationalFunction<P>::CoefficientType& a) const
{ return RationalFunction(numerator_, a*denominator_); }

template<class P>
RationalFunction<P>&
RationalFunction<P>::
operator+=(const RationalFunction& o)
{
 	numerator_ *= o.denominator_;
	numerator_ += o.numerator_*denominator_;
	denominator_ *= o.denominator_;
	return *this;
}

template<class P>
RationalFunction<P>&
RationalFunction<P>::
operator-=(const RationalFunction& o)
{
 	numerator_ *= o.denominator_;
	numerator_ -= o.numerator_*denominator_;
	denominator_ *= o.denominator_;
	return *this;
}

template<class P>
RationalFunction<P>&
RationalFunction<P>::
operator*=(const RationalFunction& o)
{
 	numerator_ *= o.numerator_;
	denominator_ *= o.denominator_;
	return *this;
}

template<class P>
RationalFunction<P>&
RationalFunction<P>::
operator/=(const RationalFunction& o)
{
 	numerator_ *= o.denominator_;
	denominator_ *= o.numerator_;
	return *this;
}

template<class P>
RationalFunction<P>&
RationalFunction<P>::
operator=(const RationalFunction& o)
{
	numerator_ = o.numerator_;
	denominator_ = o.denominator_;
	return *this;
}

#if 0
template<class P>
RationalFunction<P>&
RationalFunction<P>::
operator=(const Polynomial& o)
{
	numerator_ = o;
	denominator_ = 1;
	return *this;
}
#endif

/* Evaluation */
template<class P>
typename RationalFunction<P>::ValueType
RationalFunction<P>::
operator()(const typename RationalFunction<P>::ValueType& t) const
{ return numerator_(t)/denominator_(t); }

template<class P>
bool
RationalFunction<P>::
operator==(const RationalFunction& o) const
{ return (numerator_ == o.numerator_) && (denominator_ == o.denominator_); }

template<class P>
RationalFunction<P>&
RationalFunction<P>::
normalize()
{
#if 0
	std::cout << "This:                      " << std::flush << this << std::endl;
	std::cout << "Numerator address:         " << std::flush << &numerator_ << std::endl;
	std::cout << "Denominator address:       " << std::flush << &denominator_ << std::endl;
	std::cout << "Normalizing (numerator):   " << std::flush << numerator_ << std::endl;
	std::cout << "Normalizing (denominator): " << std::flush << denominator_ << std::endl;
#endif
	Polynomial divisor = gcd(numerator_, denominator_);
	numerator_ /= divisor;
	denominator_ /= divisor;
	return *this;
}
