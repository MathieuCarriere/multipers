#include <utilities/log.h>

template<class T>
void
UPolynomial<T>::
solve(const RationalFunctionP& rf, RootStack& stack)
{
	typedef SYNAPS::Seq<RootType>		RootSeq;

    AssertMsg(stack.empty(), "Stack must be empty before solve");

	RootSeq seq_num = SYNAPS::solve(SynapsTraits<T>::convert(rf.numerator()), Solver());
	RootSeq seq_den = SYNAPS::solve(SynapsTraits<T>::convert(rf.denominator()), Solver());

	// TODO: assert that all roots in seq_den have positive multiplicity
	// TODO: deal with multiplicities for the numerator
	for (typename RootSeq::const_reverse_iterator cur = seq_num.rbegin(); 
												  cur != seq_num.rend(); 
												  ++cur)
	{
		if (SynapsTraits<T>::multiplicity(*cur) % 2 != 0)
		{		
			if (!stack.empty() && stack.top() == *cur)			// get rid of even multiplicities
				// TODO: add logging information for this
				stack.pop();
			else
				stack.push(*cur);
		}
	}
}

template<class T>
int
UPolynomial<T>::
sign_at(const RationalFunctionP& rf, const RootType& r)
{
	return SynapsTraits<T>::sign_at(rf.numerator(), r) * SynapsTraits<T>::sign_at(rf.denominator(), r);
}

template<class T>
int
UPolynomial<T>::
sign_at_negative_infinity(const RationalFunctionP& rf)
{
	const Polynomial& num = rf.numerator();
	const Polynomial& den = rf.denominator();
	int ndegree = num.get_degree();
	int ddegree = den.get_degree();
    if (ndegree == -1) return 0;          // ndegree == -1 => num == 0, and 0 is 0 at -infinity
	return !((((ndegree + 1) % 2 == 0) ^ (num[ndegree] > 0)) ^
		     (((ddegree + 1) % 2 == 0) ^ (den[ddegree] > 0))) ? 1:-1;
}

SynapsTraits<QQ>::RootType
SynapsTraits<QQ>::
root(const CoefficientType& r)
{
	ZZ p[2] = { -SYNAPS::numerator(r), SYNAPS::denominator(r) };
	return SYNAPS::solve(SolverPolynomial(2, p), Solver(), 0); 
}

SynapsTraits<QQ>::SolverPolynomial
SynapsTraits<QQ>::
convert(const Polynomial& p)
{
	SolverPolynomial result(1, p.size() - 1);
	assert(result.size() == p.size());
	ZZ denominator_product = 1;
	for (unsigned int i = 0; i < p.size(); ++i)
		denominator_product *= SYNAPS::denominator(p[i]);
	for (unsigned int i = 0; i < p.size(); ++i)
		result[i] = SYNAPS::numerator(p[i]) * denominator_product / 
					SYNAPS::denominator(p[i]);
	return result;
}

