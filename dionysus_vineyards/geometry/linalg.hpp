template<class ValueType_>
typename LinearAlgebra<ValueType_>::ValueType 					
LinearAlgebra<ValueType_>::
determinant(const MatrixType& a)
{
	using namespace boost::numeric::ublas;
	MatrixType m = a;
	lu_factorize(m);
	return determinant(triangular_adaptor<MatrixType, upper>(m));
}

template<class ValueType_>
void
LinearAlgebra<ValueType_>::
solve(const MatrixType& a, const VectorType& b, VectorType& x)
{
	using namespace boost::numeric::ublas;
	MatrixType m = a;
	x = b;
	lu_factorize(m);
	//lu_substitute(m, x);
	inplace_solve(m,x,unit_lower_tag());
	inplace_solve(m,x,upper_tag());
}


/* Private */
template<class ValueType_>
template<class TriangularType_>
typename LinearAlgebra<ValueType_>::ValueType 					
LinearAlgebra<ValueType_>::
determinant(const boost::numeric::ublas::triangular_adaptor<MatrixType, TriangularType_>& t)
{
	ValueType res = ValueType(1);
	for (typename MatrixType::size_type i = 0; i < t.size1(); ++i)
		res *= t(i,i);
	return res;
}
