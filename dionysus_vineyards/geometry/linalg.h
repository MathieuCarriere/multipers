#ifndef __LINALG_H__
#define __LINALG_H__

#include <vector>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

template<class ValueType_>
class LinearAlgebra
{
	public:
		typedef 					ValueType_												ValueType;
		typedef 					boost::numeric::ublas::matrix<ValueType>				MatrixType;
		typedef 					boost::numeric::ublas::vector<ValueType>				VectorType;
	

		/* Currently don't need any of this */
		static ValueType 			determinant(const MatrixType& a);
		static void					solve(const MatrixType& a, const VectorType& b, VectorType& x);

	private:
		template<class TriangularType_>
		static ValueType			determinant(const boost::numeric::ublas::triangular_adaptor<MatrixType, TriangularType_>& t);
};

#include "linalg.hpp"

#endif
