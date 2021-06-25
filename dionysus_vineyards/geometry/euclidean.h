/*
 * Author: Dmitriy Morozov
 * Department of Computer Science, Duke University, 2005 -- 2007
 */

#ifndef __EUCLEDIAN_H__
#define __EUCLEDIAN_H__

#include <utility>
#include <vector>
#include <algorithm>

#include "linalg.h"
#include "number-traits.h"


/**
 * Geometric Kernel. Defines operations on geometric primitives.
 * \ingroup geometry
 */
template<class NumberType_ = double>
class Kernel
{
	public:
		typedef 					unsigned int								DimensionType;
		typedef						NumberType_									NumberType;
		typedef						LinearAlgebra<NumberType>					LinearAlgebraK;
		typedef						typename LinearAlgebraK::MatrixType			MatrixType;
		typedef						typename LinearAlgebraK::VectorType			VectorType;

		class						Point;
		class						Sphere;
		typedef						std::vector<const Point*>					PointContainer;
		

									Kernel(DimensionType dimension);

		DimensionType				dimension() const							{ return dimension_; }
		Point						origin() const								{ return origin_; }

		/** Returns matrix describing the equation of a circumsphere of points */
		Sphere						circumsphere(const PointContainer& points) const;
		
		/** Returns squared radius of the circumsphere */
		NumberType					circumradius(const PointContainer& points) const;
		
		/** Returns center of the circumsphere */
		Point						circumcenter(const PointContainer& points) const;

		/** The result is positive if points[0] lies outside the circumsphere of points,
			0 if points[0] is on the circumsphere, and negative if it's inside */
		NumberType					side_of_circumsphere(const PointContainer& points, const Point& p) const;	

	private:
		NumberType&					normalize(NumberType& n) const;
		Point&						normalize(Point& p) const;

		DimensionType				dimension_;
		Point						origin_;
};


/** 
 * Point class.
 * \ingroup geometry 
 */
template<class NumberType_>
class Kernel<NumberType_>::Point: public VectorType
{
	public:
		typedef						VectorType									Parent;
		typedef						NumberType_									NumberType;

									Point(DimensionType d):	Parent(d)			{}
		template<class Vec>			Point(const Vec& v): Parent(v)				{}
									Point(const Point& p, const NumberType& pp);


									//operator VectorType() const					{ return *this; }

		NumberType					squared_distance(const Point& p) const;

		using						Parent::size;
};


/** 
 * Sphere class.
 * \ingroup geometry
 */
template<class NumberType_>
class Kernel<NumberType_>::Sphere
{
	public:
										Sphere(const Point& center, 
											   const NumberType& squared_radius):
											center_(center), squared_radius_(squared_radius)
										{}
										
		/** The result is positive if p lies outside the sphere, 
			0 if p is on the sphere, and negative if it's inside */
		NumberType						side_of(const Point& p) const;

		const Point&					center() const								{ return center_; }
		const NumberType&				squared_radius() const						{ return squared_radius_; }

	private:
		Point							center_;
		NumberType						squared_radius_;
};


#include "euclidean.hpp"

#endif // __EUCLEDIAN_H__
