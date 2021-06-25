/* --- Point --- */
template<class NumberType_>
Kernel<NumberType_>::
Kernel(DimensionType dimension): dimension_(dimension), origin_(dimension)
{
	for (unsigned int i = 0; i < dimension; ++i)
		origin_(i) = NumberType(0);
}


template<class NumberType_>
Kernel<NumberType_>::Point::
Point(const Point& p, const NumberType& pp): Parent(p.size() + 1)
{
	using boost::numeric::ublas::vector_range;
	using boost::numeric::ublas::range;
	
	vector_range<VectorType> vr(*this, range(0, size() - 1));
	vr = p;
	(*this)(size() - 1) = pp;
}

template<class NumberType_>
typename Kernel<NumberType_>::NumberType
Kernel<NumberType_>::Point::
squared_distance(const Point& p) const
{
	return boost::numeric::ublas::inner_prod(*this - p, *this - p);
}

template<class NumberType_>
std::istream&
operator>>(std::istream& in, typename Kernel<NumberType_>::Point& p)
{
	for (unsigned int i = 0; i < p.size(); ++i)
		in >> p[i];
	return in;
}

template<class NumberType_>
std::ostream&
operator<<(std::ostream& out, typename Kernel<NumberType_>::Point& p)
{
	for (unsigned int i = 0; p.size(); ++i)
		out << p[i];
	return out;
}

/* --- Sphere --- */
template<class NumberType_>
typename Kernel<NumberType_>::NumberType
Kernel<NumberType_>::Sphere::
side_of(const Point& p) const				
{ return p.squared_distance(center_) - squared_radius(); }


/* --- Kernel --- */
template<class NumberType_>
typename Kernel<NumberType_>::Sphere
Kernel<NumberType_>::
circumsphere(const PointContainer& points) const
{
	if (points.size() == 0) return Sphere(origin(), NumberType(0));		// FIXME: should this be an assertion instead?

	using boost::numeric::ublas::inner_prod;
	using boost::numeric::ublas::vector_range;
	using boost::numeric::ublas::range;

	std::vector<Point> basis;
	basis.reserve(points.size() - 1);
	for (unsigned int i = 1; i < points.size(); ++i)
	{
		Point pt = *points[i] - *points[0];
		basis.push_back(Point(pt, inner_prod(pt, pt)));
		for (unsigned int j = 0; j < i - 1; ++j)
		{
			basis[i-1] -= basis[j] * inner_prod(basis[i-1], basis[j]) / inner_prod(basis[j], basis[j]);
			normalize(basis[i-1]);
		}
	}
	Point clifted(origin(), NumberType(1));
	for (unsigned int j = 0; j < basis.size(); ++j)
	{
		clifted -= basis[j] * inner_prod(clifted, basis[j]) / inner_prod(basis[j], basis[j]);
		normalize(clifted);
	}
	
	Point center(vector_range<VectorType>(clifted, range(0, dimension())));
	center /= NumberType_(-2) * clifted(dimension());
	NumberType squared_radius = inner_prod(center,center);
	center += *points[0];

	normalize(center);
	normalize(squared_radius);

	return Sphere(center, squared_radius);
}

template<class NumberType_>
typename Kernel<NumberType_>::NumberType
Kernel<NumberType_>::
circumradius(const PointContainer& points) const
{
	return circumsphere(points).squared_radius();
}

template<class NumberType_>
typename Kernel<NumberType_>::Point
Kernel<NumberType_>::
circumcenter(const PointContainer& points) const
{
	return circumsphere(points).center();
}

template<class NumberType_>
typename Kernel<NumberType_>::NumberType
Kernel<NumberType_>::
side_of_circumsphere(const PointContainer& points, const Point& p) const
{
	Sphere s = circumsphere(points);
	return s.side_of(p);
}

/* --- Kernel Private --- */
template<class NumberType_>
typename Kernel<NumberType_>::NumberType&
Kernel<NumberType_>::
normalize(NumberType& n) const
{
	return number_traits<typename Point::NumberType>::normalize(n);
}

template<class NumberType_>
typename Kernel<NumberType_>::Point&
Kernel<NumberType_>::
normalize(Point& p) const
{
	for (unsigned int i = 0; i < p.size(); ++i)
		number_traits<typename Point::NumberType>::normalize(p(i));
	return p;
}
