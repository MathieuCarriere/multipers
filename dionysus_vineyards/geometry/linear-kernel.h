#ifndef __LINEAR_KERNEL_H__
#define __LINEAR_KERNEL_H__

#include <stack>
#include <iostream>
#include <boost/operators.hpp>

template<class T>
class LinearKernel
{
	public:
        struct                      Function;

		typedef						T					                                RootType;
		typedef						std::stack<RootType>								RootStack;

		static void					solve(const Function& f, RootStack& stack)          { if (f.a1 != 0) stack.push(-f.a0/f.a1); }
		static RootType				root(const T& r)									{ return r; }
		static int					sign_at(const Function& f, const RootType& r)       { RootType y = f(r); if (y < 0) return -1; if (y > 0) return 1; return 0; }
		static RootType				between(const RootType& r1, const RootType& r2)		{ return (r1 + r2)/2; }
		static int					sign_at_negative_infinity(const Function& f)        { if (f.a1 < 0) return 1; if (f.a1 > 0) return -1; if (f.a0 > 0) return 1; if (f.a0 < 0) return -1; return 0; }
};

template<class T>
struct LinearKernel<T>::Function: boost::additive<Function>
{
        typedef                     T                                                   RootType;

                                    Function(RootType aa0, RootType aa1 = 0):
                                        a0(aa0), a1(aa1)                                {}

        RootType                    operator()(RootType r) const                        { return a1 * r + a0; }
        Function&                   operator+=(const Function& other)                   { a1 += other.a1; a0 += other.a0; return *this; }
        Function&                   operator-=(const Function& other)                   { a1 -= other.a1; a0 -= other.a0; return *this; }
        std::ostream&               operator<<(std::ostream& out) const                 { out << a1 << "*x + " << a0; return out; }

        RootType                    a0, a1;
};

// TODO: need to make this generic
std::ostream&                       operator<<(std::ostream& out, const LinearKernel<double>::Function& f)
{
    return f.operator<<(out);
}

#endif
