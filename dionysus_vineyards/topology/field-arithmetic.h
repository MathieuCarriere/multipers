#ifndef __FIELD_ARITHMETIC_H__
#define __FIELD_ARITHMETIC_H__

#include <vector>

class ZpField
{
    public:
        typedef     int                                             Element;

                    ZpField(Element p = 2);

        Element     id()  const                                     { return 1; }
        Element     zero()  const                                   { return 0; }
        Element     init(int a) const                               { return (a % p_ + p_) % p_; }

        Element     neg(Element a) const                            { return p_ - a; }
        Element     add(Element a, Element b) const                 { return (a+b) % p_; }

        Element     inv(Element a) const                            { return inverses_[a]; }
        Element     mul(Element a, Element b) const                 { return (a*b) % p_; }
        Element     div(Element a, Element b) const                 { return mul(a, inv(b)); }

        bool        is_zero(Element a) const                        { return (a % p_) == 0; }

    private:
        Element                 p_;
        std::vector<Element>    inverses_;
};

ZpField::
ZpField(Element p):
    p_(p), inverses_(p_)
{
    for (Element i = 1; i < p_; ++i)
        for (Element j = 1; j < p_; ++j)
            if (mul(i,j) == 1)
            {
                inverses_[i] = j;
                break;
            }
}

#if 0                   // unused example; commented out to get rid of the artificial dependence on GMP
#include <gmpxx.h>

class QField
{
    public:
        typedef     mpq_class                                       Element;

                    QField()                                        {}

        Element     id()  const                                     { return 1; }
        Element     zero()  const                                   { return 0; }

        Element     neg(Element a) const                            { return -a; }
        Element     add(Element a, Element b) const                 { return (a+b); }

        Element     inv(Element a) const                            { return id()/a; }
        Element     mul(Element a, Element b) const                 { return (a*b); }
        Element     div(Element a, Element b) const                 { return a/b; }

        bool        is_zero(Element a) const                        { return a == 0; }
};
#endif
 
#endif // __FIELD_ARITHMETIC_H__
