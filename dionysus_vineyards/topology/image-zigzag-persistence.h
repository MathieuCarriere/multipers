#ifndef __IMAGE_ZIGZAG_PERSISTENCE_H__
#define __IMAGE_ZIGZAG_PERSISTENCE_H__

#include "zigzag-persistence.h"
#include <limits>

struct SimplexSubcomplexData
{
                SimplexSubcomplexData(bool sc = false):
                    subcomplex(sc)                              {}

    bool        subcomplex;
};

template<class BirthID_ = Empty<> >
class ImageZigzagPersistence: public ZigzagPersistence<BirthID_, SimplexSubcomplexData>
{
    public:
        typedef                 BirthID_                                                    BirthID;
        typedef                 ZigzagPersistence<BirthID, SimplexSubcomplexData>           Parent;

        typedef                 typename Parent::IndexDeathPair                             IndexDeathPair;
        typedef                 typename Parent::Death                                      Death;

        typedef                 typename Parent::ZIndex                                     ZIndex;
        typedef                 typename Parent::BIndex                                     BIndex;
        typedef                 typename Parent::SimplexIndex                               SimplexIndex;
        typedef                 typename Parent::ZColumn                                    ZColumn;
        typedef                 typename Parent::BColumn                                    BColumn;
        typedef                 typename Parent::BRow                                       BRow;
        typedef                 typename Parent::CRow                                       CRow;
        typedef                 typename Parent::ZNode                                      ZNode;
        typedef                 typename Parent::BNode                                      BNode;
        typedef                 typename Parent::SimplexNode                                SimplexNode;


                                ImageZigzagPersistence():
                                    im_last(Parent::z_list.end()), cok_begin(Parent::z_list.end()),
                                    im_order_begin(std::numeric_limits<int>::min()/2),
                                    cok_order_begin(std::numeric_limits<int>::max()/2)
                                {}

        IndexDeathPair          add(ZColumn         bdry,
                                    bool            subcomplex,
                                    const BirthID&  birth = BirthID())                      { ImageZZVisitor zzv(subcomplex); return Parent::add(bdry, birth, zzv); }

        Death                   remove(SimplexIndex s, const BirthID& birth = BirthID())    { ImageZZVisitor zzv(s->subcomplex); return Parent::remove(s, birth, zzv); }


        ZIndex                  image_begin()                                               { return Parent::z_list.begin(); }
        ZIndex                  image_end()                                                 { return cok_begin; }
        BIndex                  boundary_end()                                              { return Parent::b_list.end(); }


        // Class: ImageZZVisitor
        // Contains all the tweaks to the normal zigzag algorithm to make it compute image zigzag
        class ImageZZVisitor: public Parent::ZigzagVisitor
        {
            public:
                                    ImageZZVisitor(bool sc = false):
                                        subcomplex(sc), birth_in_image(false)                   {}

                // Sets the subcomplex property of the new simplex
                SimplexIndex        new_simplex(Parent& zz);

                // Decides where to put the new column (image or cokernel)
                ZIndex              new_z_in_add(Parent& zz, const ZColumn& z, const BRow& u);

                // Checks if there is a boundary entirely in the subcomplex, and sets birth_in_image accordingly
                BIndex              select_j_in_remove(Parent& zz, const CRow& c_row);

                ZIndex              new_z_in_remove(Parent& zz);

                // Updates im_last and cok_begin if necessary
                void                erasing_z(Parent& zz, ZIndex j);

                // Determines if there is a death in the image
                Death               death(Parent& zz, ZIndex dying_z);

            private:
                ZIndex              append_in_image(Parent& zz);
                ZIndex              append_in_cokernel(Parent& zz);
                ZIndex              prepend_in_image(Parent& zz);
                ZIndex              prepend_in_cokernel(Parent& zz);
                bool                in_subcomplex(const ZColumn& z);

                bool                subcomplex, birth_in_image;
        };

        const int                   im_order_begin;
        const int                   cok_order_begin;

    private:
        using                   Parent::make_remover;
        using                   Parent::make_appender;
        using                   Parent::make_adder;

    private:
        ZIndex                  im_last;                // index of the last image cycle
        ZIndex                  cok_begin;              // index of the first cokernel cycle
};


#include "image-zigzag-persistence.hpp"

#endif
