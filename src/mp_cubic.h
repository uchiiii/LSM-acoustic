#ifndef MPCUBIC_H
#define MPCUBIC_H

#include "mp_vector.h"
#include<vector>

class mpCubic {
    public:
    mp4Vector vertices[8];
    int ext;

    mpCubic();
    mpCubic(mp4Vector v0, mp4Vector v1, mp4Vector v2, mp4Vector v3, mp4Vector v4, mp4Vector v5, mp4Vector v6, mp4Vector v7, int ext);
    mpCubic(std::vector<mp4Vector> v, int ext);
    mpCubic(mp4Vector *v, int ext);
    mpCubic(const mpCubic& c);
};

#endif