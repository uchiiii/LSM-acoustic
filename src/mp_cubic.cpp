#include "mp_cubic.h"

mpCubic::mpCubic():ext(0){}

mpCubic::mpCubic(mp4Vector v0, mp4Vector v1, mp4Vector v2, mp4Vector v3, mp4Vector v4, mp4Vector v5, mp4Vector v6, mp4Vector v7, int ext):ext(ext) {
    vertices[0] = v0;
    vertices[1] = v1;
    vertices[2] = v2;
    vertices[3] = v3;
    vertices[4] = v4;
    vertices[5] = v5;
    vertices[6] = v6;
    vertices[7] = v7;
}

mpCubic::mpCubic(std::vector<mp4Vector> v, int ext): ext(ext){
    int i = 0;
    for(auto e: v) {
        vertices[i] = e;
        i++;
    }
}

mpCubic::mpCubic(mp4Vector *v, int ext):ext(ext) {
    for(int i=0;i<=7; i++) {
        vertices[i] = v[i];
    }
}