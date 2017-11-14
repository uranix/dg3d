#pragma once

#include "point.h"

#include <string>
#include <vector>

struct tet {
    int v[4];
    int label;
};

struct face {
    int v[3];
    int label;
    point n;
    int elemleft, elemright;
};

struct mesh {
    std::vector<point> pts;
    std::vector<tet> elems;
    std::vector<face> faces;
    mesh(const std::string &fn, const bool fortran_numbering = true);
};
