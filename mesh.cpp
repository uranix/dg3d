#include "mesh.h"

#include <libmeshb7.h>

#include <unordered_map>
#include <cassert>
#include <stdexcept>

// Unique face identifier
struct faceid {
    int v[3];
    faceid(const faceid &) = default;
    faceid(int a, int b, int c) {
        int w[6];
        w[0] = w[3] = a;
        w[1] = w[4] = b;
        w[2] = w[5] = c;
        int imin;
        if (a < b)
            if (a < c)
                imin = 0;
            else
                imin = 2;
        else
            if (b < c)
                imin = 1;
            else
                imin = 2;
        for (int j = 0; j < 3; j++)
            v[j] = w[imin + j];
    }
    // Is true exactly for one of {face, flipped face}
    bool ordered() const {
        return v[1] < v[2];
    }
    const faceid flip() const {
        return faceid(v[0], v[2], v[1]);
    }
    bool operator==(const faceid &o) const {
        for (int j = 0; j < 3; j++)
            if (v[j] != o.v[j])
                return false;
        return true;
    }
};

namespace std {
    template<>
    class hash<faceid> {
        public:
        size_t operator()(const faceid &fid) const {
            size_t seed = 3;
            for (int j = 0; j < 3; j++) {
                seed ^= 0x9e3779b9 + (seed << 6) + (seed >> 2) + std::hash<int>()(fid.v[j]);
            }
            return seed;
        }
    };
}

template <class Cont, typename Key, typename Val>
void checked_insert(Cont &cont, const Key &key, const Val &v) {
    if (cont.count(key) > 0)
        throw std::logic_error("Duplicate face. Tet orientation may be inconsistent");
    cont[key] = v;
}

mesh::mesh(const std::string &fn, const bool fortran_numbering) {
    int64_t MeshIndex;
    int Version, Dimension;
    MeshIndex = GmfOpenMesh(fn.c_str(), GmfRead, &Version, &Dimension);

    if (!MeshIndex || (Version != 2) || (Dimension != 3))
        throw std::invalid_argument("Incompatible mesh");

    int nV = GmfStatKwd(MeshIndex, GmfVertices);

    pts.resize(nV);
    GmfGotoKwd(MeshIndex, GmfVertices);
    for (auto &p : pts) {
        int vlabel;
        GmfGetLin(MeshIndex, GmfVertices, &p.x, &p.y, &p.z, &vlabel);
    }

    struct faceprop {
        int elem, label;
    };
    std::unordered_map<faceid, faceprop> f2e;

    int nT = GmfStatKwd(MeshIndex, GmfTetrahedra);
    elems.resize(nT);
    GmfGotoKwd(MeshIndex, GmfTetrahedra);
    for (size_t i = 0; i < elems.size(); i++) {
        tet &t = elems[i];

        int v1, v2, v3, v4, lab;
        GmfGetLin(MeshIndex, GmfTetrahedra,
                &v1, &v2, &v3, &v4, &lab);

        v1--; v2--; v3--; v4--;

        t.v[0] = v1;
        t.v[1] = v2;
        t.v[2] = v3;
        t.v[3] = v4;

        int elem = i;
        checked_insert(f2e, faceid(v1, v2, v3), faceprop{elem, -1});
        checked_insert(f2e, faceid(v2, v1, v4), faceprop{elem, -1});
        checked_insert(f2e, faceid(v3, v4, v1), faceprop{elem, -1});
        checked_insert(f2e, faceid(v4, v3, v2), faceprop{elem, -1});

        t.label = lab;
    }

    int nB = GmfStatKwd(MeshIndex, GmfTriangles);
    GmfGotoKwd(MeshIndex, GmfTriangles);
    for (int i = 0; i < nB; i++) {
        int v1, v2, v3, lab;
        GmfGetLin(MeshIndex, GmfTriangles,
                &v1, &v2, &v3, &lab);

        v1--; v2--; v3--;

        int pos = f2e.count(faceid(v1, v2, v3));
        int neg = f2e.count(faceid(v1, v3, v2));

        if (pos == 0 && neg == 1) {
            // This is an outer boundary
            checked_insert(f2e, faceid(v1, v2, v3), faceprop{-1, lab});
            f2e[faceid(v1, v3, v2)].label = lab;
            continue;
        }

        if (pos == 1 && neg == 1) {
            // This is a slit boundary
            f2e[faceid(v1, v2, v3)].label = lab;
            f2e[faceid(v1, v3, v2)].label = lab;
            continue;
        }

        throw std::logic_error("Tet faces and boundary faces orientation mismatch, pos = " +
                std::to_string(pos)+ ", neg = " + std::to_string(neg) + ", lab = " +
                std::to_string(lab));
    }

    for (const auto &kv : f2e) {
        const faceid &id = kv.first;
        const faceprop &prop = kv.second;
        if (id.ordered()) {
            faces.push_back(face());
            face &f = faces.back();
            for (int i = 0; i < 3; i++)
                f.v[i] = id.v[i];
            f.label = prop.label;
            f.elemleft = prop.elem;

            const faceprop &flip = f2e.at(id.flip());
            if (f.label != flip.label)
                throw std::logic_error("Label of flipped face mismatch face label");
            f.elemright = flip.elem;
            // Might be true for some boundary conditions, but not possible for this mesh format
            if (f.elemleft == f.elemright)
                throw std::logic_error("Elements on the left and on the right are the same");

            // TODO: Compute normal
        }
    }

}
