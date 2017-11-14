#include <libmeshb7.h>
#include <iostream>
#include <vector>
#include <fstream>

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "USAGE: ./gmf2vtk <input>.mesh[b] <output>.vtk" << std::endl;
        return 1;
    }

    int64_t MeshIndex;
    int Version, Dimension;
    MeshIndex = GmfOpenMesh(argv[1], GmfRead, &Version, &Dimension);

    std::cout << "Version = " << Version << ", Dimension = " << Dimension << std::endl;

    if (!MeshIndex || (Version != 2) || (Dimension != 3)) {
        return 2;
    }

    int nbv = GmfStatKwd(MeshIndex, GmfVertices);
    int nbf = GmfStatKwd(MeshIndex, GmfTriangles);
    int nbt = GmfStatKwd(MeshIndex, GmfTetrahedra);

    std::cout << "#V/F/T = " << nbv << " / " << nbf << " / " << nbt << std::endl;

    typedef double point[3];
    typedef int face[3];
    typedef int elem[4];

    std::vector<point> vrt(nbv);
    std::vector<face> bnd(nbf);
    std::vector<elem> tet(nbt);

    std::vector<int> bndlab(nbf);
    std::vector<int> tetlab(nbt);

    GmfGotoKwd(MeshIndex, GmfVertices);
    for (int i = 0; i < nbv; i++) {
        int vlabel;
        GmfGetLin(MeshIndex, GmfVertices, &vrt[i][0], &vrt[i][1], &vrt[i][2], &vlabel);
    }

    GmfGotoKwd(MeshIndex, GmfTriangles);
    for (int i = 0; i < nbf; i++)
        GmfGetLin(MeshIndex, GmfTriangles, &bnd[i][0], &bnd[i][1], &bnd[i][2], &bndlab[i]);

    GmfGotoKwd(MeshIndex, GmfTetrahedra);
    for (int i = 0; i < nbt; i++)
        GmfGetLin(MeshIndex, GmfTetrahedra, &tet[i][0], &tet[i][1], &tet[i][2], &tet[i][3], &tetlab[i]);

    GmfCloseMesh(MeshIndex);

    std::ofstream g(argv[2]);

    g << "# vtk DataFile Version 3.0\n";
    g << "nocomments\n";
    g << "ASCII\n";
    g << "DATASET UNSTRUCTURED_GRID\n";

    g << "POINTS " << nbv << " float\n";

    for (int i = 0; i < nbv; i++) {
        double px = vrt[i][0];
        double py = vrt[i][1];
        double pz = vrt[i][2];
        g << px << " " << py << " " << pz << "\n";
    }

    g << "CELLS " << nbt + nbf << " " << 5*nbt + 4*nbf << "\n";

    for (int i = 0; i < nbt; i++) {
        int v1 = tet[i][0] - 1;
        int v2 = tet[i][1] - 1;
        int v3 = tet[i][2] - 1;
        int v4 = tet[i][3] - 1;
        g << "4 " << v1 << " " << v2 << " " << v3 << " " << v4 << "\n";
    }

    for (int i = 0; i < nbf; i++) {
        int v1 = bnd[i][0] - 1;
        int v2 = bnd[i][1] - 1;
        int v3 = bnd[i][2] - 1;
        g << "3 " << v1 << " " << v2 << " " << v3 << "\n";
    }

    g << "CELL_TYPES " << nbt + nbf << "\n";
    for (int i = 0; i < nbt; i++)
        g << "10\n";
    for (int i = 0; i < nbf; i++)
        g << "5\n";

    g << "POINT_DATA " << nbv << "\n";
    g << "CELL_DATA " << nbt + nbf << "\n";

    g << "SCALARS mat int\nLOOKUP_TABLE default\n";
    for (int i = 0; i < nbt; i++) {
        g << tetlab[i] << "\n";
    }
    for (int i = 0; i < nbf; i++) {
        g << bndlab[i] << "\n";
    }

    return 0;
}
