#include "mesh.h"

#include <iostream>

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "USAGE: " << argv[0] << " <mesh.meshb>" << std::endl;
        return 1;
    }
    mesh m(argv[1]);
    return 0;
}
