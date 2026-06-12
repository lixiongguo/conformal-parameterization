#include "ComplexPolyField.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

struct Args {
    std::string inputObj;
    std::string outputDir;
    int N = 4;
};

bool parseArgs(int argc, char** argv, Args& args) {
    if (argc < 3) return false;
    args.inputObj = argv[1];
    args.outputDir = argv[2];
    for (int i = 3; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--n" && i + 1 < argc) {
            args.N = std::max(1, std::atoi(argv[++i]));
        } else {
            std::cerr << "Unknown option: " << a << "\n";
            return false;
        }
    }
    return true;
}

void usage() {
    std::cout
        << "Usage:\n"
        << "  complex_poly_field.exe <input.obj> <output_dir> [--n <N>]\n\n"
        << "Computes N-RoSy field via complex polynomial (Diamanti et al. SGP 2014).\n";
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parseArgs(argc, argv, args)) {
        usage();
        return 1;
    }

    Mesh mesh;
    if (!mesh.read(args.inputObj)) {
        std::cerr << "Failed to read mesh: " << args.inputObj << "\n";
        return 2;
    }

    ComplexPolyField field(mesh, args.N);
    field.compute();

    std::filesystem::create_directories(args.outputDir);

    const std::string faceCsv =
        (std::filesystem::path(args.outputDir) / "face_nrosy_field.csv").string();
    const std::string singCsv =
        (std::filesystem::path(args.outputDir) / "singularity_detected.csv").string();

    {
        std::ofstream out(faceCsv);
        out << "face,theta,dir0x,dir0y,dir0z,cx,cy,cz\n";
        for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); ++f) {
            if (f->isBoundary()) continue;
            const int fi = field.localFaceIndex(f->index);
            if (fi < 0) continue;

            const Eigen::Vector3d d = field.direction(fi, 0);
            HalfEdgeCIter h = f->he;
            const Eigen::Vector3d c =
                (h->vertex->position + h->next->vertex->position +
                 h->next->next->vertex->position) /
                3.0;
            out << f->index << "," << field.theta(fi) << ","
                << d.x() << "," << d.y() << "," << d.z() << ","
                << c.x() << "," << c.y() << "," << c.z() << "\n";
        }
    }

    {
        std::ofstream out(singCsv);
        out << "vertex,index_Nrosy\n";
        for (int vi = 0; vi < field.numVertices(); ++vi) {
            const int w = field.windingIndex(vi);
            if (w == -1) continue;
            if (w != 0) {
                out << vi << "," << w << "\n";
            }
        }
    }

    std::cout << "Done.\n"
              << "  faces: " << field.numFaces() << "  N=" << field.N() << "\n"
              << "  field: " << faceCsv << "\n"
              << "  singularities: " << singCsv << "\n";
    return 0;
}
