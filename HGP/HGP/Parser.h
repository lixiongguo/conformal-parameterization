#pragma once

#include "CGAL/CGAL_Mesh.h"
#include <vector>

class Vec2 {
public:
    double v[2];

    Vec2() {}
    Vec2(double x, double y);
    ~Vec2() {}

    double& operator[](int index);
    double operator[](int index) const;
};

class Vec3 {
public:
    double v[3];

    Vec3() {}
    Vec3(double x, double y, double z);
    ~Vec3() {}

    double& operator[](int index);
    double operator[](int index) const;
};

class Cone {
public:
    int posidx;
    double x;

    Cone() {}
    Cone(int p, double xx);
    ~Cone() {}
};

class Seam {
public:
    int faceidx, posidx, rot;

    Seam() {}
    Seam(int f, int p, int r);
    ~Seam() {}
};

class MeshBuffer {
public:
    MeshBuffer();
    ~MeshBuffer();

    std::vector<int> idx_sizes, sample_sizes, sample_row_sizes;
    std::vector<int> idx_uv, idx_pos, idx_nor;
    std::vector<Vec3> positions, normals;
    std::vector<Vec2> uvs;
    std::vector<Seam> seams;
    std::vector<Cone> cones;
};

class Parser {
public:
    Parser() {}
    ~Parser() {}

    static bool loadOBJ(const char* filename, MeshBuffer* mb, MeshBuffer* smb);
    static bool loadVectorField(const char* filename, Mesh& cgalMesh);
    static void setMeshAdditionalData(Mesh& cgalMesh, MeshBuffer& m1);
};
