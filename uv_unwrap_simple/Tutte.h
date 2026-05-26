#ifndef TUTTE_H
#define TUTTE_H

#include "Parameterization.h"

enum class TutteBoundary { CIRCLE, SQUARE };

class Tutte : public Parameterization {
public:
    Tutte(Mesh& mesh0, TutteBoundary boundaryShape = TutteBoundary::CIRCLE);
    
    void parameterize() override;

private:
    // fix boundary vertices to circle or square
    void pinBoundary();
    
    // find all boundary vertices in order
    void findBoundaryLoop(std::vector<int>& boundaryVerts) const;
    
    // solve Laplacian system for interior vertices
    void solveLaplacian(const std::vector<int>& boundaryVerts);
    
    TutteBoundary m_shape;
};

#endif
