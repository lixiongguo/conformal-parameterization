#ifndef TUTTE_H
#define TUTTE_H

#include "Parameterization.h"

enum class TutteBoundary { CIRCLE, SQUARE };
enum class TutteWeight { COTAN, UNIFORM };

class Tutte : public Parameterization {
public:
    Tutte(Mesh& mesh0,
          TutteBoundary boundaryShape = TutteBoundary::CIRCLE,
          TutteWeight weight = TutteWeight::COTAN);
    
    void parameterize() override;

private:
    // fix boundary vertices to circle or square
    void pinBoundary();
    
    // find all boundary vertices in order
    void findBoundaryLoop(std::vector<int>& boundaryVerts) const;
    
    // solve Laplacian system for interior vertices (cotan or uniform weights)
    void solveLaplacian(const std::vector<int>& boundaryVerts);
    
    TutteBoundary m_shape;
    TutteWeight m_weight;
};

#endif
