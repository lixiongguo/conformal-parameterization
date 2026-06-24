#ifndef TUTTE_H
#define TUTTE_H

#include "Parameterization.h"

enum class TutteBoundary { CIRCLE, SQUARE, FREE };
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
    
    // pin two vertices for free-boundary harmonic map (removes rigid DOFs)
    void pinTwoVertices(std::vector<int>& pinVerts);
    
    // find all boundary vertices in order
    void findBoundaryLoop(std::vector<int>& boundaryVerts) const;
    
    // solve Laplacian for vertices not in fixedVerts (cotan or uniform weights)
    void solveLaplacian(const std::vector<int>& fixedVerts);
    
    TutteBoundary m_shape;
    TutteWeight m_weight;
};

#endif
