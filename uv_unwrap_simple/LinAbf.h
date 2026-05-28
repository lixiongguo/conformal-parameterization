#ifndef LIN_ABF_H
#define LIN_ABF_H

#include "Parameterization.h"

class LinAbf : public Parameterization {
public:
    explicit LinAbf(Mesh& mesh0, int maxIterations = 6);

    void parameterize() override;

private:
    int m_maxIterations;
};

#endif

