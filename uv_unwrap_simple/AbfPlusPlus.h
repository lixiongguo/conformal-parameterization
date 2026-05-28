#ifndef ABF_PLUS_PLUS_H
#define ABF_PLUS_PLUS_H

#include "Parameterization.h"

class AbfPlusPlus : public Parameterization {
public:
    explicit AbfPlusPlus(Mesh& mesh0, int maxIterations = 12);

    void parameterize() override;

private:
    int m_maxIterations;
};

#endif

