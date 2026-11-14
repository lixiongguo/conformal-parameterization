#include "GlobalFieldsParameterization.h"

GlobalFieldsParameterization::GlobalFieldsParameterization(Mesh& mesh0)
    : Parameterization(mesh0)
    , hasExternalField_(false)
{
}

void GlobalFieldsParameterization::setCrossField(const Eigen::MatrixXd& dirs)
{
    faceDirs_        = dirs;
    hasExternalField_ = (dirs.cols() == 3 && dirs.rows() > 0);
}

const Eigen::MatrixXd& GlobalFieldsParameterization::getCrossField() const
{
    return faceDirs_;
}

bool GlobalFieldsParameterization::hasExternalField() const
{
    return hasExternalField_;
}
