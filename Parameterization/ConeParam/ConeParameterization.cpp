#include "ConeParameterization.h"

ConeParameterization::ConeParameterization(Mesh& mesh0)
    : Parameterization(mesh0)
{
}

void ConeParameterization::setConeSingulars(const std::vector<int>& coneIdx,
                                            const std::vector<double>& coneAngles)
{
    coneSingulars.clear();
    for (size_t i = 0; i < coneIdx.size() && i < coneAngles.size(); i++) {
        coneSingulars[coneIdx[i]] = coneAngles[i];
    }
}

bool ConeParameterization::isConeVertex(int vertexIndex) const
{
    return coneSingulars.find(vertexIndex) != coneSingulars.end();
}

double ConeParameterization::getConeTargetAngle(int vertexIndex) const
{
    auto it = coneSingulars.find(vertexIndex);
    return (it != coneSingulars.end()) ? it->second : 0.0;
}

const std::unordered_map<int, double>& ConeParameterization::getConeSingulars() const
{
    return coneSingulars;
}

void ConeParameterization::clearConeSingulars()
{
    coneSingulars.clear();
}
