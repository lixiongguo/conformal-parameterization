#ifndef CONE_PARAMETERIZATION_H
#define CONE_PARAMETERIZATION_H

#include "Parameterization.h"
#include <unordered_map>
#include <vector>

/**
 * Extended base class for parameterization algorithms that support
 * cone singularities (target angle sum at selected vertices).
 *
 * Cone singularity: at vertex v, the sum of incident angles is set to a
 * prescribed value Theta(v) instead of the default 2π (interior) or π (boundary).
 *
 * Common values:
 *   Theta = 4π  → cone with deficit -2π (saddle-like)
 *   Theta =  π  → cone with deficit  π (boundary-like interior)
 *   Theta = 6π  → cone with deficit -4π, etc.
 *
 * Usage:
 *   1. Construct the subclass.
 *   2. Call setConeSingulars(idxVec, angleVec) before parameterize().
 *   3. Call parameterize().
 */
class ConeParameterization : public Parameterization {
public:
    // constructor
    ConeParameterization(Mesh& mesh0);

    // destructor
    virtual ~ConeParameterization() {}

    /**
     * Set cone singularities.
     * @param coneIdx    vertex indices of cone points
     * @param coneAngles target angle sums at those vertices (e.g. 4π, 6π)
     *
     * Override this if additional processing is needed (e.g. RicciFlow
     * converts angles to target curvatures).
     */
    virtual void setConeSingulars(const std::vector<int>& coneIdx,
                                  const std::vector<double>& coneAngles);

    /** Check whether a vertex is designated as a cone singularity. */
    bool isConeVertex(int vertexIndex) const;

    /** Get the target angle sum for a cone vertex.
     *  Returns 0 if the vertex is not a cone. */
    double getConeTargetAngle(int vertexIndex) const;

    /** Get read-only access to the full cone singularities map. */
    const std::unordered_map<int, double>& getConeSingulars() const;

    /** Clear all cone singularities. */
    void clearConeSingulars();

protected:
    // vertex_index → target_angle_sum (e.g. 4π, 6π)
    std::unordered_map<int, double> coneSingulars;
};

#endif
