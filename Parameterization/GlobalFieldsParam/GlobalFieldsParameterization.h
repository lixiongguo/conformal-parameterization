#ifndef GLOBAL_FIELDS_PARAMETERIZATION_H
#define GLOBAL_FIELDS_PARAMETERIZATION_H

#include "Parameterization.h"
#include <Eigen/Core>

/**
 * Extended base class for global field-aligned parameterization algorithms
 * that require a 4-RoSy cross / frame field as input.
 *
 * The cross field is a per-face unit direction (nF × 3) that represents
 * one of the four rotationally-equivalent directions of a 4-RoSy field.
 *
 * Typical pipeline:
 *   1. Construct the subclass (QuadCover, MIQQuad, PGP, etc.)
 *   2. Call setCrossField(dirs) with a per-face direction field, or let the
 *      subclass estimate it automatically (e.g. from principal curvature).
 *   3. Call parameterize().
 *
 * Algorithms supported:
 *   - QuadCover  (Kälberer et al. 2007)  — branch cover + global param
 *   - MIQQuad    (Bommes et al. 2009)   — mixed-integer quadrangulation
 *   - PGP        (Ray et al. 2006)      — periodic global parameterization
 */
class GlobalFieldsParameterization : public Parameterization {
public:
    // constructor
    GlobalFieldsParameterization(Mesh& mesh0);

    // destructor
    virtual ~GlobalFieldsParameterization() {}

    /**
     * Set an external per-face cross / direction field.
     * @param dirs  (nFaces × 3) matrix, one unit-length direction per face.
     *              If no external field is provided, the subclass will
     *              estimate one automatically (e.g. from principal curvature).
     */
    virtual void setCrossField(const Eigen::MatrixXd& dirs);

    /** Get the current per-face direction field. */
    const Eigen::MatrixXd& getCrossField() const;

    /** Whether an external cross field has been provided. */
    bool hasExternalField() const;

protected:
    // per-face direction field (nF × 3), one unit vector per face
    // derived classes should populate this before parameterize()
    Eigen::MatrixXd faceDirs_;

    // true if setCrossField() was called with valid data
    bool hasExternalField_;
};

#endif
