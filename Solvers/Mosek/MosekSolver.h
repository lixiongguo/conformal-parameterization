#ifndef MOSEK_SOLVER_H
#define MOSEK_SOLVER_H

#ifndef MSKCONST
#define MSKCONST const
#endif

#include "mosek.h"
#include <functional>

namespace MosekSolver {
    struct MeshHandle {
        typedef std::function<void(double&, const double*)> ComputeEnergy;
        typedef std::function<void(double*, const double*)> ComputeGradient;
        typedef std::function<void(double*, const double*)> ComputeHessian;
        typedef std::function<void(int*)> BuildGradientSparsity;
        typedef std::function<void(int*, int*)> BuildHessianSparsity;

        MeshHandle(int numGrad_, int numHess_):
            numGrad(numGrad_),
            numHess(numHess_)
        {}

        int numGrad;
        int numHess;
        ComputeEnergy computeEnergy;
        ComputeGradient computeGradient;
        ComputeHessian computeHessian;
        BuildGradientSparsity buildGradientSparsity;
        BuildHessianSparsity buildHessianSparsity;
    };

    enum ProblemType {
        QO,   // quadratic optimization
        GECO  // general convex (requires legacy MSK_putnlfunc, not in MOSEK 11+)
    };

    class Solver {
    public:
        Solver();
        ~Solver();

        bool initialize(MSKint32t variables_, MSKint32t constraints_,
                        MSKint32t numanz_, MSKint32t numqcnz_, MSKint32t numqnz_);

        bool solve(ProblemType type);
        void reset();

        MSKboundkeye *bkc, *bkx;
        MSKint32t *ptrb, *ptre, *asub, *qcsubi, *qcsubj, *qsubi, *qsubj;
        MSKrealt *blc, *buc, *blx, *bux, *aval, *qcval, *qval, *c, *xx;
        MeshHandle *handle;

    private:
        void initializeEnvironment();
        bool createTask(MSKtask_t& task);
        void allocateMemory();
        bool setSolution();
        void deallocateMemory();

        MSKenv_t env;
        MSKtask_t task;
        MSKrescodee r;
        MSKint32t variables, constraints, numanz, numqcnz, numqnz;
    };
}

#endif
