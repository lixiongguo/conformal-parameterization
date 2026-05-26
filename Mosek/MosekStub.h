/**
 * MOSEK Stub for WebAssembly compilation
 * Provides minimal type definitions to allow compilation without the MOSEK library.
 * Circle Patterns and Trust Region methods will not be available in WASM build.
 */

#ifndef MOSEK_STUB_H
#define MOSEK_STUB_H

// MOSEK type stubs
typedef int MSKint32t;
typedef double MSKrealt;
typedef int MSKboundkeye;
typedef int MSKrescodee;
typedef void* MSKenv_t;
typedef void* MSKtask_t;

// MOSEK bound key constants
#define MSK_BK_FX 0  // Fixed
#define MSK_BK_FR 1  // Free
#define MSK_BK_LO 2  // Lower bound
#define MSK_BK_RA 3  // Range
#define MSK_BK_UP 4  // Upper bound

// MOSEK infinity constant
#define MSK_INFINITY 1e30

namespace MosekSolver {
    enum ProblemType { QO, GECO };

    class Solver {
    public:
        Solver() {}
        ~Solver() {}

        bool initialize(MSKint32t variables_, MSKint32t constraints_,
                        MSKint32t numanz_, MSKint32t numqcnz_, MSKint32t numqnz_) {
            return false;  // Always fail - MOSEK not available in WASM
        }

        bool solve(ProblemType type) {
            return false;
        }

        void reset() {}

        // Stub memory allocations (no-op)
        MSKboundkeye *bkc = nullptr, *bkx = nullptr;
        MSKint32t *ptrb = nullptr, *ptre = nullptr, *asub = nullptr;
        MSKint32t *qcsubi = nullptr, *qcsubj = nullptr, *qsubi = nullptr, *qsubj = nullptr;
        MSKrealt *blc = nullptr, *buc = nullptr, *blx = nullptr, *bux = nullptr;
        MSKrealt *aval = nullptr, *qcval = nullptr, *qval = nullptr, *c = nullptr, *xx = nullptr;

    private:
        MSKenv_t env = nullptr;
        MSKtask_t task = nullptr;
        MSKrescodee r;
        MSKint32t variables = 0, constraints = 0, numanz = 0, numqcnz = 0, numqnz = 0;
    };
}

#endif // MOSEK_STUB_H
