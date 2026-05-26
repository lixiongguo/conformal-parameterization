#include "MosekSolver.h"
#include <cstdio>
#include <cstring>

using namespace MosekSolver;

Solver::Solver():
    bkc(NULL), bkx(NULL),
    ptrb(NULL), ptre(NULL), asub(NULL), qsubi(NULL), qsubj(NULL),
    blc(NULL), buc(NULL), blx(NULL), bux(NULL), aval(NULL), qval(NULL), c(NULL), xx(NULL),
    qcsubi(NULL), qcsubj(NULL), qcval(NULL),
    handle(NULL),
    env(NULL),
    task(NULL),
    r(MSK_RES_OK),
    variables(0),
    constraints(0),
    numanz(0),
    numqcnz(0),
    numqnz(0)
{
    initializeEnvironment();
}

Solver::~Solver()
{
    reset();
    if (env) MSK_deleteenv(&env);
}

static void MSKAPI printError(MSKuserhandle_t /*handle*/, const char* str)
{
    std::fprintf(stderr, "%s", str);
}

void Solver::initializeEnvironment()
{
    r = MSK_makeenv(&env, NULL);
    if (r == MSK_RES_OK) {
        MSK_linkfunctoenvstream(env, MSK_STREAM_LOG, NULL, printError);
    }
}

bool Solver::createTask(MSKtask_t& outTask)
{
    r = MSK_maketask(env, constraints, variables, &outTask);
    return r == MSK_RES_OK;
}

void Solver::allocateMemory()
{
    bkc = new MSKboundkeye[constraints];
    blc = new MSKrealt[constraints];
    buc = new MSKrealt[constraints];
    bkx = new MSKboundkeye[variables];
    blx = new MSKrealt[variables];
    bux = new MSKrealt[variables];
    ptrb = new MSKint32t[variables];
    ptre = new MSKint32t[variables];
    asub = new MSKint32t[numanz];
    aval = new MSKrealt[numanz];
    qcsubi = new MSKint32t[numqcnz];
    qcsubj = new MSKint32t[numqcnz];
    qcval = new MSKrealt[numqcnz];
    qsubi = new MSKint32t[numqnz];
    qsubj = new MSKint32t[numqnz];
    qval = new MSKrealt[numqnz];
    c = new MSKrealt[variables];
    xx = new MSKrealt[variables];

    // Default: empty constraints 0 = 0 (imaginary boundary faces etc.)
    for (MSKint32t i = 0; i < constraints; ++i) {
        bkc[i] = MSK_BK_FX;
        blc[i] = 0.0;
        buc[i] = 0.0;
    }
    std::memset(ptrb, 0, sizeof(MSKint32t) * variables);
    std::memset(ptre, 0, sizeof(MSKint32t) * variables);
    std::memset(c, 0, sizeof(MSKrealt) * variables);
}

bool Solver::initialize(MSKint32t variables_, MSKint32t constraints_,
                        MSKint32t numanz_, MSKint32t numqcnz_, MSKint32t numqnz_)
{
    if (r != MSK_RES_OK) return false;

    variables = variables_;
    constraints = constraints_;
    numanz = numanz_;
    numqcnz = numqcnz_;
    numqnz = numqnz_;

    if (!createTask(task)) return false;
    allocateMemory();
    return true;
}

static void logMosekError(MSKrescodee code, const char* where)
{
    char sym[64] = {0};
    char desc[256] = {0};
    MSK_getcodedesc(code, sym, desc);
    std::fprintf(stderr, "MosekSolver [%s]: %s — %s\n", where, sym, desc);
}

bool Solver::setSolution()
{
    MSKsolstae solsta;
    r = MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
    if (r != MSK_RES_OK) {
        logMosekError(r, "getsolsta");
        return false;
    }

    switch (solsta) {
        case MSK_SOL_STA_OPTIMAL:
        case MSK_SOL_STA_PRIM_FEAS:
        case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
        case MSK_SOL_STA_UNKNOWN:
            r = MSK_getxx(task, MSK_SOL_ITR, xx);
            if (r == MSK_RES_OK) return true;
            logMosekError(r, "getxx");
            return false;
        case MSK_SOL_STA_PRIM_INFEAS_CER: {
            MSKprostae prosta;
            if (MSK_getprosta(task, MSK_SOL_ITR, &prosta) == MSK_RES_OK) {
                std::fprintf(stderr, "MosekSolver: primal infeasible (prosta=%d).\n", (int)prosta);
            } else {
                std::fprintf(stderr, "MosekSolver: primal infeasible (angle QP).\n");
            }
            MSK_writedata(task, "mosek_angle_qp.ptf");
            return false;
        }
        default:
            std::fprintf(stderr, "MosekSolver: solution status = %d\n", (int)solsta);
            return false;
    }
}

bool Solver::solve(ProblemType type)
{
    r = MSK_inputdata(task,
                      constraints, variables,
                      constraints, variables,
                      c, 0.0,
                      ptrb, ptre, asub, aval,
                      bkc, blc, buc,
                      bkx, blx, bux);
    if (r != MSK_RES_OK) {
        logMosekError(r, "inputdata");
        return false;
    }

    switch (type) {
        case QO:
            if (numqnz > 0) {
                r = MSK_putqobj(task, numqnz, qsubi, qsubj, qval);
            }
            if (r == MSK_RES_OK && numqcnz > 0) {
                r = MSK_putqconk(task, 0, numqcnz, qcsubi, qcsubj, qcval);
            }
            break;

        case GECO:
            std::fprintf(stderr,
                "MosekSolver: GECO (nonlinear) problems require MOSEK 7 API; "
                "use QO for CirclePatterns angle optimization.\n");
            return false;

        default:
            return false;
    }
    if (r != MSK_RES_OK) return false;

    MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    MSK_putintparam(task, MSK_IPAR_INFEAS_REPORT_AUTO, MSK_ON);

    r = MSK_optimize(task);
    if (r != MSK_RES_OK) {
        logMosekError(r, "optimize");
        return false;
    }

    return setSolution();
}

void Solver::deallocateMemory()
{
    delete [] bkc; delete [] blc; delete [] buc;
    delete [] bkx; delete [] blx; delete [] bux;
    delete [] ptrb; delete [] ptre; delete [] asub; delete [] aval;
    delete [] qcsubi; delete [] qcsubj; delete [] qcval;
    delete [] qsubi; delete [] qsubj; delete [] qval; delete [] c;
    delete [] xx;
    bkc = bkx = NULL;
    ptrb = ptre = asub = qcsubi = qcsubj = qsubi = qsubj = NULL;
    blc = buc = blx = bux = aval = qcval = qval = c = xx = NULL;
}

void Solver::reset()
{
    deallocateMemory();
    if (task) {
        MSK_deletetask(&task);
        task = NULL;
    }
}
