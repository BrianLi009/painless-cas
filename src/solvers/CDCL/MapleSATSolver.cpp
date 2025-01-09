// MapleSAT includes
#include <mapleSAT/mapleSAT/core/Dimacs.h>
#include <mapleSAT/mapleSAT/simp/SimpSolver.h>
#include <mapleSAT/mapleSAT/utils/System.h>
#include <mapleSAT/mapleSAT/mtl/Vec.h>

#include "utils/Logger.hpp"
#include "utils/Parameters.hpp"
#include "utils/System.hpp"

#include "solvers/CDCL/MapleSATSolver.hpp"

#include "containers/ClauseDatabaseSingleBuffer.hpp"

#include <iomanip>

using namespace MapleSAT;

// Macros for minisat literal representation conversion
#define MINI_LIT(lit) lit > 0 ? mkLit(lit - 1, false) : mkLit((-lit) - 1, true)

#define INT_LIT(lit) sign(lit) ? -(var(lit) + 1) : (var(lit) + 1)

static void
makeMiniVec(ClauseExchangePtr cls, vec<Lit>& mcls)
{
	for (size_t i = 0; i < cls->size; i++) {
		mcls.push(MINI_LIT(cls->lits[i]));
	}
}

void
cbkMapleSATExportClause(void* issuer, unsigned int lbd, vec<Lit>& cls)
{
	MapleSATSolver* mp = (MapleSATSolver*)issuer;

	ClauseExchangePtr ncls = ClauseExchange::create(cls.size(), lbd, mp->getSharingId());

	for (int i = 0; i < cls.size(); i++) {
		ncls->lits[i] = INT_LIT(cls[i]);
	}
	/* filtering defined by a sharing strategy */
	mp->exportClause(ncls);
}

Lit
cbkMapleSATImportUnit(void* issuer)
{
	MapleSATSolver* mp = (MapleSATSolver*)issuer;

	Lit l = lit_Undef;

	ClauseExchangePtr cls;

	if (mp->unitsToImport->getOneClause(cls) == false)
		return l;

	l = MINI_LIT(cls->lits[0]);

	return l;
}

bool
cbkMapleSATImportClause(void* issuer, unsigned int* lbd, vec<Lit>& mcls)
{
	MapleSATSolver* mp = (MapleSATSolver*)issuer;

	ClauseExchangePtr cls;

	if (mp->m_clausesToImport->getOneClause(cls) == false)
	{
		mp->m_clausesToImport->shrinkDatabase();
		return false;
	}
	makeMiniVec(cls, mcls);

	*lbd = cls->lbd;

	return true;
}

MapleSATSolver::MapleSATSolver(int id, const std::shared_ptr<ClauseDatabase>& clauseDB)
	: SolverCdclInterface(id, clauseDB, SolverCdclType::MAPLESAT)
	, clausesToAdd(__globalParameters__.defaultClauseBufferSize)
{
	solver = new SimpSolver();

	this->unitsToImport = std::make_unique<ClauseDatabaseSingleBuffer>(__globalParameters__.defaultClauseBufferSize);

	solver->cbkExportClause = cbkMapleSATExportClause;
	solver->cbkImportClause = cbkMapleSATImportClause;
	solver->cbkImportUnit = cbkMapleSATImportUnit;
	solver->issuer = this;

	initializeTypeId<MapleSATSolver>();
}

MapleSATSolver::MapleSATSolver(const MapleSATSolver& other,
									 int id,
									 const std::shared_ptr<ClauseDatabase>& clauseDB)
	: SolverCdclInterface(id, clauseDB, SolverCdclType::MAPLESAT)
	, clausesToAdd(__globalParameters__.defaultClauseBufferSize)
{
	solver = new SimpSolver(*(other.solver));

	this->unitsToImport = std::make_unique<ClauseDatabaseSingleBuffer>(__globalParameters__.defaultClauseBufferSize);

	solver->cbkExportClause = cbkMapleSATExportClause;
	solver->cbkImportClause = cbkMapleSATImportClause;
	solver->cbkImportUnit = cbkMapleSATImportUnit;
	solver->issuer = this;

	initializeTypeId<MapleSATSolver>();
}

MapleSATSolver::~MapleSATSolver()
{
	delete solver;
}

// void MapleSATSolver::loadFormula(const char *filename)
// {
//    gzFile in = gzopen(filename, "rb");

//    parse_DIMACS(in, *solver);

//    gzclose(in);

//    return true;
// }

// Get the number of variables of the formula
unsigned int
MapleSATSolver::getVariablesCount()
{
	return solver->nVars();
}

// Get a variable suitable for search splitting
int
MapleSATSolver::getDivisionVariable()
{
	return (rand() % getVariablesCount()) + 1;
}

// Set initial phase for a given variable
void
MapleSATSolver::setPhase(const unsigned int var, const bool phase)
{
	solver->setPolarity(var - 1, phase ? true : false);
}

// Bump activity for a given variable
void
MapleSATSolver::bumpVariableActivity(const int var, const int times)
{
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void
MapleSATSolver::setSolverInterrupt()
{
	stopSolver = true;

	solver->interrupt();

	LOGDEBUG1("Asking MapleSAT (%d, %u) to end", this->getSolverId(), this->getSolverTypeId());
}

void
MapleSATSolver::unsetSolverInterrupt()
{
	stopSolver = false;

	solver->clearInterrupt();
}

// Diversify the solver
void
MapleSATSolver::diversify(const SeedGenerator& getSeed)
{
	/* TODO enhance this diversification */

	int seed = this->getSolverTypeId();
	if (seed == ID_XOR) {
		solver->GE = true;
	} else {
		solver->GE = false;
	}

	if (seed % 2) {
		solver->VSIDS = false;
	} else {
		solver->VSIDS = true;
	}

	if (seed % 4 >= 2) {
		solver->verso = false;
	} else {
		solver->verso = true;
	}
}

// Solve the formula with a given set of assumptions
// return 10 for SAT, 20 for UNSAT, 0 for UNKNOWN
SatResult
MapleSATSolver::solve(const std::vector<int>& cube)
{
	unsetSolverInterrupt();

	std::vector<ClauseExchangePtr> tmp;

	tmp.clear();
	clausesToAdd.getClauses(tmp);

	for (size_t ind = 0; ind < tmp.size(); ind++) {
		vec<Lit> mcls;
		makeMiniVec(tmp[ind], mcls);

		if (solver->addClause(mcls) == false) {
			LOG0(" unsat when adding cls");
			return SatResult::UNSAT;
		}
	}

	vec<Lit> miniAssumptions;
	for (size_t ind = 0; ind < cube.size(); ind++) {
		miniAssumptions.push(MINI_LIT(cube[ind]));
	}

	lbool res = solver->solveLimited(miniAssumptions);

	if (res == l_True)
		return SatResult::SAT;

	if (res == l_False)
		return SatResult::UNSAT;

	return SatResult::UNKNOWN;
}

void
MapleSATSolver::loadFormula(const char* filename)
{
	gzFile in = gzopen(filename, "rb");
	parse_DIMACS(in, *solver);
	gzclose(in);
}

void
MapleSATSolver::addClause(ClauseExchangePtr clause)
{
	clausesToAdd.addClause(clause);

	setSolverInterrupt();
}

bool
MapleSATSolver::importClause(const ClauseExchangePtr& clause)
{
	assert(clause->size > 0);

	if (clause->size == 1) {
		unitsToImport->addClause(clause);
	} else {
		m_clausesToImport->addClause(clause);
	}
	return true;
}

void
MapleSATSolver::addClauses(const std::vector<ClauseExchangePtr>& clauses)
{
	clausesToAdd.addClauses(clauses);

	setSolverInterrupt();
}

void
MapleSATSolver::addInitialClauses(const std::vector<simpleClause>& clauses, unsigned int nbVars)
{
	for (size_t ind = 0; ind < clauses.size(); ind++) {
		vec<Lit> mcls;

		for (size_t i = 0; i < clauses[ind].size(); i++) {
			int lit = clauses[ind][i];
			int var = abs(lit);

			while (solver->nVars() < var) {
				solver->newVar();
			}

			mcls.push(MINI_LIT(lit));
		}

		if (solver->addClause(mcls) == false) {
			LOG0(" unsat when adding initial cls");
			return;
		}
	}

	LOG2("The Maple Solver %d loaded all the clauses", this->getSolverId());
}

void
MapleSATSolver::importClauses(const std::vector<ClauseExchangePtr>& clauses)
{
	for (auto cls : clauses) {
		importClause(cls);
	}
}

void
MapleSATSolver::printStatistics()
{
	SolvingCdclStatistics stats;

	stats.conflicts = solver->conflicts;
	stats.propagations = solver->propagations;
	stats.restarts = solver->starts;
	stats.decisions = solver->decisions;

	std::cout << std::left << std::setw(15) << ("| MC" + std::to_string(this->getSolverTypeId())) << std::setw(20)
			  << ("| " + std::to_string(stats.conflicts)) << std::setw(20)
			  << ("| " + std::to_string(stats.propagations)) << std::setw(17) << ("| " + std::to_string(stats.restarts))
			  << std::setw(20) << ("| " + std::to_string(stats.decisions)) << std::setw(20) << "|"
			  << "\n";
}

void
MapleSATSolver::printWinningLog()
{
	this->SolverCdclInterface::printWinningLog();
	LOGSTAT("The winner is MapleSAT(%d, %u) ", this->getSolverId(), this->getSolverTypeId());
}

std::vector<int>
MapleSATSolver::getModel()
{
	std::vector<int> model;

	for (int i = 0; i < solver->nVars(); i++) {
		if (solver->model[i] != l_Undef) {
			int lit = solver->model[i] == l_True ? i + 1 : -(i + 1);
			model.push_back(lit);
		}
	}

	return model;
}

std::vector<int>
MapleSATSolver::getFinalAnalysis()
{
	std::vector<int> outCls;

	for (int i = 0; i < solver->conflict.size(); i++) {
		outCls.push_back(INT_LIT(solver->conflict[i]));
	}

	return outCls;
}

std::vector<int>
MapleSATSolver::getSatAssumptions()
{
	std::vector<int> outCls;
	vec<Lit> lits;
	solver->getAssumptions(lits);
	for (int i = 0; i < lits.size(); i++) {
		outCls.push_back(-(INT_LIT(lits[i])));
	}
	return outCls;
};

void
MapleSATSolver::setStrengthening(bool b)
{
	solver->setStrengthening(b);
}