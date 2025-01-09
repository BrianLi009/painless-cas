#pragma once

#include "containers/ClauseBuffer.hpp"
#include "containers/ClauseDatabase.hpp"
#include "utils/Threading.hpp"

#include "SolverCdclInterface.hpp"

#define MAPLESAT_

#define ID_SYM 0
#define ID_XOR 1

struct parameter
{
	int tier1;
	int chrono;
	int stable;
	int walkinitially;
	int target;
	int phase;
	int heuristic;
	int margin;
	int ccanr;
	int targetinc;
};

// Some forward declatarations for MapleSAT
namespace MapleSAT {
class SimpSolver;
class Lit;
template<class T>
class vec;
}

/// @ingroup solving_cdcl
/// Instance of a MapleSAT solver
class MapleSATSolver : public SolverCdclInterface
{
  public:
	/// Get the number of variables of the current resolution.
	unsigned int getVariablesCount();

	/// Get a variable suitable for search splitting.
	int getDivisionVariable();

	/// Set initial phase for a given variable.
	void setPhase(const unsigned int var, const bool phase);

	/// Bump activity of a given variable.
	void bumpVariableActivity(const int var, const int times);

	/// Interrupt resolution, solving cannot continue until interrupt is unset.
	void setSolverInterrupt();

	/// Remove the SAT solving interrupt request.
	void unsetSolverInterrupt();

	/// Solve the formula with a given cube.
	SatResult solve(const std::vector<int>& cube);

	/// Add a permanent clause to the formula.
	void addClause(ClauseExchangePtr clause);

	/// Add a list of permanent clauses to the formula.
	void addClauses(const std::vector<ClauseExchangePtr>& clauses);

	/// Add a list of initial clauses to the formula.
	void addInitialClauses(const std::vector<simpleClause>& clauses, unsigned int nbVars) override;

	/// Load formula from a given dimacs file, return false if failed.
	void loadFormula(const char* filename) override;

	/// Add a learned clause to the formula.
	bool importClause(const ClauseExchangePtr& clause);

	/// Add a list of learned clauses to the formula.
	void importClauses(const std::vector<ClauseExchangePtr>& clauses);

	/// Get solver statistics.
	void printStatistics();

	void printWinningLog() override;

	/// Return the model in case of SAT result.
	std::vector<int> getModel();

	/// Native diversification.
	/// @param noise added for some radomness
	/// @param family to enable some parameters dependent on the family
	void diversify(const SeedGenerator& getSeed);

	/// Constructor.
	MapleSATSolver(int id, const std::shared_ptr<ClauseDatabase>& clauseDB);

	/// Copy constructor.
	MapleSATSolver(const MapleSATSolver& other, int id, const std::shared_ptr<ClauseDatabase>& clauseDB);
	/// Destructor.
	virtual ~MapleSATSolver();

	std::vector<int> getFinalAnalysis();

	std::vector<int> getSatAssumptions();

	void setStrengthening(bool b);

	void setParameter(parameter p) {};

	void initshuffle(int id) {};

  protected:
	/// Pointer to a MapleSAT solver.
	MapleSAT::SimpSolver* solver;

	/// Buffer used to import clauses (units included).

	std::unique_ptr<ClauseDatabase> unitsToImport;

	/// Buffer used to add permanent clauses.
	ClauseBuffer clausesToAdd;

	/// Used to stop or continue the resolution.
	std::atomic<bool> stopSolver;

	/// Callback to export/import clauses.
	friend MapleSAT::Lit cbkMapleSATImportUnit(void*);
	friend bool cbkMapleSATImportClause(void*, unsigned int*, MapleSAT::vec<MapleSAT::Lit>&);
	friend void cbkMapleSATExportClause(void*, unsigned int, MapleSAT::vec<MapleSAT::Lit>&);

  public:
	;
};
