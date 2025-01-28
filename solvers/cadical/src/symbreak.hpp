#include "internal.hpp"
#include <set>
#include <unordered_set>

#define l_False 0
#define l_True 1
#define l_Undef 2

#define MAX(X,Y) ((X) > (Y)) ? (X) : (Y)
#define MIN(X,Y) ((X) > (Y)) ? (Y) : (X)

const int MAXORDER = 40;

// Add hash values array declaration
extern const unsigned long hash_values[];

class SymmetryBreaker : CaDiCaL::ExternalPropagator {
    CaDiCaL::Solver * solver;
    std::vector<std::vector<int>> new_clauses;
    std::deque<std::vector<int>> current_trail;
    int * assign;
    bool * fixed;
    int * colsuntouched;
    int n = 0;
    long sol_count = 0;
    int num_edge_vars = 0;
    std::unordered_set<unsigned long> canonical_hashes[MAXORDER];
    std::unordered_set<unsigned long> solution_hashes;
    long total_perms[MAXORDER] = {};     
    long subgraph_count[MAXORDER] = {};  

public:
    // Only declare the constructor and destructor
    SymmetryBreaker(CaDiCaL::Solver * s, int order);
    ~SymmetryBreaker();
    void notify_assignment(int lit, bool is_fixed);
    void notify_new_decision_level ();
    void notify_backtrack (size_t new_level);
    bool cb_check_found_model (const std::vector<int> & model);
    bool cb_has_external_clause ();
    int cb_add_external_clause_lit ();
    int cb_decide ();
    int cb_propagate ();
    int cb_add_reason_clause_lit (int plit);
    bool is_canonical(int k, int p[], int& x, int& y, int& i);
    void print_tracking_stats();
};