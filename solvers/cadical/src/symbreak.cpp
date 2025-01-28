#include "symbreak.hpp"
#include <iostream>
#include "hash_values.h"
#include <iomanip>
#include <algorithm>

FILE * canonicaloutfile = NULL;
FILE * noncanonicaloutfile = NULL;
FILE * exhaustfile = NULL;
FILE * musoutfile = NULL;

// The kth entry estimates the number of permuations needed to show canonicity in order (k+1)
long perm_cutoff[MAXORDER] = {0, 0, 0, 0, 0, 0, 20, 50, 125, 313, 783, 1958, 4895, 12238, 30595, 76488, 191220, 478050, 1195125, 2987813, 7469533, 18673833, 46684583};
long canon = 0;
long noncanon = 0;
double canontime = 0;
double noncanontime = 0;
long canonarr[MAXORDER] = {};
long noncanonarr[MAXORDER] = {};
double canontimearr[MAXORDER] = {};
double noncanontimearr[MAXORDER] = {};
#ifdef PERM_STATS
long canon_np[MAXORDER] = {};
long noncanon_np[MAXORDER] = {};
#endif
long muscount = 0;
long muscounts[17] = {};
double mustime = 0;

const double TIME_PRINT_INTERVAL = 10;  // Print every 10 seconds
double last_time_checkpoint = 0;

std::vector<int> permutation_counts[MAXORDER];

SymmetryBreaker::SymmetryBreaker(CaDiCaL::Solver * s, int order) : solver(s) {
    if (order == 0) {
        std::cout << "c Need to provide order to use programmatic code" << std::endl;
        return;
    }
    
    std::cout << "c Pseudo-check is " << (solver->pseudocheck ? "enabled" : "disabled") << std::endl;
    
    n = order;
    num_edge_vars = n*(n-1)/2;
    assign = new int[num_edge_vars];
    fixed = new bool[num_edge_vars];
    colsuntouched = new int[n];
    solver->connect_external_propagator(this);
    for (int i = 0; i < num_edge_vars; i++) {
        assign[i] = l_Undef;
        fixed[i] = false;
    }
    std::cout << "c Running orderly generation on order " << n << " (" << num_edge_vars << " edge variables)" << std::endl;
    current_trail.push_back(std::vector<int>());
    for (int i = 0; i < num_edge_vars; i++) {
        solver->add_observed_var(i+1);
    }
}

SymmetryBreaker::~SymmetryBreaker () {
    if (n != 0) {
        solver->disconnect_external_propagator ();
        delete [] assign;
        delete [] colsuntouched;
        delete [] fixed;
        printf("Number of solutions   : %ld\n", sol_count);
        printf("Canonical subgraphs   : %-12" PRIu64 "   (%.0f /sec)\n", canon, canon/canontime);
        for(int i=2; i<n; i++) {
            printf("          order %2d    : %-12" PRIu64 "   (%.0f /sec)\n", i+1, canonarr[i], canonarr[i]/canontimearr[i]);
        }
        printf("Noncanonical subgraphs: %-12" PRIu64 "   (%.0f /sec)\n", noncanon, noncanon/noncanontime);
        for(int i=2; i<n; i++) {
            printf("          order %2d    : %-12" PRIu64 "   (%.0f /sec)\n", i+1, noncanonarr[i], noncanonarr[i]/noncanontimearr[i]);
        }
        printf("Canonicity checking   : %g s\n", canontime);
        printf("Noncanonicity checking: %g s\n", noncanontime);
        printf("Total canonicity time : %g s\n", canontime+noncanontime);
        
        if (solver->tracking) {
            printf("\nPermutation Statistics:\n");
            for(int i=2; i<n; i++) {
                if (subgraph_count[i] > 0) {
                    printf("Order %2d: avg perms = %.1f, max perms = %ld, total perms = %ld\n",
                        i+1,
                        (double)total_perms[i] / subgraph_count[i],
                        max_perms[i],
                        total_perms[i]);
                }
            }
            printf("\n");
        }
    }
}

void SymmetryBreaker::notify_assignment(int lit, bool is_fixed) {
    assign[abs(lit)-1] = (lit > 0 ? l_True : l_False);
    if (is_fixed) {
        fixed[abs(lit)-1] = true;
    } else {
        current_trail.back().push_back(lit);
    }
}

void SymmetryBreaker::notify_new_decision_level () {
    current_trail.push_back(std::vector<int>());
}

void SymmetryBreaker::notify_backtrack (size_t new_level) {
    while (current_trail.size() > new_level + 1) {
        for (const auto& lit: current_trail.back()) {
            const int x = abs(lit) - 1;
            // Don't remove literals that have been fixed
            if(fixed[x])
                continue;
            assign[x] = l_Undef;
            const int col = 1+(-1+sqrt(1+8*x))/2;
            for(int i=col; i<n; i++)
                colsuntouched[i] = false;
        }
        current_trail.pop_back();
    }
}

bool SymmetryBreaker::cb_check_found_model (const std::vector<int> & model) {
    assert(model.size() == num_edge_vars);
    sol_count += 1;

#ifdef VERBOSE
    std::cout << "c New solution was found: ";
#endif
    std::vector<int> clause;
    for (const auto& lit: model) {
#ifdef VERBOSE
        if (lit > 0) {
            std::cout << lit << " ";
        }
#endif
        clause.push_back(-lit);
    }
#ifdef VERBOSE
    std::cout << std::endl;
#endif
    new_clauses.push_back(clause);
    solver->add_trusted_clause(clause);
    if(solver->permoutfile != NULL) {
        fprintf(solver->permoutfile, "Complete solution\n");
    }

    return false;
}

bool SymmetryBreaker::cb_has_external_clause () {
    if(!new_clauses.empty())
        return true;

    long hash = 0;

    // Initialize i to be the first column that has been touched since the last analysis
    int i = 2;
    for(; i < n; i++) {
        if(!colsuntouched[i])
            break;
    }
    // Ensure variables are defined and update current graph hash
    for(int j = 0; j < i*(i-1)/2; j++) {
        if(assign[j] == l_Undef)
            return false;
        else if(assign[j] == l_True)
            hash += hash_values[j];
    }
    for(; i < n; i++) {
        // Ensure variables are defined and update current graph hash
        for(int j = i*(i-1)/2; j < i*(i+1)/2; j++) {
            if(assign[j]==l_Undef) {
                return false;
            }
            if(assign[j]==l_True) {
                hash += hash_values[j];
            }
        }
        colsuntouched[i] = true;

        // Check if current graph hash has been seen
        if(canonical_hashes[i].find(hash)==canonical_hashes[i].end())
        {
            // Found a new subgraph of order i+1 to test for canonicity
            const double before = CaDiCaL::absolute_process_time();
            // Run canonicity check
            int p[i+1]; // Permutation on i+1 vertices
            int x, y;   // These will be the indices of first adjacency matrix entry that demonstrates noncanonicity
            int mi;     // This will be the index of the maximum defined entry of p
            bool ret = (hash == 0) ? true : is_canonical(i+1, p, x, y, mi, solver->pseudocheck && i < n-1);

            const double after = CaDiCaL::absolute_process_time();

            // If subgraph is canonical
            if (ret) {
                canon++;
                canontime += (after-before);
                canonarr[i]++;
                canontimearr[i] += (after-before);
                canonical_hashes[i].insert(hash);

                if(canonicaloutfile != NULL) {
                    fprintf(canonicaloutfile, "a ");
                    for(int j = 0; j < i*(i+1)/2; j++)
                        fprintf(canonicaloutfile, "%s%d ", assign[j]==l_True ? "" : "-", j+1);
                    fprintf(canonicaloutfile, "0\n");
                    fflush(canonicaloutfile);
                }
            }
            // If subgraph is not canonical then block it
            else {
                noncanon++;
                noncanontime += (after-before);
                noncanonarr[i]++;
                noncanontimearr[i] += (after-before);
                new_clauses.push_back(std::vector<int>());

                if(solver->minclause) {
                    // Generate a blocking clause smaller than the naive blocking clause
                    new_clauses.back().push_back(-(x*(x-1)/2+y+1));
                    const int px = MAX(p[x], p[y]);
                    const int py = MIN(p[x], p[y]);
                    new_clauses.back().push_back(px*(px-1)/2+py+1);
                    for(int ii=0; ii < x+1; ii++) {
                        for(int jj=0; jj < ii; jj++) {
                            if(ii==x && jj==y) {
                                break;
                            }
                            const int pii = MAX(p[ii], p[jj]);
                            const int pjj = MIN(p[ii], p[jj]);
                            if(ii==pii && jj==pjj) {
                                continue;
                            } else if(assign[ii*(ii-1)/2+jj] == l_True) {
                                new_clauses.back().push_back(-(ii*(ii-1)/2+jj+1));
                            } else if (assign[pii*(pii-1)/2+pjj] == l_False) {
                                new_clauses.back().push_back(pii*(pii-1)/2+pjj+1);
                            }
                        }
                    }
                } else {
                    // Generate the naive blocking clause
                    for(int jj = 0; jj < i*(i+1)/2; jj++) {
                        new_clauses.back().push_back((jj+1) * (assign[jj]==l_True ? -1 : 1));
                    }
                }
                solver->add_trusted_clause(new_clauses.back());

                if(noncanonicaloutfile != NULL) {
                    const int size = new_clauses.back().size();
                    for(int j = 0; j < size; j++)
                        fprintf(noncanonicaloutfile, "%d ", new_clauses.back()[j]);
                    fprintf(noncanonicaloutfile, "0\n");
                    fflush(noncanonicaloutfile);
                }

                if(solver->permoutfile != NULL) {
                    for(int j = 0; j <= mi; j++)
                        fprintf(solver->permoutfile, "%s%d", j == 0 ? "" : " ", p[j]);
                    fprintf(solver->permoutfile, "\n");
                }

                return true;
            }
        }
    }

    // No programmatic clause generated
    return false;
}

int SymmetryBreaker::cb_add_external_clause_lit () {
    if (new_clauses.empty()) return 0;
    else {
        assert(!new_clauses.empty());
        size_t clause_idx = new_clauses.size() - 1;
        if (new_clauses[clause_idx].empty()) {
            new_clauses.pop_back();
            return 0;
        }

        int lit = new_clauses[clause_idx].back();
        new_clauses[clause_idx].pop_back();
        return lit;
    }
}

int SymmetryBreaker::cb_decide () { return 0; }
int SymmetryBreaker::cb_propagate () { return 0; }
int SymmetryBreaker::cb_add_reason_clause_lit (int plit) {
    (void)plit;
    return 0;
};

// Returns true when the k-vertex subgraph (with adjacency matrix M) is canonical
// M is determined by the current assignment to the first k*(k-1)/2 variables
// If M is noncanonical, then p, x, and y will be updated so that
// * p will be a permutation of the rows of M so that row[i] -> row[p[i]] generates a matrix smaller than M (and therefore demonstrates the noncanonicity of M)
// * (x,y) will be the indices of the first entry in the permutation of M demonstrating that the matrix is indeed lex-smaller than M
// * i will be the maximum defined index defined in p
bool SymmetryBreaker::is_canonical(int k, int p[], int& x, int& y, int& i, bool opt_pseudo_test = false) {
    int pl[k]; // pl[k] contains the current list of possibilities for kth vertex (encoded bitwise)
    int pn[k+1]; // pn[k] contains the initial list of possibilities for kth vertex (encoded bitwise)
    pl[0] = (1 << k) - 1;
    pn[0] = (1 << k) - 1;
    i = 0;
    int last_x = 0;
    int last_y = 0;

    int np = 1;
    int limit = INT32_MAX;

    // If pseudo-test enabled then stop test if it is taking over 10 times longer than average
    if(opt_pseudo_test && k >= 7) {
        limit = 10*perm_cutoff[k-1];
    }

    while(np < limit) {
        // If no possibilities for ith vertex then backtrack
        if(pl[i]==0) {
            // Backtrack to vertex that has at least two possibilities
            while((pl[i] & (pl[i] - 1)) == 0) {
                if(last_x > p[i]) {
                    last_x = p[i];
                    last_y = 0;
                }
                i--;
                if(i==-1) {
                    // Update tracking stats before returning
                    if (solver->tracking) {
                        total_perms[k-1] += np;
                        if (np > max_perms[k-1]) {
                            max_perms[k-1] = np;
                        }
                        subgraph_count[k-1]++;
                        permutation_counts[k-1].push_back(np);
                        
                        double current_time = CaDiCaL::absolute_process_time();
                        double time_since_checkpoint = current_time - last_time_checkpoint;
                        if (time_since_checkpoint >= TIME_PRINT_INTERVAL) {
                            print_tracking_stats();
                            last_time_checkpoint = current_time;
                        }
                    }
                    return true;
                }
            }
            // Remove p[i] as a possibility from the ith vertex
            pl[i] = pl[i] & ~(1 << p[i]);
        }

        p[i] = log2(pl[i] & -pl[i]); // Get index of rightmost high bit
        pn[i+1] = pn[i] & ~(1 << p[i]); // List of possibilities for (i+1)th vertex

        // If pseudo-test enabled then stop shortly after the first row is no longer fixed
        if(i == 0 && p[i] == 1 && opt_pseudo_test && k < n) {
            limit = np + 100;
        }

        // Check if the entry on which to begin lex-checking needs to be updated
        if(last_x > p[i]) {
            last_x = p[i];
            last_y = 0;
        }
        if(i == last_x) {
            last_y = 0;
        }

        // Determine if the permuted matrix p(M) is lex-smaller than M
        bool lex_result_unknown = false;
        x = last_x == 0 ? 1 : last_x;
        y = last_y;
        int j;
        for(j=last_x*(last_x-1)/2+last_y; j<k*(k-1)/2; j++) {
            if(x > i) {
                // Unknown if permutation produces a larger or smaller matrix
                lex_result_unknown = true;
                break;
            }
            const int px = MAX(p[x], p[y]);
            const int py = MIN(p[x], p[y]);
            const int pj = px*(px-1)/2 + py;
            if(assign[j] == l_False && assign[pj] == l_True) {
                // Permutation produces a larger matrix; stop considering
                break;
            }
            if(assign[j] == l_True && assign[pj] == l_False) {
                // Update tracking stats before returning
                if (solver->tracking) {
                    total_perms[k-1] += np;
                    if (np > max_perms[k-1]) {
                        max_perms[k-1] = np;
                    }
                    subgraph_count[k-1]++;
                    permutation_counts[k-1].push_back(np);
                    
                    double current_time = CaDiCaL::absolute_process_time();
                    double time_since_checkpoint = current_time - last_time_checkpoint;
                    if (time_since_checkpoint >= TIME_PRINT_INTERVAL) {
                        print_tracking_stats();
                        last_time_checkpoint = current_time;
                    }
                }
                return false;
            }

            y++;
            if(x==y) {
                x++;
                y = 0;
            }
        }
        last_x = x;
        last_y = y;

        if(lex_result_unknown) {
            // Lex result is unknown; need to define p[i] for another i
            i++;
            pl[i] = pn[i];
        }
        else {
            np++;  // Count this as a complete permutation check
            // Remove p[i] as a possibility from the ith vertex
            pl[i] = pl[i] & ~(1 << p[i]);
        }
    }

    // Update tracking stats before final return
    if (solver->tracking) {
        total_perms[k-1] += np;
        if (np > max_perms[k-1]) {
            max_perms[k-1] = np;
        }
        subgraph_count[k-1]++;
        permutation_counts[k-1].push_back(np);
        
        double current_time = CaDiCaL::absolute_process_time();
        double time_since_checkpoint = current_time - last_time_checkpoint;
        if (time_since_checkpoint >= TIME_PRINT_INTERVAL) {
            print_tracking_stats();
            last_time_checkpoint = current_time;
        }
    }
    return true;
}

void SymmetryBreaker::print_tracking_stats() {
    std::cout << "\nIntermediate Permutation Statistics:\n";
    for(int i=2; i<n; i++) {
        if (subgraph_count[i] > 0) {
            // Calculate median
            std::sort(permutation_counts[i].begin(), permutation_counts[i].end());
            double median = (permutation_counts[i].size() % 2 == 0) ? 
                (permutation_counts[i][permutation_counts[i].size()/2 - 1] + permutation_counts[i][permutation_counts[i].size()/2]) / 2.0 : 
                permutation_counts[i][permutation_counts[i].size()/2];

            // Calculate 25th and 75th percentiles
            double p25 = permutation_counts[i][(int)(0.25 * (permutation_counts[i].size() - 1))];
            double p75 = permutation_counts[i][(int)(0.75 * (permutation_counts[i].size() - 1))];

            std::cout << "Order " << std::setw(2) << i+1 
                     << ": avg perms = " << std::fixed << std::setprecision(1) 
                     << (double)total_perms[i] / subgraph_count[i]
                     << ", max perms = " << max_perms[i]
                     << ", median perms = " << median
                     << ", 25th percentile = " << p25
                     << ", 75th percentile = " << p75
                     << ", total perms = " << total_perms[i] << "\n";
        }
    }
    std::cout << "\n";
    std::cout.flush();
}