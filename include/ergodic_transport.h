#ifndef ERGODIC_TRANSPORT_H
#define ERGODIC_TRANSPORT_H

#include <stddef.h>
#include <stdbool.h>

/* Maximum number of states */
#define ET_MAX_STATES 64

/* ========== Markov Chain ========== */

typedef struct {
    double P[ET_MAX_STATES][ET_MAX_STATES]; /* transition matrix */
    int n;                                   /* number of states */
} TransitionMatrix;

/* Build a transition matrix from flat row-major array */
void markov_init(TransitionMatrix *tm, int n, const double *flat);

/* Find stationary distribution pi via power iteration (Perron-Frobenius).
   Returns 0 on success, -1 if no convergence. */
int markov_stationary_distribution(const TransitionMatrix *tm, double pi[ET_MAX_STATES]);

/* Mixing time: steps to reach ||P^k(i,.) - pi||_1 < epsilon */
int markov_mixing_time(const TransitionMatrix *tm, double epsilon, int max_iter);

/* Seed the internal RNG. Call once before simulate if you want a specific seed.
   If never called, the default seed is 42. Does NOT reset on each simulate call. */
void markov_seed(unsigned int seed);

/* Simulate N steps from state s0. trajectory[i] = state at step i.
   Each call produces a different trajectory (RNG state persists across calls). */
void markov_simulate(const TransitionMatrix *tm, int s0, int N, int trajectory[]);

/* ========== Ergodic Theory ========== */

/* Time average: (1/N) sum_{t=0}^{N-1} f(trajectory[t]) */
double ergodic_time_average(const int trajectory[], int N, const double f[], int n_states);

/* Ensemble average: sum_i pi_i * f(i) */
double ergodic_ensemble_average(const double pi[], const double f[], int n);

/* Check ergodicity: |time_avg - ensemble_avg| < tolerance */
bool ergodic_check(const int trajectory[], int N, const double f[],
                   const double pi[], int n_states, double tolerance);

/* Birkhoff bound: how many samples needed for accuracy eps with confidence.
   Uses Chebyshev-like bound: N >= var_f / (eps^2 * delta). */
int birkhoff_bound(const double pi[], const double f[], int n,
                   double eps, double delta);

/* Empirical occupation measure from trajectory. mu[i] = fraction of time in state i */
void occupation_measure(const int trajectory[], int N, double mu[], int n_states);

/* ========== Budget Planning ========== */

typedef struct {
    double resources[ET_MAX_STATES]; /* resource level at each state */
    int n;                            /* number of states */
} BudgetState;

typedef struct {
    double allocation[ET_MAX_STATES]; /* budget allocation per state */
    int n;
} BudgetPlan;

/* Predict long-run resource consumption using ergodic theorem */
double predict_consumption(const double pi[], const BudgetState *bs);

/* Check if budget covers predicted consumption. Returns margin (positive = safe). */
double budget_adequacy(const BudgetPlan *plan, const double pi[], const BudgetState *bs);

/* Extra budget needed for non-stationary periods.
   Uses: mean_cost + n_sigma * std_cost, where n_sigma = sqrt(mixing_time)
   to account for concentration during the mixing period.
   Returns the additional budget beyond the stationary prediction. */
double budget_safety_margin(const TransitionMatrix *tm, const BudgetState *bs,
                            double epsilon);

/* ========== Ergodicity Checks ========== */

/* Check if a Markov chain is ergodic (irreducible + aperiodic).
   Returns true if ergodic, false otherwise.
   Optionally fills reason[256] with a human-readable explanation. */
bool et_is_ergodic(const TransitionMatrix *tm, char reason[256]);

/* ========== Wasserstein Distance & Control ========== */

/* Wasserstein-1 distance between two distributions on 1D line (states 0..n-1).
   Uses the CDF formula: W1 = sum |F(x) - G(x)| */
double wasserstein_1_distance(const double mu[], const double nu[], int n);

/* Drift: W1 distance between current empirical measure and target (stationary). */
double consumption_drift(const int trajectory[], int N,
                         const double pi[], int n_states);

/* Compute control correction: how much to adjust budget at each state */
void control_correction(const double current_mu[], const double target_pi[],
                        const double cost[],
                        double correction[], int n);

/* Simple LQR step: linear-quadratic regulator for budget adjustment.
   K = gain, returns adjusted allocation. */
void lqr_step(const double allocation[], const double target[],
              const double cost[], double gain,
              double result[], int n);

#endif /* ERGODIC_TRANSPORT_H */
