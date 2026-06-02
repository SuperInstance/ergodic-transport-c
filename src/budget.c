#include "ergodic_transport.h"
#include <math.h>
#include <float.h>

double predict_consumption(const double pi[], const BudgetState *bs) {
    double total = 0.0;
    for (int i = 0; i < bs->n; i++)
        total += pi[i] * bs->resources[i];
    return total;
}

double budget_adequacy(const BudgetPlan *plan, const double pi[], const BudgetState *bs) {
    double budget_total = 0.0;
    for (int i = 0; i < plan->n; i++)
        budget_total += plan->allocation[i];

    double predicted = predict_consumption(pi, bs);
    return budget_total - predicted;
}

double budget_safety_margin(const TransitionMatrix *tm, const BudgetState *bs,
                            double epsilon) {
    /*
     * Safety margin for non-stationary transients.
     *
     * During the mixing period (first t_mix steps), the chain hasn't converged
     * to its stationary distribution. The consumption can deviate from the
     * long-run prediction. This function estimates how much extra budget to set
     * aside.
     *
     * Derivation:
     *   Let C_t = cost at step t (random, depends on chain state).
     *   The total cost during mixing is S = sum_{t=0}^{t_mix-1} C_t.
     *   E[S] = t_mix * mean_cost  (approximately, since chain isn't stationary yet).
     *   Var(S) ≈ t_mix * var_cost  (assuming weak dependence; exact bound would
     *          use the spectral gap, but this is a conservative approximation).
     *
     *   A concentration bound (Chebyshev) gives:
     *     P(S - E[S] >= k * sqrt(t_mix * var_cost)) <= 1/k^2
     *
     *   We pick k = sqrt(t_mix) so the safety margin is:
     *     margin = mean_cost * t_mix + sqrt(t_mix) * sqrt(t_mix * var_cost)
     *            = mean_cost * t_mix + t_mix * std_cost
     *            = t_mix * (mean_cost + std_cost)
     *
     *   But the caller already accounts for the stationary cost, so the EXTRA
     *   margin beyond stationary prediction is:
     *     extra = t_mix * (mean_cost + n_sigma * std_cost) - t_mix * mean_cost
     *           = t_mix * n_sigma * std_cost
     *
     *   where n_sigma = sqrt(t_mix / t_mix) = 1... Actually let's keep it simpler.
     *   We return the margin for the mixing period beyond the stationary prediction:
     *     margin = t_mix * std_cost
     *
     *   This has correct dimensions: [steps] * [cost/state] = [cost].
     */
    int mixing = markov_mixing_time(tm, epsilon, 10000);
    if (mixing < 0) return -1.0;

    /* Compute mean and std of cost under stationary distribution */
    double pi[ET_MAX_STATES];
    if (markov_stationary_distribution(tm, pi) != 0) return -1.0;

    double mean_cost = 0.0;
    for (int i = 0; i < bs->n; i++)
        mean_cost += pi[i] * bs->resources[i];

    double var_cost = 0.0;
    for (int i = 0; i < bs->n; i++)
        var_cost += pi[i] * (bs->resources[i] - mean_cost) * (bs->resources[i] - mean_cost);

    double std_cost = sqrt(var_cost);

    /* Extra budget for mixing period = mixing_steps * std_cost
       Dimensional analysis: [steps] * [cost] = [cost] ✓ */
    return (double)mixing * std_cost;
}
