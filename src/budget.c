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
    /* Extra budget = mixing_time * max deviation from predicted consumption */
    int mixing = markov_mixing_time(tm, epsilon, 10000);
    if (mixing < 0) return -1.0;

    double max_cost = 0.0;
    for (int i = 0; i < bs->n; i++)
        if (bs->resources[i] > max_cost) max_cost = bs->resources[i];

    /* During mixing period, worst case is spending max_cost every step */
    return mixing * max_cost * epsilon;
}
