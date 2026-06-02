#include "ergodic_transport.h"
#include <math.h>
#include <string.h>

double ergodic_time_average(const int trajectory[], int N, const double f[], int n_states) {
    (void)n_states;
    double sum = 0.0;
    for (int t = 0; t < N; t++)
        sum += f[trajectory[t]];
    return sum / N;
}

double ergodic_ensemble_average(const double pi[], const double f[], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++)
        sum += pi[i] * f[i];
    return sum;
}

bool ergodic_check(const int trajectory[], int N, const double f[],
                   const double pi[], int n_states, double tolerance) {
    double ta = ergodic_time_average(trajectory, N, f, n_states);
    double ea = ergodic_ensemble_average(pi, f, n_states);
    return fabs(ta - ea) < tolerance;
}

int birkhoff_bound(const double pi[], const double f[], int n,
                   double eps, double delta) {
    /* Mean */
    double mu = 0.0;
    for (int i = 0; i < n; i++) mu += pi[i] * f[i];
    /* Variance */
    double var = 0.0;
    for (int i = 0; i < n; i++)
        var += pi[i] * (f[i] - mu) * (f[i] - mu);
    /* Chebyshev: P(|avg - mu| >= eps) <= var / (N * eps^2)
       Need var / (N * eps^2) <= delta => N >= var / (eps^2 * delta) */
    if (eps <= 0 || delta <= 0) return -1;
    int N = (int)ceil(var / (eps * eps * delta));
    return N < 1 ? 1 : N;
}

void occupation_measure(const int trajectory[], int N, double mu[], int n_states) {
    for (int i = 0; i < n_states; i++) mu[i] = 0.0;
    for (int t = 0; t < N; t++)
        mu[trajectory[t]] += 1.0;
    for (int i = 0; i < n_states; i++)
        mu[i] /= N;
}
