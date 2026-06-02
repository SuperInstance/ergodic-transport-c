#include "ergodic_transport.h"
#include <math.h>
#include <string.h>

double wasserstein_1_distance(const double mu[], const double nu[], int n) {
    /* W1 on ordered states 0..n-1 via CDF: sum of |F_mu(i) - F_nu(i)| */
    double dist = 0.0;
    double cdf_mu = 0.0, cdf_nu = 0.0;
    for (int i = 0; i < n; i++) {
        cdf_mu += mu[i];
        cdf_nu += nu[i];
        dist += fabs(cdf_mu - cdf_nu);
    }
    return dist;
}

double consumption_drift(const int trajectory[], int N,
                         const double pi[], int n_states) {
    double mu[ET_MAX_STATES];
    occupation_measure(trajectory, N, mu, n_states);
    return wasserstein_1_distance(mu, pi, n_states);
}

void control_correction(const double current_mu[], const double target_pi[],
                        const double cost[],
                        double correction[], int n) {
    /* Correction proportional to deviation from target, weighted by cost */
    for (int i = 0; i < n; i++) {
        double deviation = current_mu[i] - target_pi[i];
        correction[i] = -deviation * cost[i];
    }
}

void lqr_step(const double allocation[], const double target[],
              const double cost[], double gain,
              double result[], int n) {
    /* Simple proportional controller: u = -K * (x - target)
       Adjust allocation toward target proportional to cost deviation */
    for (int i = 0; i < n; i++) {
        double error = allocation[i] - target[i] * cost[i];
        result[i] = allocation[i] - gain * error;
        if (result[i] < 0) result[i] = 0.0;
    }
}
