#include "ergodic_transport.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

/* Deterministic RNG for reproducibility */
static double rng_state = 42.0;

static double rng_next(void) {
    rng_state = fmod(rng_state * 1103515245.0 + 12345.0, 2147483647.0);
    if (rng_state < 0) rng_state += 2147483647.0;
    return rng_state / 2147483647.0;
}

void markov_init(TransitionMatrix *tm, int n, const double *flat) {
    tm->n = n;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            tm->P[i][j] = flat[i * n + j];
}

int markov_stationary_distribution(const TransitionMatrix *tm, double pi[ET_MAX_STATES]) {
    int n = tm->n;
    /* Initialize uniform */
    for (int i = 0; i < n; i++) pi[i] = 1.0 / n;

    for (int iter = 0; iter < 100000; iter++) {
        double new_pi[ET_MAX_STATES];
        for (int j = 0; j < n; j++) {
            new_pi[j] = 0.0;
            for (int i = 0; i < n; i++)
                new_pi[j] += pi[i] * tm->P[i][j];
        }
        /* Normalize */
        double sum = 0.0;
        for (int j = 0; j < n; j++) sum += new_pi[j];
        for (int j = 0; j < n; j++) new_pi[j] /= sum;

        /* Check convergence */
        double diff = 0.0;
        for (int j = 0; j < n; j++)
            diff += fabs(new_pi[j] - pi[j]);
        memcpy(pi, new_pi, n * sizeof(double));
        if (diff < 1e-12) return 0;
    }
    return -1;
}

int markov_mixing_time(const TransitionMatrix *tm, double epsilon, int max_iter) {
    int n = tm->n;
    double pi[ET_MAX_STATES];
    if (markov_stationary_distribution(tm, pi) != 0) return -1;

    /* Start from worst-case (deterministic at state 0) */
    double dist[ET_MAX_STATES];
    for (int i = 0; i < n; i++) dist[i] = 0.0;
    dist[0] = 1.0;

    for (int k = 0; k < max_iter; k++) {
        /* Multiply dist * P */
        double new_dist[ET_MAX_STATES];
        for (int j = 0; j < n; j++) {
            new_dist[j] = 0.0;
            for (int i = 0; i < n; i++)
                new_dist[j] += dist[i] * tm->P[i][j];
        }
        memcpy(dist, new_dist, n * sizeof(double));

        /* Total variation distance */
        double tv = 0.0;
        for (int j = 0; j < n; j++)
            tv += fabs(dist[j] - pi[j]);
        tv *= 0.5;
        if (tv < epsilon) return k + 1;
    }
    return -1; /* didn't converge */
}

void markov_simulate(const TransitionMatrix *tm, int s0, int N, int trajectory[]) {
    int n = tm->n;
    int state = s0;
    rng_state = 42.0; /* reset for reproducibility */

    for (int t = 0; t < N; t++) {
        trajectory[t] = state;
        double r = rng_next();
        double cum = 0.0;
        int next = state; /* fallback */
        for (int j = 0; j < n; j++) {
            cum += tm->P[state][j];
            if (r <= cum) { next = j; break; }
        }
        state = next;
    }
}
