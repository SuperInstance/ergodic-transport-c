#include "ergodic_transport.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>

/* Deterministic RNG for reproducibility */
static double rng_state = 42.0;

static double rng_next(void) {
    rng_state = fmod(rng_state * 1103515245.0 + 12345.0, 2147483647.0);
    if (rng_state < 0) rng_state += 2147483647.0;
    return rng_state / 2147483647.0;
}

void markov_seed(unsigned int seed) {
    rng_state = (double)(seed % 2147483647u);
    if (rng_state == 0.0) rng_state = 1.0;
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
    /* NOTE: rng_state persists across calls — do NOT reset here.
       Call markov_seed() before the first simulate call if you need
       a specific starting seed. */

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

/* ========== Ergodicity Check ========== */

/* GCD helper */
static int gcd(int a, int b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) { int t = b; b = a % b; a = t; }
    return a;
}

/* BFS reachability from a given state. visited[i] = 1 if reachable. */
static void bfs_reachable(const TransitionMatrix *tm, int start, bool visited[ET_MAX_STATES]) {
    int queue[ET_MAX_STATES];
    int head = 0, tail = 0;
    memset(visited, 0, ET_MAX_STATES * sizeof(bool));
    visited[start] = true;
    queue[tail++] = start;
    while (head < tail) {
        int s = queue[head++];
        for (int j = 0; j < tm->n; j++) {
            if (tm->P[s][j] > 0.0 && !visited[j]) {
                visited[j] = true;
                queue[tail++] = j;
            }
        }
    }
}

/* Compute period of state s in the chain. Period = GCD of lengths of all cycles through s.
   Uses BFS to find shortest path distances from s back to s modulo various lengths. */
static int state_period(const TransitionMatrix *tm, int s) {
    int n = tm->n;
    /* dist[i] = shortest distance from s to i, or -1 if unreachable */
    int dist[ET_MAX_STATES];
    memset(dist, -1, sizeof(int) * n);
    dist[s] = 0;

    /* BFS */
    int queue[ET_MAX_STATES * ET_MAX_STATES];
    int head = 0, tail = 0;
    queue[tail++] = s;

    int g = 0; /* GCD accumulator */

    while (head < tail) {
        int u = queue[head++];
        for (int j = 0; j < n; j++) {
            if (tm->P[u][j] > 0.0) {
                if (dist[j] < 0) {
                    dist[j] = dist[u] + 1;
                    queue[tail++] = j;
                } else {
                    /* Found a cycle through s: s -> ... -> u -> j -> ... -> s
                       The cycle length through this edge is dist[u] + 1 + (dist back from j to s).
                       For simplicity, any return path gives cycle length = dist[u] + 1 + dist[j]
                       (since dist[j] is the shortest s->j, and we found another path to j). */
                    int cycle_len = dist[u] + 1 + dist[j];
                    g = gcd(g, cycle_len);
                    if (g == 1) return 1; /* early exit: aperiodic */
                }
            }
        }
    }

    return g == 0 ? 1 : g;
}

bool et_is_ergodic(const TransitionMatrix *tm, char reason[256]) {
    int n = tm->n;

    if (n <= 0) {
        if (reason) snprintf(reason, 256, "chain has no states");
        return false;
    }

    /* Check irreducibility: from state 0, all states must be reachable,
       and from every state, state 0 must be reachable. */
    bool reachable_from_0[ET_MAX_STATES];
    bfs_reachable(tm, 0, reachable_from_0);

    for (int i = 0; i < n; i++) {
        if (!reachable_from_0[i]) {
            if (reason) snprintf(reason, 256,
                "not irreducible: state %d not reachable from state 0", i);
            return false;
        }
    }

    /* Also check reverse reachability (state 0 reachable from all) */
    /* Build transpose graph and BFS from 0 */
    bool can_reach_0[ET_MAX_STATES];
    /* Use transpose BFS */
    int queue[ET_MAX_STATES];
    int head = 0, tail = 0;
    memset(can_reach_0, 0, ET_MAX_STATES * sizeof(bool));
    can_reach_0[0] = true;
    queue[tail++] = 0;
    while (head < tail) {
        int s = queue[head++];
        for (int i = 0; i < n; i++) {
            if (tm->P[i][s] > 0.0 && !can_reach_0[i]) {
                can_reach_0[i] = true;
                queue[tail++] = i;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        if (!can_reach_0[i]) {
            if (reason) snprintf(reason, 256,
                "not irreducible: state 0 not reachable from state %d", i);
            return false;
        }
    }

    /* Check aperiodicity: period of state 0 (all states in same class have same period) */
    int period = state_period(tm, 0);
    if (period != 1) {
        if (reason) snprintf(reason, 256,
            "chain has period %d (not aperiodic)", period);
        return false;
    }

    if (reason) snprintf(reason, 256, "chain is ergodic (irreducible + aperiodic)");
    return true;
}
