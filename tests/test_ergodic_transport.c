#include "ergodic_transport.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TOLERANCE 1e-6
#define LOOSE_TOL 0.05

#define TEST(name) do { printf("  %-55s", #name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (!(cond)) { FAIL(msg); return; } else { PASS(); } } while(0)

/* ===== Helper: 2-state symmetric chain ===== */
static void make_symmetric_2state(TransitionMatrix *tm) {
    double P[] = {0.7, 0.3,
                  0.3, 0.7};
    markov_init(tm, 2, P);
}

/* ===== MARKOV TESTS ===== */

static void test_stationary_symmetric(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    double pi[ET_MAX_STATES];
    int rc = markov_stationary_distribution(&tm, pi);
    TEST("stationary_symmetric");
    CHECK(rc == 0 && fabs(pi[0] - 0.5) < TOLERANCE && fabs(pi[1] - 0.5) < TOLERANCE,
          "expected [0.5, 0.5]");
}

static void test_stationary_asymmetric(void) {
    TransitionMatrix tm;
    double P[] = {0.9, 0.1,
                  0.2, 0.8};
    markov_init(&tm, 2, P);
    double pi[ET_MAX_STATES];
    int rc = markov_stationary_distribution(&tm, pi);
    TEST("stationary_asymmetric");
    /* pi = [2/3, 1/3] */
    CHECK(rc == 0 && fabs(pi[0] - 2.0/3.0) < TOLERANCE && fabs(pi[1] - 1.0/3.0) < TOLERANCE,
          "expected [2/3, 1/3]");
}

static void test_stationary_3state(void) {
    TransitionMatrix tm;
    double P[] = {0.5, 0.3, 0.2,
                  0.1, 0.6, 0.3,
                  0.2, 0.2, 0.6};
    markov_init(&tm, 3, P);
    double pi[ET_MAX_STATES];
    int rc = markov_stationary_distribution(&tm, pi);
    TEST("stationary_3state");
    /* Verify pi * P = pi and sums to 1 */
    double sum = 0;
    int ok = (rc == 0);
    for (int i = 0; i < 3; i++) sum += pi[i];
    ok = ok && fabs(sum - 1.0) < TOLERANCE;
    for (int j = 0; j < 3; j++) {
        double ppi = 0;
        for (int i = 0; i < 3; i++) ppi += pi[i] * P[i*3+j];
        ok = ok && fabs(ppi - pi[j]) < 1e-8;
    }
    CHECK(ok, "pi*P != pi or sum != 1");
}

static void test_stationary_single_state(void) {
    TransitionMatrix tm;
    double P[] = {1.0};
    markov_init(&tm, 1, P);
    double pi[ET_MAX_STATES];
    int rc = markov_stationary_distribution(&tm, pi);
    TEST("stationary_single_state");
    CHECK(rc == 0 && fabs(pi[0] - 1.0) < TOLERANCE, "expected [1.0]");
}

static void test_stationary_absorbing(void) {
    TransitionMatrix tm;
    double P[] = {1.0, 0.0,
                  0.5, 0.5};
    markov_init(&tm, 2, P);
    double pi[ET_MAX_STATES];
    int rc = markov_stationary_distribution(&tm, pi);
    TEST("stationary_absorbing");
    CHECK(rc == 0 && fabs(pi[0] - 1.0) < TOLERANCE && fabs(pi[1]) < TOLERANCE,
          "expected [1.0, 0.0]");
}

static void test_mixing_time_symmetric(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    int mt = markov_mixing_time(&tm, 0.01, 10000);
    TEST("mixing_time_symmetric");
    CHECK(mt > 0, "mixing time should be positive");
}

static void test_mixing_time_single_state(void) {
    TransitionMatrix tm;
    double P[] = {1.0};
    markov_init(&tm, 1, P);
    int mt = markov_mixing_time(&tm, 0.01, 10000);
    TEST("mixing_time_single_state");
    CHECK(mt >= 0, "single state should have instant mixing");
}

static void test_simulate_length(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    int traj[100];
    markov_simulate(&tm, 0, 100, traj);
    TEST("simulate_length");
    int ok = 1;
    for (int t = 0; t < 100; t++)
        if (traj[t] < 0 || traj[t] >= 2) ok = 0;
    CHECK(ok, "all states in range");
}

static void test_simulate_start_state(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    int traj[1];
    markov_simulate(&tm, 1, 1, traj);
    TEST("simulate_start_state");
    CHECK(traj[0] == 1, "first state should be s0");
}

/* ===== ERGODIC TESTS ===== */

static void test_time_average_simple(void) {
    int traj[] = {0, 1, 0, 1, 0, 1, 0, 1};
    double f[] = {1.0, 2.0};
    double avg = ergodic_time_average(traj, 8, f, 2);
    TEST("time_average_simple");
    CHECK(fabs(avg - 1.5) < TOLERANCE, "expected 1.5");
}

static void test_ensemble_average_symmetric(void) {
    double pi[] = {0.5, 0.5};
    double f[] = {1.0, 3.0};
    double avg = ergodic_ensemble_average(pi, f, 2);
    TEST("ensemble_average_symmetric");
    CHECK(fabs(avg - 2.0) < TOLERANCE, "expected 2.0");
}

static void test_ensemble_average_single(void) {
    double pi[] = {1.0};
    double f[] = {5.0};
    double avg = ergodic_ensemble_average(pi, f, 1);
    TEST("ensemble_average_single");
    CHECK(fabs(avg - 5.0) < TOLERANCE, "expected 5.0");
}

static void test_ergodic_check_convergence(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    int traj[10000];
    markov_simulate(&tm, 0, 10000, traj);
    double pi[ET_MAX_STATES];
    markov_stationary_distribution(&tm, pi);
    double f[] = {1.0, 2.0};
    bool ok = ergodic_check(traj, 10000, f, pi, 2, 0.05);
    TEST("ergodic_check_convergence");
    CHECK(ok, "time avg should converge to ensemble avg");
}

static void test_ergodic_check_constant(void) {
    /* Constant function always passes */
    int traj[] = {0, 1, 1, 0};
    double pi[] = {0.5, 0.5};
    double f[] = {3.0, 3.0};
    bool ok = ergodic_check(traj, 4, f, pi, 2, 0.001);
    TEST("ergodic_check_constant_function");
    CHECK(ok, "constant function should always pass");
}

static void test_birkhoff_bound_positive(void) {
    double pi[] = {0.5, 0.5};
    double f[] = {0.0, 1.0};
    int N = birkhoff_bound(pi, f, 2, 0.1, 0.05);
    TEST("birkhoff_bound_positive");
    CHECK(N > 0, "bound should be positive");
}

static void test_birkhoff_bound_tighter_needs_more(void) {
    double pi[] = {0.5, 0.5};
    double f[] = {0.0, 1.0};
    int N1 = birkhoff_bound(pi, f, 2, 0.1, 0.05);
    int N2 = birkhoff_bound(pi, f, 2, 0.01, 0.05);
    TEST("birkhoff_bound_tighter_needs_more");
    CHECK(N2 > N1, "smaller epsilon needs more samples");
}

static void test_occupation_measure_uniform(void) {
    int traj[] = {0, 1, 0, 1, 0, 1, 0, 1};
    double mu[ET_MAX_STATES];
    occupation_measure(traj, 8, mu, 2);
    TEST("occupation_measure_uniform");
    CHECK(fabs(mu[0] - 0.5) < TOLERANCE && fabs(mu[1] - 0.5) < TOLERANCE,
          "expected [0.5, 0.5]");
}

static void test_occupation_measure_biased(void) {
    int traj[] = {0, 0, 0, 1};
    double mu[ET_MAX_STATES];
    occupation_measure(traj, 4, mu, 2);
    TEST("occupation_measure_biased");
    CHECK(fabs(mu[0] - 0.75) < TOLERANCE && fabs(mu[1] - 0.25) < TOLERANCE,
          "expected [0.75, 0.25]");
}

static void test_occupation_measure_single(void) {
    int traj[] = {0, 0, 0};
    double mu[ET_MAX_STATES];
    occupation_measure(traj, 3, mu, 1);
    TEST("occupation_measure_single_state");
    CHECK(fabs(mu[0] - 1.0) < TOLERANCE, "expected [1.0]");
}

/* ===== BUDGET TESTS ===== */

static void test_predict_consumption(void) {
    double pi[] = {0.5, 0.5};
    BudgetState bs = {.resources = {10.0, 20.0}, .n = 2};
    double pred = predict_consumption(pi, &bs);
    TEST("predict_consumption");
    CHECK(fabs(pred - 15.0) < TOLERANCE, "expected 15.0");
}

static void test_predict_consumption_skewed(void) {
    double pi[] = {0.75, 0.25};
    BudgetState bs = {.resources = {10.0, 30.0}, .n = 2};
    double pred = predict_consumption(pi, &bs);
    TEST("predict_consumption_skewed");
    CHECK(fabs(pred - 15.0) < TOLERANCE, "expected 15.0");
}

static void test_budget_adequacy_sufficient(void) {
    double pi[] = {0.5, 0.5};
    BudgetState bs = {.resources = {10.0, 20.0}, .n = 2};
    BudgetPlan plan = {.allocation = {8.0, 12.0}, .n = 2}; /* total = 20, pred = 15 */
    double margin = budget_adequacy(&plan, pi, &bs);
    TEST("budget_adequacy_sufficient");
    CHECK(margin > 0, "should have positive margin");
}

static void test_budget_adequacy_insufficient(void) {
    double pi[] = {0.5, 0.5};
    BudgetState bs = {.resources = {10.0, 20.0}, .n = 2};
    BudgetPlan plan = {.allocation = {3.0, 3.0}, .n = 2}; /* total = 6, pred = 15 */
    double margin = budget_adequacy(&plan, pi, &bs);
    TEST("budget_adequacy_insufficient");
    CHECK(margin < 0, "should have negative margin");
}

static void test_budget_adequacy_exact(void) {
    double pi[] = {0.5, 0.5};
    BudgetState bs = {.resources = {10.0, 20.0}, .n = 2};
    BudgetPlan plan = {.allocation = {5.0, 10.0}, .n = 2}; /* total = 15, pred = 15 */
    double margin = budget_adequacy(&plan, pi, &bs);
    TEST("budget_adequacy_exact");
    CHECK(fabs(margin) < TOLERANCE, "should be zero margin");
}

static void test_budget_safety_margin(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    BudgetState bs = {.resources = {10.0, 20.0}, .n = 2};
    double margin = budget_safety_margin(&tm, &bs, 0.01);
    TEST("budget_safety_margin");
    CHECK(margin >= 0, "safety margin should be non-negative");
}

/* ===== WASSERSTEIN TESTS ===== */

static void test_wasserstein_same_dist(void) {
    double mu[] = {0.5, 0.5};
    double nu[] = {0.5, 0.5};
    double d = wasserstein_1_distance(mu, nu, 2);
    TEST("wasserstein_same_distribution");
    CHECK(fabs(d) < TOLERANCE, "distance should be zero");
}

static void test_wasserstein_different(void) {
    double mu[] = {1.0, 0.0};
    double nu[] = {0.0, 1.0};
    double d = wasserstein_1_distance(mu, nu, 2);
    TEST("wasserstein_shifted");
    CHECK(fabs(d - 1.0) < TOLERANCE, "expected distance 1.0");
}

static void test_wasserstein_3state(void) {
    double mu[] = {1.0, 0.0, 0.0};
    double nu[] = {0.0, 0.0, 1.0};
    double d = wasserstein_1_distance(mu, nu, 3);
    TEST("wasserstein_3state_shift");
    CHECK(fabs(d - 2.0) < TOLERANCE, "expected distance 2.0");
}

static void test_wasserstein_symmetric(void) {
    double mu[] = {0.8, 0.2};
    double nu[] = {0.3, 0.7};
    double d1 = wasserstein_1_distance(mu, nu, 2);
    double d2 = wasserstein_1_distance(nu, mu, 2);
    TEST("wasserstein_symmetry");
    CHECK(fabs(d1 - d2) < TOLERANCE, "W1 should be symmetric");
}

static void test_consumption_drift_positive(void) {
    int traj[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
    double pi[] = {0.5, 0.5};
    double drift = consumption_drift(traj, 10, pi, 2);
    TEST("consumption_drift_positive");
    CHECK(drift > 0, "drift should be positive for biased trajectory");
}

static void test_consumption_drift_zero(void) {
    int traj[] = {0, 1, 0, 1, 0, 1, 0, 1};
    double pi[] = {0.5, 0.5};
    double drift = consumption_drift(traj, 8, pi, 2);
    TEST("consumption_drift_zero");
    CHECK(fabs(drift) < TOLERANCE, "drift should be zero for matching dist");
}

/* ===== CONTROL TESTS ===== */

static void test_control_correction_direction(void) {
    double mu[] = {0.7, 0.3};
    double pi[] = {0.5, 0.5};
    double cost[] = {1.0, 1.0};
    double corr[2];
    control_correction(mu, pi, cost, corr, 2);
    TEST("control_correction_direction");
    /* state 0 over-represented, correction should be negative */
    CHECK(corr[0] < 0 && corr[1] > 0, "correction should reduce over-representation");
}

static void test_control_correction_zero(void) {
    double mu[] = {0.5, 0.5};
    double pi[] = {0.5, 0.5};
    double cost[] = {1.0, 1.0};
    double corr[2];
    control_correction(mu, pi, cost, corr, 2);
    TEST("control_correction_zero_drift");
    CHECK(fabs(corr[0]) < TOLERANCE && fabs(corr[1]) < TOLERANCE,
          "no drift => no correction");
}

static void test_lqr_step_basic(void) {
    double alloc[] = {10.0, 20.0};
    double target[] = {0.5, 0.5};
    double cost[] = {10.0, 20.0};
    double result[2];
    lqr_step(alloc, target, cost, 0.5, result, 2);
    TEST("lqr_step_basic");
    /* With gain 0.5, allocation should move toward target*cost */
    CHECK(result[0] >= 0 && result[1] >= 0, "allocations should be non-negative");
}

static void test_lqr_step_at_target(void) {
    double alloc[] = {5.0, 10.0};
    double target[] = {0.5, 0.5};
    double cost[] = {10.0, 20.0};
    double result[2];
    lqr_step(alloc, target, cost, 0.5, result, 2);
    TEST("lqr_step_at_target");
    /* alloc = target * cost => no error => no change */
    CHECK(fabs(result[0] - 5.0) < TOLERANCE && fabs(result[1] - 10.0) < TOLERANCE,
          "at target should not change");
}

static void test_lqr_step_converges(void) {
    double alloc[] = {10.0, 20.0};
    double target[] = {0.5, 0.5};
    double cost[] = {10.0, 20.0};
    /* Apply LQR multiple times */
    for (int i = 0; i < 50; i++) {
        double result[2];
        lqr_step(alloc, target, cost, 0.3, result, 2);
        alloc[0] = result[0];
        alloc[1] = result[1];
    }
    TEST("lqr_step_converges");
    CHECK(fabs(alloc[0] - 5.0) < 0.1 && fabs(alloc[1] - 10.0) < 0.1,
          "should converge to target*cost");
}

/* ===== INTEGRATION TESTS ===== */

static void test_full_pipeline(void) {
    TransitionMatrix tm;
    double P[] = {0.8, 0.15, 0.05,
                  0.1, 0.7, 0.2,
                  0.1, 0.2, 0.7};
    markov_init(&tm, 3, P);
    double pi[ET_MAX_STATES];
    int rc = markov_stationary_distribution(&tm, pi);
    TEST("full_pipeline_stationary");
    if (rc != 0) { FAIL("no convergence"); return; }
    PASS();
    
    BudgetState bs = {.resources = {5.0, 10.0, 20.0}, .n = 3};
    double pred = predict_consumption(pi, &bs);
    TEST("full_pipeline_prediction");
    CHECK(pred > 0 && pred < 20.0, "prediction in valid range");
}

static void test_wasserstein_triangle(void) {
    double a[] = {1.0, 0.0, 0.0};
    double b[] = {0.0, 1.0, 0.0};
    double c[] = {0.0, 0.0, 1.0};
    double dab = wasserstein_1_distance(a, b, 3);
    double dac = wasserstein_1_distance(a, c, 3);
    double dbc = wasserstein_1_distance(b, c, 3);
    TEST("wasserstein_triangle_inequality");
    CHECK(dac <= dab + dbc + TOLERANCE, "triangle inequality should hold");
}

static void test_simulate_different_trajectories(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    int traj1[100], traj2[100];
    markov_simulate(&tm, 0, 100, traj1);
    markov_simulate(&tm, 0, 100, traj2);
    int differ = 0;
    for (int t = 1; t < 100; t++) /* skip t=0, always = s0 */
        if (traj1[t] != traj2[t]) differ++;
    TEST("simulate_different_trajectories");
    CHECK(differ > 0, "two simulate calls should produce different trajectories");
}

static void test_markov_seed_reproducible(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    int traj1[100], traj2[100];
    markov_seed(1234);
    markov_simulate(&tm, 0, 100, traj1);
    markov_seed(1234);
    markov_simulate(&tm, 0, 100, traj2);
    int same = 1;
    for (int t = 0; t < 100; t++)
        if (traj1[t] != traj2[t]) same = 0;
    TEST("markov_seed_reproducible");
    CHECK(same, "same seed should produce identical trajectories");
}

static void test_is_ergodic_irreducible(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    char reason[256];
    bool ok = et_is_ergodic(&tm, reason);
    TEST("is_ergodic_irreducible_aperiodic");
    CHECK(ok, reason);
}

static void test_is_ergodic_3state(void) {
    TransitionMatrix tm;
    double P[] = {0.5, 0.3, 0.2,
                  0.1, 0.6, 0.3,
                  0.2, 0.2, 0.6};
    markov_init(&tm, 3, P);
    char reason[256];
    bool ok = et_is_ergodic(&tm, reason);
    TEST("is_ergodic_3state");
    CHECK(ok, reason);
}

static void test_not_ergodic_reducible(void) {
    /* State 1 cannot reach state 0 */
    TransitionMatrix tm;
    double P[] = {0.5, 0.5,
                  0.0, 1.0};
    markov_init(&tm, 2, P);
    char reason[256];
    bool ok = et_is_ergodic(&tm, reason);
    TEST("not_ergodic_reducible");
    CHECK(!ok, "reducible chain should not be ergodic");
}

static void test_not_ergodic_periodic(void) {
    /* Period-2 chain: deterministic flip */
    TransitionMatrix tm;
    double P[] = {0.0, 1.0,
                  1.0, 0.0};
    markov_init(&tm, 2, P);
    char reason[256];
    bool ok = et_is_ergodic(&tm, reason);
    TEST("not_ergodic_periodic");
    CHECK(!ok, "period-2 chain should not be ergodic");
}

static void test_budget_safety_margin_dimensions(void) {
    /* With uniform costs, std=0, so safety margin should be 0 */
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    BudgetState bs = {.resources = {10.0, 10.0}, .n = 2};
    double margin = budget_safety_margin(&tm, &bs, 0.01);
    TEST("budget_safety_margin_zero_for_uniform_cost");
    CHECK(fabs(margin) < TOLERANCE, "uniform costs => zero safety margin");
}

static void test_budget_safety_margin_positive_for_varying_cost(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    BudgetState bs = {.resources = {10.0, 30.0}, .n = 2};
    double margin = budget_safety_margin(&tm, &bs, 0.01);
    TEST("budget_safety_margin_positive_for_varying_cost");
    CHECK(margin > 0, "varying costs should give positive safety margin");
}

static void test_budget_safety_margin_monotonic(void) {
    /* Higher cost variance should give larger margin */
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    BudgetState bs_low = {.resources = {10.0, 20.0}, .n = 2};
    BudgetState bs_high = {.resources = {10.0, 40.0}, .n = 2};
    double m_low = budget_safety_margin(&tm, &bs_low, 0.01);
    double m_high = budget_safety_margin(&tm, &bs_high, 0.01);
    TEST("budget_safety_margin_monotonic_in_variance");
    CHECK(m_high > m_low, "higher cost variance should give larger margin");
}

static void test_occupation_converges_to_pi(void) {
    TransitionMatrix tm;
    make_symmetric_2state(&tm);
    double pi[ET_MAX_STATES];
    markov_stationary_distribution(&tm, pi);
    int traj[50000];
    markov_simulate(&tm, 0, 50000, traj);
    double mu[ET_MAX_STATES];
    occupation_measure(traj, 50000, mu, 2);
    TEST("occupation_converges_to_pi");
    CHECK(fabs(mu[0] - pi[0]) < LOOSE_TOL && fabs(mu[1] - pi[1]) < LOOSE_TOL,
          "empirical should converge to stationary");
}

static void test_periodic_chain_ergodic(void) {
    /* A chain that's nearly periodic - still ergodic check should work with enough samples */
    TransitionMatrix tm;
    double P[] = {0.01, 0.99,
                  0.99, 0.01};
    markov_init(&tm, 2, P);
    double pi[ET_MAX_STATES];
    markov_stationary_distribution(&tm, pi);
    int traj[100000];
    markov_simulate(&tm, 0, 100000, traj);
    double f[] = {1.0, 2.0};
    bool ok = ergodic_check(traj, 100000, f, pi, 2, 0.05);
    TEST("periodic_chain_ergodic_check");
    CHECK(ok, "near-periodic chain should still be ergodic with enough samples");
}

/* ===== MAIN ===== */

int main(void) {
    printf("\n=== Markov Chain Tests ===\n");
    test_stationary_symmetric();
    test_stationary_asymmetric();
    test_stationary_3state();
    test_stationary_single_state();
    test_stationary_absorbing();
    test_mixing_time_symmetric();
    test_mixing_time_single_state();
    test_simulate_length();
    test_simulate_start_state();

    printf("\n=== Ergodic Theory Tests ===\n");
    test_time_average_simple();
    test_ensemble_average_symmetric();
    test_ensemble_average_single();
    test_ergodic_check_convergence();
    test_ergodic_check_constant();
    test_birkhoff_bound_positive();
    test_birkhoff_bound_tighter_needs_more();
    test_occupation_measure_uniform();
    test_occupation_measure_biased();
    test_occupation_measure_single();

    printf("\n=== Budget Planning Tests ===\n");
    test_predict_consumption();
    test_predict_consumption_skewed();
    test_budget_adequacy_sufficient();
    test_budget_adequacy_insufficient();
    test_budget_adequacy_exact();
    test_budget_safety_margin();

    printf("\n=== Wasserstein & Control Tests ===\n");
    test_wasserstein_same_dist();
    test_wasserstein_different();
    test_wasserstein_3state();
    test_wasserstein_symmetric();
    test_consumption_drift_positive();
    test_consumption_drift_zero();
    test_control_correction_direction();
    test_control_correction_zero();
    test_lqr_step_basic();
    test_lqr_step_at_target();
    test_lqr_step_converges();

    printf("\n=== Integration Tests ===\n");
    test_full_pipeline();
    test_wasserstein_triangle();
    test_occupation_converges_to_pi();
    test_periodic_chain_ergodic();

    printf("\n=== RNG Fix Tests ===\n");
    test_simulate_different_trajectories();
    test_markov_seed_reproducible();

    printf("\n=== Ergodicity Check Tests ===\n");
    test_is_ergodic_irreducible();
    test_is_ergodic_3state();
    test_not_ergodic_reducible();
    test_not_ergodic_periodic();

    printf("\n=== Budget Formula Fix Tests ===\n");
    test_budget_safety_margin_dimensions();
    test_budget_safety_margin_positive_for_varying_cost();
    test_budget_safety_margin_monotonic();

    printf("\n========================================\n");
    printf("Passed: %d  Failed: %d  Total: %d\n",
           tests_passed, tests_failed, tests_passed + tests_failed);
    printf("========================================\n");

    return tests_failed > 0 ? 1 : 0;
}
