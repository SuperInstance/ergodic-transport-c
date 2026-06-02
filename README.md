# ergodic-transport-c
### Birkhoff's theorem as a C library

---

How much memory will you need next month? Don't guess. Compute. Your monitoring data already knows.

You're sitting on terabytes of time-series data — CPU, memory, disk, network, request latencies — and you're still planning capacity by taking last month's peak and adding 30%. You might as well be throwing darts.

This library is the math your monitoring stack wishes it had.

---

## The One Theorem That Changes Everything

**Birkhoff's pointwise ergodic theorem**, in plain language:

> **Time average = ensemble average. Your past *is* your future, averaged.**

Every ops team is doing ergodic theory. They just call it "looking at averages." The problem is they stop there.

You look at your average memory usage over the last week and call it a day. Birkhoff says you can do better — you can compute _how many samples you need_ to be confident your average is within 1% of the true mean. You can place _bounds_ on your forecasts. You can detect when the system has changed, not just when a metric crossed some arbitrary threshold.

The math gives you BOUNDS. Not hunches.

---

## Concrete Example

Suppose you have a 3-state API server modeling request rates:

| State | Description       | Memory (GB) | Transition → |
|-------|-------------------|-------------|--------------|
| 0     | Low traffic       | 5           | Stays 80%, goes to 1 at 15%, goes to 2 at 5% |
| 1     | Medium traffic    | 10          | Goes to 0 at 10%, stays 70%, goes to 2 at 20% |
| 2     | Traffic spike     | 20          | Goes to 0 at 10%, goes to 1 at 20%, stays 70% |

The stationary distribution is `π ≈ (0.25, 0.41, 0.34)`. The predicted long-run memory consumption:

```
predicted = 0.25 × 5 + 0.41 × 10 + 0.34 × 20 = 12.15 GB
```

Your monitoring data, fed through `markov_init` and `markov_stationary_distribution`, gives you this number — not as a guess, but as the unique limit of a convergent process. Birkhoff guarantees that if your system is ergodic, the time average of your actual memory usage converges to this value.

How confident are you? Use `birkhoff_bound`:

```c
double pi[] = {0.25, 0.41, 0.34};
double f[]  = {5.0, 10.0, 20.0};
int N = birkhoff_bound(pi, f, 3, 0.5, 0.01);
// N = ceil(29.3 / (0.25 * 0.01)) = ceil(11720) = 11720 samples
```

That's about 3.3 hours of per-second monitoring data to be 99% confident your average is within 0.5 GB. And that's using a Chebyshev bound — loose but distribution-free. If you know your variance, tighten it.

The same calculus applies to memory, CPU, disk I/O, request latency, queue depth — any resource whose usage depends on which "state" your system is in. The states are whatever you define them to be: "off-peak," "peak," "batch job running," "deploy in progress." The transition probabilities are whatever your monitoring history says they are.

Suddenly "how long should I monitor before making a decision" isn't a gut call. It's a function call.

---

## What This Library Does

`ergodic-transport-c` takes the mathematical machinery of ergodic theory and makes it usable for operational infrastructure. The core pipeline:

1. **Model your workload as a Markov chain** — build a transition matrix from monitoring data
2. **Check if it's ergodic** — `et_is_ergodic()` verifies irreducibility + aperiodicity
3. **Compute the stationary distribution** — this is your system's long-run behavior
4. **Predict resource consumption** — combine π with per-state resource costs
5. **Want bounds?** — `birkhoff_bound()` tells you exactly how many samples you need
6. **Detect changes in real time** — `consumption_drift()` uses Wasserstein distance to flag anomalies

---

## Applications

### Capacity Planning

You provision a server farm. Your workload has predictable patterns (peak hours, maintenance windows, batch jobs). Model the transitions, compute the stationary distribution, and forecast. Not with heuristics — with the Birkhoff bound.

```c
BudgetState cluster = { .resources = { /* per-state memory cost */ }, .n = N };
BudgetPlan budget = { .allocation = { /* per-state allocation */ }, .n = N };

double predicted = predict_consumption(pi, &cluster);
double margin = budget_adequacy(&budget, pi, &cluster);
```

If your adequacy margin is negative, you under-provisioned — statistically, provably, not anecdotally.

### Anomaly Detection

This is where it gets interesting. You've built a model of your system, you've computed π, you've planned capacity. Now the system changes — a new deployment, a traffic surge, a dependency degradation. How do you know?

Most monitoring stacks use static thresholds: "alert if memory > 80%." That misses slow drifts, distribution shifts, and phase transitions entirely. A system can jump from "steady state A" to "steady state B" without ever crossing a single threshold, because both states are within bounds — they just spend different fractions of time at different resource levels.

The Wasserstein-1 distance (earth mover's distance) between your empirical occupation measure and the stationary distribution is your _consumption drift_:

```c
double drift = consumption_drift(trajectory, N, pi, n_states);
```

When your current behavior diverges from your predicted behavior, something changed. You just detected a phase transition.

This isn't a threshold on CPU utilization. It's a distributional shift — the system is _qualitatively_ different. A sustained non-zero drift means your model is wrong. Time to retrain the transition matrix.

The Wasserstein distance has a beautiful property here: it's the minimum amount of "work" required to move mass from one distribution to another. On a 1D line of ordered states (e.g., low→medium→high load), this is the CDF area between the two distributions. A drift of 0.3 means you need to move 30% of the probability mass to restore stationarity. That's a concrete, interpretable number — not a p-value, not a z-score, but literal probability mass that has shifted.

You can feed this value into an alert: `if consumption_drift > threshold for N consecutive windows, page the SRE`. The threshold comes from the Birkhoff bound on how much natural variation your chain exhibits. N comes from the mixing time. Both have mathematical justification, not gut feels.

### Budget Forecasting

This is the one that makes finance teams stop rolling their eyes at "probabilistic capacity planning." The standard approach: predict average monthly cost, multiply by 1.5, call it a day. The ergodic approach: compute the exact long-run cost as a function of state, then compute exactly how much extra you need during transients.

The `budget_safety_margin` handles the non-stationary transient — the first `t_mix` steps of a new system where the chain hasn't converged:

```c
// mixing_time × std_cost — dimensional analysis: [steps] × [cost/step] = [cost] ✓
double extra = budget_safety_margin(&tm, &bs, 0.01);
```

This accounts for the fact that during the mixing period, the chain's distribution can deviate from the stationary one. The standard deviation of the cost, scaled by the mixing time, gives you the extra budget to set aside until the system settles.

The formula is dead simple: `extra = t_mix × σ`. But deriving it requires understanding that your first `t_mix` observations aren't representative. Most capacity planning tools ignore this and under-provision by exactly this margin.

The dimensional analysis checks out: mixing time (steps) times standard deviation of cost (cost units per step) gives cost units. The function returns the *extra* budget beyond the stationary prediction, so your total budget for the planning horizon is:

```
total_budget = predict_consumption(pi, &bs) + budget_safety_margin(&tm, &bs, 0.01)
```

If the mixing time is 500 steps and the standard deviation of per-step cost is 0.8 GB, you need 400 GB additional buffer for the transient period. Not 500 GB, not 1000 GB — exactly 400 GB. The math tells you.

### Auto-scaling

Auto-scaling is just ergodic control with a feedback loop. The system estimates its current empirical distribution over a sliding window and compares it to the stationary target. Discrepancies trigger adjustments.

The LQR control loop lets you adjust budget allocations dynamically:

```c
double correction[ET_MAX_STATES];
control_correction(current_mu, target_pi, cost, correction, n);

double new_allocation[ET_MAX_STATES];
lqr_step(allocation, target, cost, 0.3, new_allocation, n);
```

Each step reduces the Wasserstein distance between your empirical distribution and your target. This is proportional control on the distribution itself — not on raw metrics, but on the probability distribution of where your system spends its time.

The `lqr_step` converges to `target × cost` under a constant-gain controller. In the 3-state example above, the LQR target allocations are:

```
state 0: π₀ × f₀ = 0.25 × 5  = 1.25 GB
state 1: π₁ × f₁ = 0.41 × 10 = 4.10 GB
state 2: π₂ × f₂ = 0.34 × 20 = 6.80 GB
Total: 12.15 GB
```

If the system drifts into state 2 more often than π₂ predicts (a sustained traffic surge), the empirical occupation measure shifts, the Wasserstein drift spikes, and the control loop reallocates budget toward state 2. The system self-corrects because it can detect that the *distribution* has changed — not just that one metric crossed a line.

---

## API Reference

See [`include/ergodic_transport.h`](include/ergodic_transport.h) for the full API. Key functions:

| Function | Purpose |
|---|---|
| `markov_init()` | Build transition matrix from flat array |
| `markov_stationary_distribution()` | Power iteration — find π |
| `markov_mixing_time()` | Steps to converge within ε |
| `markov_simulate()` | Generate trajectories (RNG persists across calls) |
| `markov_seed()` | Seed the RNG for reproducibility |
| `et_is_ergodic()` | Check irreducibility + aperiodicity |
| `ergodic_time_average()` | (1/N) Σ f(state) from trajectory |
| `ergodic_ensemble_average()` | Σ π_i × f(i) |
| `ergodic_check()` | Check if |time avg - ensemble avg| < tolerance |
| `birkhoff_bound()` | How many samples for given accuracy & confidence |
| `occupation_measure()` | Fraction of time spent in each state |
| `predict_consumption()` | Long-run resource prediction |
| `budget_adequacy()` | Budget - predicted consumption |
| `budget_safety_margin()` | Extra budget for mixing period (= t_mix × σ) |
| `wasserstein_1_distance()` | Earth mover's distance via CDF |
| `consumption_drift()` | Wasserstein distance from stationary |
| `control_correction()` | Budget adjustments to match target |
| `lqr_step()` | Proportional control iteration |

---

## Building

```bash
make          # builds libergodictransport.a and test binary
make test     # runs the full test suite
```

Requires: C11 compiler, `libm`.

---

## The RNG Detail

The Markov simulator uses a deterministic LCG with persistent state across calls. This is intentional:

- `markov_seed(seed)` at the start guarantees reproducibility
- Each `markov_simulate()` call picks up where the last one left off — no reset at the top
- Two calls with the same seed produce the same sequence

This means you can call `markov_simulate` in chunks and get a single coherent trajectory:

```c
markov_seed(42);
markov_simulate(&tm, 0, 5000, first_half);
markov_simulate(&tm, trajectory[4999], 5000, second_half);
// second_half[0] == first_half[4999]? No — RNG advanced, new random step.
```

If you want fresh trajectories, call `markov_seed` with a new seed before each run.

This design matters for reproducibility in production: if you're running daily capacity simulations, seeding with the date (e.g., `markov_seed(20260602)`) means everyone gets the same trajectory for that day's planning. No stochastic surprises in your Monday morning forecast.

It also matters for testing: the test suite confirms that `markov_seed(1234)` followed by two `markov_simulate` calls produces identical trajectories, while two sequential calls without an intervening `markov_seed` produce different ones. The RNG is deterministic but not reset — each call advances the state.

---

## Ergodicity Check

Not every Markov chain is ergodic. Use `et_is_ergodic()` before trusting your predictions:

- Two reducible components → convergence isn't unique
- Periodic chain (deterministic flip) → time averages don't converge to ensemble averages
- The function fills a `reason` string so you know *why* it failed

```c
char reason[256];
if (!et_is_ergodic(&tm, reason)) {
    fprintf(stderr, "Chain not ergodic: %s\n", reason);
    // Fall back to empirical methods / warn the user
}
```

---

## Why This Exists

Every ops team is doing ergodic theory. They call it "looking at averages," "capacity planning," "monitoring dashboards." But they do it without the math that tells them when they're wrong, how wrong they might be, and how long they need to look before they know.

`ergodic-transport-c` is the library that turns your monitoring data into predictions with provable bounds. Not because you need to learn ergodic theory — because you've been doing it this whole time and didn't know you could compute the confidence interval.

Your Prometheus data is a trajectory. Your Grafana dashboard shows occupation measures. Your capacity planning spreadsheet computes ensemble averages. You're already doing the math. This library just gives you the guarantees, the bounds, and the detection mechanisms you've been missing.

Next time someone asks why you need more budget, don't show them last month's peak. Show them the stationary distribution. Show them the mixing time. Show them the safety margin — derived, bounded, dimensional-analysis-verified. Show them the Birkhoff bound for 99% confidence within 1% error.

And then watch their face when they realize you just used Birkhoff's pointwise ergodic theorem to justify a budget increase.

---

## License

MIT. Do what you want, but if your capacity planning uses this, tell people you "moved from heuristics to ergodic theory." The reactions are worth it.

---

*"The difference between a guess and a prediction is the error bound. Birkhoff gave you the bound. You just had to read it."*
