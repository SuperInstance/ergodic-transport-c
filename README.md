# ergodic-transport-c

Ergodic theory + optimal transport for agent budget planning. Uses Birkhoff's ergodic theorem, Perron-Frobenius stationary distributions, and Wasserstein distance to predict and control agent resource consumption over time.

## Architecture

### Markov Chain (`src/markov.c`)
- Transition matrix with power iteration for stationary distribution (Perron-Frobenius)
- Mixing time estimation via total variation distance
- Stochastic simulation with deterministic RNG for reproducibility

### Ergodic Theory (`src/ergodic.c`)
- Time averages via Birkhoff's ergodic theorem: `(1/N) Σ f(X_t)`
- Ensemble averages: `Σ π_i f(i)` 
- Convergence verification and sample complexity bounds
- Empirical occupation measure computation

### Budget Planning (`src/budget.c`)
- Long-run resource consumption prediction using ergodic theorem
- Budget adequacy checking with margin calculation
- Safety margins for non-stationary mixing periods

### Wasserstein Control (`src/wasserstein_control.c`)
- W₁ distance via CDF formula between occupation measures
- Consumption drift detection
- Control corrections proportional to state deviations
- Linear-quadratic regulator (LQR) for budget adjustments

## Building

```bash
make          # build library and tests
make test     # run test suite
make clean    # clean artifacts
```

## Requirements

- C11 compiler (gcc/clang)
- libm (math library)

## Test Suite

35 tests covering Markov chains, ergodic convergence, budget planning, Wasserstein distance, and LQR control.

## License

MIT
