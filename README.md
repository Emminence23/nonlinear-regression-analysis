# Nonlinear Regression Analysis

Course project for **MTH 686 (Nonlinear Regression)** at IIT Kanpur, under Prof. Debasis Kundu.

The task: fit three competing nonlinear regression families to a data set of 70 observations. For each, find the least-squares estimators and explain how they were found, choose the best model, estimate σ², give Fisher-information confidence intervals, check the residuals and normality, and plot the fits ([problem statement](docs/problem_statement.pdf)). This repository implements the whole analysis in R, including its own Levenberg–Marquardt solver and Osborne's modified Prony algorithm. **[Read the report (PDF)](report/report.pdf)**, a revised and extended version of the six-page report submitted in November 2025.

$$
\begin{aligned}
\text{Model 1:}\quad & y = \alpha_0 + \alpha_1 e^{\beta_1 t} + \alpha_2 e^{\beta_2 t} + \varepsilon \\
\text{Model 2:}\quad & y = (\alpha_0 + \alpha_1 t)/(\beta_0 + \beta_1 t) + \varepsilon \\
\text{Model 3:}\quad & y = \textstyle\sum_{k=0}^{4} \beta_k t^k + \varepsilon
\end{aligned}
$$

![Data and fitted curves](figures/fits.png)

## Results

| Model | p | RSS | σ̂² | AIC | BIC | Notes |
|---|---:|---:|---:|---:|---:|---|
| Model 1, as given | 5 | → 8.4282 | | | | **no least-squares estimate exists**: rates merge, amplitudes → ±∞ |
| Model 1, confluent limit | 4 | 8.4282 | 0.1277 | **60.47** | 71.71 | closure of Model 1; attains its infimum |
| **Model 1, one exponential** | 3 | 8.7086 | 0.1300 | 60.76 | **69.75** | **selected** |
| Decaying wave (extension) | 5 | 8.4262 | 0.1296 | 62.45 | 75.94 | best curve of Model 1's general form; wave is noise (p = 0.22) |
| Model 2 (β₀ = 1) | 3 | 9.1039 | 0.1359 | 63.87 | 72.86 | identifiable but weakly so |
| Model 3 | 5 | 10.9099 | 0.1678 | 80.53 | 94.02 | |

**Selected model:** ŷ(t) = 0.761 + 3.730 e^(−0.372 t), σ̂² = 0.130. The second exponential term is not supported (F = 2.20, p = 0.14). The residuals pass Shapiro–Wilk (p = 0.68), Anderson–Darling and Jarque–Bera, and a runs test shows no systematic lack of fit.

## What the analysis found

**1. Model 1 has no interior least-squares estimate.** The profiled RSS keeps decreasing as the two rates merge (β₂ − β₁ → 0) while α₁ = −α₂ → ±∞. Levenberg–Marquardt follows that path indefinitely, and R's `nls(algorithm = "plinear")` fails. The infimum is attained only by the confluent limit α₀ + (α₁ + γt)e^{βt}. A naive fit therefore reports huge standard errors (about 7 × 10⁵ for the amplitudes where Levenberg–Marquardt stops, and growing the longer it runs), because J'J is numerically singular (condition number about 6 × 10¹²).

<p align="center"><img src="figures/model1_rss_surface.png" width="520" alt="Profiled RSS surface of Model 1"></p>

**2. An independent check confirms it.** Osborne's modified Prony algorithm reparametrises the problem by the coefficients of the characteristic polynomial of the exponential recurrence. In those coordinates the global optimum has **complex-conjugate roots** 0.473 ± 0.140i, a damped oscillation that a real-rate model can't represent. The constrained optimum therefore sits on the double-root boundary, which matches the variable-projection result.

**3. Model 2 as written is not identifiable.** Scaling (α₀, α₁, β₀, β₁) by any constant leaves the curve unchanged, so the Fisher information is singular. Fixing β₀ = 1 restores identifiability. Even so, the profile-likelihood intervals are unbounded above (the data mainly determine α₁/β₁), which symmetric Wald intervals hide.

**4. The "wave" is noise.** Fitting the decaying wave a₀ + e^{ρt}(c cos ωt + d sin ωt) directly gives exactly Osborne's RSS, a third independent check. It improves on one exponential by only 0.28. On 1000 fake data sets made from the one-exponential model (so with no wave at all), the wave did at least that well 22% of the time. The improvement is what noise alone produces.

**5. The method checks out on simulated data.** On 1000 fake data sets per model, where the right answer is known:

| | right model chosen | 95% intervals actually cover | two-exponential fit breaks down |
|---|---:|---:|---:|
| true model: one exponential | 76% | 95–96% | 61% |
| true model: rational | 75% | 85–87% | 33% |
| true model: quartic | 50% | 93–94% | 82% |

The quartic is a linear model, so its intervals must cover exactly 95%. With 1000 data sets the estimate itself varies by about ±1.4 points; a separate run of 10,000 data sets gives 94.7–95.3%, which confirms the simulation itself. The selected model's intervals can be trusted, while the rational model's fall short. Model 1's breakdown on data set 94 was typical, not bad luck.

## Method

- **Variable projection (partial linearity).** Each model is linear in its amplitudes once the nonlinear parameters are fixed. Profiling out the linear part in closed form, S*(β) = yᵀ(I − P_X(β))y, turns Model 1's five-dimensional non-convex problem into a two-dimensional one that is scanned exhaustively on a grid.
- **Nonlinear regression (NLR) refinement with Levenberg–Marquardt** of all parameters jointly from the best grid point, with analytic Jacobians and convergence and gradient-norm checks ([`R/fit.R`](R/fit.R)).
- **Osborne's modified Prony algorithm** for sums of exponentials, extended to Model 1's constant term by fixing a root of the characteristic polynomial at 1, and run multi-start ([`R/osborne.R`](R/osborne.R)). The derivation corrects a sign in the course handout (not included here), and `run_analysis.R` checks the implementation against finite differences.
- **Inference.** σ̂² = RSS/(n − p). Fisher-information (Wald) intervals are cross-checked against profile-likelihood intervals and 2000-replicate residual-bootstrap intervals. Identifiability is measured by the condition number of J'J. Model selection uses AIC, AICc and BIC plus an extra-sum-of-squares F test for nested models ([`R/inference.R`](R/inference.R)).
- **Simulation checks.** Parametric simulation from each fitted model, refitting with the full pipeline: model-selection accuracy, interval coverage, and a simulation-based test for the wave ([`R/simulation.R`](R/simulation.R)).
- **Diagnostics.** Residuals against t and fitted values, Q–Q plots, Shapiro–Wilk, Anderson–Darling, Jarque–Bera and Wald–Wolfowitz runs tests ([`R/diagnostics.R`](R/diagnostics.R)).

## Reproduce

Requires R ≥ 4.3 with `ggplot2`, `patchwork`, `nortest` and `knitr`, plus [Quarto](https://quarto.org) with TinyTeX for the PDF.

```bash
Rscript -e 'install.packages(c("ggplot2", "patchwork", "nortest", "knitr"))'
Rscript simulate.R                  # simulation checks, all cores (≈6 min on 8 cores)
Rscript run_analysis.R              # results/*.csv, figures/*.png, report tables (≈1 min)
quarto render report/report.qmd     # report/report.pdf
```

Every estimate, test result and table in the report is read from the files `run_analysis.R` writes; only method settings (grid sizes, number of starts) are written into the text. The simulation gives the same numbers on any number of cores. Package versions used are recorded in [`results/session_info.txt`](results/session_info.txt).

## Layout

```
data/set-94.dat          the data (t = 1..70)
R/models.R               mean functions and analytic Jacobians
R/fit.R                  variable-projection grid searches, Levenberg–Marquardt
R/osborne.R              Osborne's modified Prony algorithm
R/inference.R            σ², information criteria, Wald / profile / bootstrap CIs, F test
R/diagnostics.R          residual tests and figures
R/simulation.R           simulation checks
run_analysis.R           entry point
simulate.R               simulation entry point (results cached in results/simulation/)
report/report.qmd        report source  →  report/report.pdf
results/                 CSV outputs
figures/                 PNG figures
docs/                    problem statement
```
