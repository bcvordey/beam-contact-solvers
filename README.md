# Computational Solvers for Beam-Obstacle Contact Under the Damped Normal Compliance Condition

This repository contains MATLAB solvers, convergence studies, parameter analyses, and spectral tools for simulating the dynamics of Euler-Bernoulli beams in contact with obstacles under the Damped Normal Compliance (DNC) condition. The code implements the finite element and finite difference schemes described in the peer-reviewed publication:

> Saylor, G.; Shillor, M.; Vordey, C. "Model and Simulations of Contact Between a Vibrating Beam and an Obstacle Using the Damped Normal Compliance Condition." *Axioms* 2025, 14(12), 866.

and in the doctoral dissertation:

> Vordey, C. "Dynamics and Vibrations of an Euler-Bernoulli Beam in Contact with Obstacles under the Damped Normal Compliance Condition." Ph.D. Dissertation, Oakland University, 2026.

## Purpose

Classical contact models in structural mechanics assume either rigid obstacles (the Signorini condition, which cannot capture energy loss) or purely elastic contact response (the standard Normal Compliance condition, which does not account for damping). The DNC condition replaces these idealizations with a physically realistic formulation that includes both a stiffness term and an explicit damping term at the contact interface. This means the model captures what occurs during real impact: energy is dissipated through heat, deformation, and internal friction rather than being recovered elastically.

The solvers in this repository allow researchers and engineers to simulate beam-obstacle contact under the DNC condition, explore how stiffness and damping parameters govern penetration depth, rebound dynamics, and settling behavior, and reproduce the computational results reported in the published paper. All code is publicly available so that other researchers can validate, extend, and build upon these methods without deriving the numerical schemes or coding them from scratch.

## Repository Structure

### `solvers/`

Production solvers for beam-obstacle contact problems. Each solver can be run with user-defined parameters to explore different contact configurations.

- `beam_1obs_implicit_fd.m` -- Implicit finite difference solver for a beam with unilateral contact at the tip under the DNC condition.
- `beam_1obs_fem.m` -- Finite element solver using Hermite cubic elements and Newmark time integration for a beam with a single obstacle.
- `beam_distributed_obstacle_fd.m` -- Implicit finite difference solver for a beam with a distributed obstacle over the spatial domain.
- `beam_2obs_fd.m` -- Finite difference solver for a beam with two bilateral obstacles (upper and lower barriers).
- `beam_2obs_fem.m` -- Finite element solver for two-obstacle bilateral contact.
- `beam_3d_bilateral_fem.m` -- Three-dimensional beam finite element solver with bilateral contact at the free end.
- `rod_contact_fd.m` -- Finite difference solver for a simplified one-dimensional rod (wave equation) with DNC contact, included for comparison with the fourth-order beam model.

### `convergence_analysis/`

Convergence studies that validate the numerical schemes. The published paper demonstrates convergence at a rate higher than 1, establishing the reliability of the solvers. These scripts reproduce and extend those results.

- `spatial_convergence_1obs.m` -- Spatial discretization convergence under mesh refinement.
- `temporal_convergence_1obs.m` -- Temporal discretization convergence under time-step refinement.
- `coupled_convergence_1obs.m` -- Joint space-time convergence analysis.
- `fem_convergence_demo.m` -- Finite element convergence order demonstration.
- `beam_contact_convergence.m` -- General beam-obstacle contact convergence study.
- `distributed_obstacle_convergence.m` -- Convergence analysis for the distributed obstacle configuration.

### `comparison_studies/`

Scripts comparing numerical methods, solver formulations, and parameter regimes.

- `fem_vs_fd_vanishing_viscosity.m` -- Convergence study comparing finite element and finite difference solutions as viscosity approaches zero.
- `compare_beta_effects_1obs_fd.m` -- Finite difference study of damping parameter (beta) effects on single-obstacle contact dynamics.
- `compare_beta_effects_1obs_fem.m` -- Finite element study of damping parameter effects.
- `compare_beta_effects_distributed.m` -- Damping parameter effects on distributed obstacle contact.
- `compare_beta_psi_effects.m` -- Combined damping and obstacle profile parameter comparisons for distributed contact.

### `parameter_studies/`

Systematic exploration of the DNC parameter space. These scripts generate the quantitative data underlying the engineering applications described in the published paper, showing how stiffness (kappa), damping (beta), and beam properties (alpha) interact to determine contact behavior.

- `alpha_kappa_stiffness_maps.m` -- Heatmaps of contact metrics as functions of beam stiffness and obstacle stiffness.
- `contact_threshold_surface.m` -- Three-dimensional surface plots identifying contact and no-contact boundaries in parameter space.
- `restitution_phase_diagrams.m` -- Phase diagrams for impact restitution and damping relationships.
- `impact_dynamics_analysis.m` -- Detailed analysis of impact dynamics and effective restitution coefficients, demonstrating why the classical coefficient of restitution is not an intrinsic material property.
- `compare_kappa_effects.m` -- Effects of obstacle stiffness on two-obstacle bilateral contact.
- `beam_2obs_restitution_beta.m` -- Restitution and damping studies for bilateral obstacle configurations.
- `beam_2obs_restitution.m` -- Restitution patterns in two-obstacle beam contact.
- `compare_beta_kappa_restitution.m` -- Multi-parameter sweep across damping, stiffness, and restitution.
- `force_influence_study.m` -- Effects of external forcing on contact dynamics.
- `alpha_kappa_grid.m` -- Grid-based parameter exploration across beam stiffness and obstacle stiffness.
- `collapse_plot.m` -- Visualization of parameter collapse phenomena in the DNC model.
- `fft_vs_alpha2.m` -- Fast Fourier Transform analysis as a function of beam stiffness.
- `run_alpha2_kappa.m` -- Batch runs across the full beam stiffness and obstacle stiffness parameter space.
- `beam_multiple_parameters.m` -- Multi-parameter DNC beam studies.
- `dnc_final.m` -- Optimized parameter configurations for the DNC method.

### `spectral_analysis/`

Frequency-domain analysis of beam vibrations under DNC contact. The published paper shows that the DNC condition modifies the vibration frequency spectrum in ways that match physical intuition. These scripts produce the supporting spectral data.

- `fft_frequency_vs_beta.m` -- FFT spectra as a function of damping parameter beta.
- `fft_frequency_vs_kappa.m` -- FFT spectra as a function of obstacle stiffness kappa.
- `fft_contact_vs_free.m` -- Comparison of frequency content with and without obstacle contact.
- `frequency_response_maps.m` -- Two-dimensional frequency maps across the parameter space.
- `fft_free_beam.m` -- FFT of free beam vibration (baseline reference without contact).
- `plot_fft_vs_beta_twoobs.m` -- FFT plots for two-obstacle systems as a function of damping.
- `plot_fft_vs_kappa_twoobs.m` -- FFT plots for two-obstacle systems as a function of stiffness.
- `beam_distributed_fft.m` -- FFT analysis for the distributed obstacle configuration.

### `reference_solutions/`

Analytical solutions and validation benchmarks. These provide the reference data against which numerical convergence is measured.

- `clamped_free_eigenvalues.m` -- Natural frequencies and mode shapes for the clamped-free Euler-Bernoulli beam.
- `clamped_free_eigenvalues2.m` -- Alternative eigenvalue computation method.
- `steady_state_1obs.m` -- Analytical steady-state solution for single-obstacle contact.
- `steady_state_2obs.m` -- Steady-state solution for two-obstacle contact.
- `steady_state_rod.m` -- Steady-state rod solution.
- `compare_numerical_exact.m` -- Validation of numerical solutions against analytical benchmarks.
- `compare_numerical_exact3.m` -- Extended validation studies.
- `solve_steady_beam_contact_vectorf.m` -- Steady-state solver with vector forcing.

### `examples/`

Introductory scripts demonstrating solver usage for common configurations.

- `free_beam_vibration.m` -- Free vibration of an Euler-Bernoulli beam without contact (baseline case).
- `beam_with_contact_basic.m` -- Basic beam with DNC contact at a single obstacle.
- `forced_beam.m` -- Beam with external forcing and no contact.
- `forced_beam_with_contact.m` -- Forced beam with DNC obstacle contact.
- `fem_tip_trace_example.m` -- Finite element tip displacement trace over time.
- `beam_contact_bar_example.m` -- Contact with a bar-type obstacle.
- `plot_tip_with_pluck.m` -- Initial pluck excitation example.
- `rod_fd_solver.m` -- Basic rod contact solver example.
- `beam_1obs_implicit_fd_alt.m` -- Alternative single-obstacle finite difference implementation.
- `beam_contact_force_example.m` -- Contact force visualization example.

### `animations/`

Scripts generating animations and GIF files of beam contact dynamics under varying DNC parameters. These animations are deployed as a public visualization page at [https://bcvordey.github.io/beam-animations/](https://bcvordey.github.io/beam-animations/).

- `beam_1obs_animation.m` -- Animation for single-obstacle beam contact.
- `beam_1obs_gif.m` -- GIF generation for single-obstacle dynamics.
- `beam_2obs_animation.m` -- Animation for two-obstacle bilateral contact.
- `beam_2obs_animation2.m` -- Alternative two-obstacle animation with different parameter settings.
- `beam_2obs_fdd.m` -- Finite-difference-based two-obstacle animation.
- `beam_2obs_fdt.m` -- Time-domain two-obstacle animation.
- `index.html` -- Web-based animation viewer for local preview.

### `utilities/`

Helper functions for plotting, obstacle definition, and auxiliary analysis.

- `make_contact_audio.m` -- Sonification of contact forces, generating WAV audio files that represent contact dynamics for different obstacle types.
- `plot_beam_configuration.m` -- Steady-state beam configuration plotting.
- `plot_beam_steady.m` -- Steady-state beam visualization.
- `plot_steady_beam_distributed.m` -- Plotting for distributed obstacle configurations.
- `obstacle_profile.m` -- Function defining obstacle geometry.
- `distributed_obs_profile.m` -- Profile definition for distributed obstacles.
- `obs.m` -- Obstacle utility functions.
- `psi_to_latex_rhs.m` -- LaTeX output generation for obstacle profile functions.

## Getting Started

**Prerequisites:** MATLAB with standard toolboxes. No additional commercial toolboxes are required.

**Running a solver:** Open any script in `solvers/` and run it directly. Physical and numerical parameters are defined at the top of each script and can be modified to explore different contact scenarios.

**Key parameters:**

- `kappa` -- Obstacle stiffness. Controls the force response at the contact interface. Values in the published paper range from 0.1 (soft contact) to 1000 (hard contact approaching the Signorini limit).
- `beta` -- Damping coefficient. Controls the rate of energy dissipation during contact. Values in the published paper range from 0.1 to 5.
- `alpha` -- Beam stiffness parameter (alpha squared equals EI divided by rho times A, where EI is the flexural rigidity and rho A is the mass per unit length).
- `delta` -- Contact method selector (delta = 1 activates the DNC formulation).
- `Nx`, `Nt` -- Number of spatial and temporal grid points.

**Outputs:** Most scripts generate plots of beam displacement, contact force, tip trajectory, and energy evolution. Animation scripts produce GIF files or MATLAB movie objects. The `make_contact_audio.m` utility generates WAV files that sonify the contact force signal.

## Contact Configurations

The solvers support five contact configurations:

- **Single obstacle:** Unilateral contact at the beam tip or at a specified location along the beam.
- **Two obstacles:** Bilateral constraints with upper and lower barriers, modeling confined vibration.
- **Distributed obstacle:** Contact over a continuous region of the beam length, modeling surface-to-surface interaction.
- **Three-dimensional beam:** Spatial beam with bilateral contact at the free end.
- **Rod (wave equation):** Simplified one-dimensional wave equation with DNC contact, included for comparison with the fourth-order Euler-Bernoulli beam model.

## Numerical Methods

- **Finite Difference (FD):** Implicit and explicit schemes for spatial and temporal discretization of the beam PDE with DNC contact boundary conditions.
- **Finite Element (FEM):** Hermite cubic elements for spatial discretization with Newmark-beta time integration, providing higher-order accuracy.
- **Fast Fourier Transform (FFT):** Spectral analysis of beam displacement and contact force signals in the frequency domain.
- **Analytical solutions:** Closed-form steady-state solutions and eigenvalue computations for validation benchmarks.

## Engineering Relevance

The parameter studies and solver capabilities in this repository support research and design applications in several domains:

- **MEMS and semiconductor reliability**, where micro-scale contact interfaces undergo millions of loading cycles and device performance depends on energy dissipation at contact surfaces.
- **Transportation safety**, where impact dynamics in crash detection systems, brake disc-pad contact, and suspension components involve finite stiffness and damping at contact interfaces.
- **Infrastructure resilience**, where bridge expansion joints, pavement-soil contact, and railway wheel-rail interfaces experience repetitive impact loading that can be modeled using the DNC framework.
- **Robotics and prosthetic devices**, where controlled contact requires optimized stiffness and damping parameters to achieve safe, efficient manipulation and locomotion.

## Companion Repositories

- **Python solvers:** [https://github.com/bcvordey/dnc-beam-python](https://github.com/bcvordey/dnc-beam-python)
- **Beam-obstacle contact animations:** [https://github.com/bcvordey/beam-animations](https://github.com/bcvordey/beam-animations)

## Citation

If you use this code in your research, please cite:

```
Saylor, G.; Shillor, M.; Vordey, C. "Model and Simulations of Contact Between
a Vibrating Beam and an Obstacle Using the Damped Normal Compliance Condition."
Axioms 2025, 14(12), 866.
```

## Author

Cornelius Bright Vordey
Ph.D. Candidate, Applied Mathematical Sciences
Oakland University, Rochester, Michigan
ORCID: [0009-0005-1644-9827](https://orcid.org/0009-0005-1644-9827)