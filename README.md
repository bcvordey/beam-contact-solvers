# Dissertation MATLAB Solvers for Beam Contact Problems

This repository contains MATLAB solvers and analysis scripts for simulating beam contact problems, specifically focusing on Euler-Bernoulli beams with unilateral or bilateral contact constraints. These solvers implement various numerical methods (Finite Difference, Finite Element, FFT-based) to solve dynamic contact problems with obstacles.

The code is organized into folders for easy navigation and public use. Each solver can be run with modified parameters to explore different contact scenarios.

## Repository Structure

### `solvers/`
Core production solvers for beam contact problems. These are the main reusable functions/scripts.

- `beam_1obs_implicit_fd.m`: Implicit finite difference solver for 1D beam with unilateral contact at the tip (δ=1, DNC method).
- `beam_1obs_fem.m`: Finite element solver (Hermite elements, Newmark time integration) for beam with single obstacle.
- `beam_distributed_obstacle_fd.m`: Implicit FD solver for beam with distributed obstacle over the domain.
- `beam_2obs_fd.m`: FD solver for beam with two bilateral obstacles.
- `beam_2obs_fem.m`: FEM solver for beam with two obstacles.
- `beam_3d_bilateral_fem.m`: 3D beam FEM solver with bilateral contact at free end.
- `rod_contact_fd.m`: Simplified 1D rod (wave equation) with contact, using FD methods.

### `comparison_studies/`
Scripts comparing different numerical methods, parameters, or approaches.

- `fem_vs_fd_vanishing_viscosity.m`: Convergence study comparing FEM and FD as viscosity → 0.
- `compare_beta_effects_1obs_fd.m`: FD study of damping (β) effects on 1-obstacle beam.
- `compare_beta_effects_1obs_fem.m`: FEM study of β effects.
- `compare_beta_effects_distributed.m`: β effects on distributed obstacle.
- `compare_beta_psi_effects.m`: β and ψ parameter comparisons for distributed contact.

### `convergence_analysis/`
Error analysis and convergence studies for numerical methods.

- `spatial_convergence_1obs.m`: Spatial discretization convergence (Δx refinement).
- `temporal_convergence_1obs.m`: Temporal discretization convergence (Δt refinement).
- `coupled_convergence_1obs.m`: Joint space-time convergence analysis.
- `fem_convergence_demo.m`: FEM convergence order demonstration.
- `beam_contact_convergence.m`: General beam contact convergence study.
- `distributed_obstacle_convergence.m`: Convergence for distributed obstacles.

### `parameter_studies/`
Exploration of parameter spaces (damping, stiffness, restitution, etc.).

- `alpha_kappa_stiffness_maps.m`: Heatmaps of contact metrics vs beam stiffness (α²) and obstacle stiffness (κ).
- `contact_threshold_surface.m`: 3D surface plotting contact/no-contact boundaries.
- `restitution_phase_diagrams.m`: Phase diagrams for impact restitution and damping.
- `impact_dynamics_analysis.m`: Detailed analysis of impact dynamics and restitution coefficients.
- `compare_kappa_effects.m`: Effects of obstacle stiffness (κ) on 2-obstacle system.
- `beam_2obs_restitution_beta.m`: Restitution and β studies for 2 obstacles.
- `beam_2obs_restitution.m`: Restitution patterns in 2-obstacle beam.
- `compare_beta_kappa_restitution.m`: Multi-parameter sweep (β, κ, restitution).
- `force_influence_study.m`: How external forcing affects contact dynamics.
- `alpha_kappa_grid.m`: Grid-based parameter exploration for α and κ.
- `collapse_plot.m`: Visualization of parameter collapse phenomena.
- `fft_vs_alpha2.m`: FFT analysis vs beam stiffness.
- `run_alpha2_kappa.m`: Batch runs for α² × κ parameter space.
- `beam_multiple_parameters.m`: Multi-parameter DNC beam studies.
- `dnc_final.m`: Final parameter optimization for DNC method.

### `spectral_analysis/`
Frequency domain analysis using FFT and spectral methods.

- `fft_frequency_vs_beta.m`: FFT spectra vs damping parameter β.
- `fft_frequency_vs_kappa.m`: FFT spectra vs stiffness κ.
- `fft_contact_vs_free.m`: Comparison of frequency content with/without contact.
- `frequency_response_maps.m`: 2D frequency maps for parametric studies.
- `fft_free_beam.m`: FFT of free beam vibration.
- `plot_fft_vs_beta_twoobs.m`: FFT plots for 2-obstacle system vs β.
- `plot_fft_vs_kappa_twoobs.m`: FFT plots for 2-obstacle system vs κ.
- `beam_distributed_fft.m`: FFT analysis for distributed obstacle.

### `reference_solutions/`
Analytical solutions and validation against numerical methods.

- `clamped_free_eigenvalues.m`: Natural frequencies and mode shapes for clamped-free beam.
- `clamped_free_eigenvalues2.m`: Alternative eigenvalue computation.
- `steady_state_1obs.m`: Analytical steady-state solution for 1-obstacle beam.
- `steady_state_2obs.m`: Steady-state for 2-obstacle beam.
- `steady_state_rod.m`: Steady-state rod solution.
- `compare_numerical_exact.m`: Validation of numerical vs analytical solutions.
- `compare_numerical_exact3.m`: Extended validation studies.
- `solve_steady_beam_contact_vectorf.m`: Steady-state solver with vector forcing.

### `utilities/`
Helper functions, plotting utilities, and audio generation.

- `make_contact_audio.m`: Sonify contact forces (generate WAV files for different obstacle types).
- `plot_beam_configuration.m`: Plot steady-state beam configurations.
- `plot_beam_steady.m`: Steady-state beam plotting.
- `plot_steady_beam_distributed.m`: Plotting for distributed obstacles.
- `obstacle_profile.m`: Function defining obstacle shapes.
- `distributed_obs_profile.m`: Profile for distributed obstacles.
- `obs.m`: Obstacle utility functions.
- `psi_to_latex_rhs.m`: LaTeX output for ψ functions.

### `examples/`
Basic example scripts demonstrating solver usage.

- `free_beam_vibration.m`: Simple free vibration (no contact).
- `beam_with_contact_basic.m`: Basic beam with DNC contact.
- `forced_beam.m`: Beam with external forcing (no contact).
- `forced_beam_with_contact.m`: Forced beam with contact.
- `fem_tip_trace_example.m`: FEM tip displacement trace.
- `beam_contact_bar_example.m`: Contact with bar obstacle.
- `plot_tip_with_pluck.m`: Initial pluck excitation example.
- `rod_fd_solver.m`: Basic rod solver example.
- `beam_1obs_implicit_fd_alt.m`: Alternative 1-obstacle FD implementation.
- `beam_contact_force_example.m`: Force-influenced contact example.

### `animations/`
Scripts generating animations and GIFs of beam dynamics.

- `beam_1obs_animation.m`: Animation for 1-obstacle beam.
- `beam_1obs_gif.m`: GIF generation for 1-obstacle.
- `beam_2obs_animation.m`: Animation for 2-obstacle beam.
- `beam_2obs_animation2.m`: Alternative 2-obstacle animation.
- `beam_2obs_fdd.m`: FD-based 2-obstacle animation.
- `beam_2obs_fdt.m`: Time-domain 2-obstacle animation.
- `index.html`: Web-based animation viewer.
- `README.md`: Animation documentation.

## Getting Started

1. **Prerequisites**: MATLAB with basic toolboxes.
2. **Running a Solver**: Open any `.m` file in `solvers/` and run it. Modify parameters at the top of the script.
3. **Parameters**: Common parameters include:
   - `beta`: Damping coefficient
   - `kappa`: Obstacle stiffness
   - `alpha`: Beam stiffness parameter (α² = EI/ρA)
   - `delta`: Contact method (δ=1 for DNC)
   - `Nx`, `Nt`: Spatial/temporal grid points
4. **Outputs**: Most scripts generate plots. Animations require MATLAB's plotting functions.
5. **Audio**: Run `utilities/make_contact_audio.m` to generate sound files.

## Problem Types

- **1 Obstacle**: Unilateral contact at beam tip.
- **2 Obstacles**: Bilateral constraints (upper/lower barriers).
- **Distributed Obstacle**: Contact over beam length.
- **3D Beam**: Spatial beam with 3D contact.
- **Rod**: Simplified 1D wave equation with contact.

## Numerical Methods

- **FD**: Finite Difference (explicit/implicit)
- **FEM**: Finite Element (Hermite elements, Newmark)
- **FFT**: Spectral analysis for frequency domain
- **Analytical**: Closed-form steady-state solutions

## Citation

If using this code, please cite the dissertation or contact the author for proper attribution.
