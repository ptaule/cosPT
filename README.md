# CosPT

CosPT computes loop corrections to the power spectrum or bispectrum in
cosmological perturbation theory using the algorithms outlined in
[1304.1546](http://www.arxiv.org/abs/1304.1546),
[1309.3308](http://www.arxiv.org/abs/1309.3308),
[1507.06665](http://www.arxiv.org/abs/1507.06665),
[1907.10729](http://www.arxiv.org/abs/1907.10729) and
[2008.00013](http://www.arxiv.org/abs/2008.00013).

In addition to the standard Einstein-de-Sitter Standard Perturbation Theory
(EdS-SPT) kernels, the program can compute numerically kernels with a general
equation of motion, see [2008.00013](http://www.arxiv.org/abs/2008.00013) for
further explanation and example usage for a cosmology with massive neutrinos.

Anybody is free to use the software, but please cite at least
[2008.00013](http://www.arxiv.org/abs/2008.00013) in publications.

## Dependencies

- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) used to solve ODEs numerically
- [The CUBA library](http://www.feynarts.de/cuba/) library for Monte Carlo integration
- [libconfig](https://hyperrealm.github.io/libconfig/) for reading configuration files

## Building and running the program

The binary can be compiled by running `make` (with the `-j` option to
build in parallel). It is placed in the build directory and can be run as
```
cosPT configuration_file [--option argument]
```
where a configuration file is given as argument (see below). The following
command line options can be given (see also `cosPT --help`):

- `-a`, `--k_a_idx`: Index of the external wavenumber to compute in `k_a_grid`.
- `-b`, `--k_b_idx`, `-c`, `--k_c_idx`: Similarly for the second and third external wavenumbers (bispectrum).
- `-n`, `--n_evals`: Maximum number of evaluations for the Monte Carlo integration.
- `-p`, `--n_cores`: Number of cores/threads to spawn for CUBA, in addition to the master.
- `-P`, `--stdout_mode`: Output results to stdout instead of to a file.
- `-v`, `--verbosity`: Verbosity level for standard output.

## Configuration file

The program settings are set in the configuration file (by convention a
.cfg-file) that is given as argument to the program (some settings can also be
set by command line options, see `cosPT --help`).

The settings should be written in the form
```
integer = value
float   = value
boolean = true/false
string  = "value"
array   = [a,b,c]
list    = (a,b,c)
group   = {
    integer        = value,
    string         = "value",
    list_of_arrays = ([1,2], [3,4])
    # etc.
}
```
with a newline or semicolon separating settings. Lines starting with \# are comments.
Unrecognized settings are reported as an error.

##### Options that can be set in the configuration file:

- `loops` [integer, **required**]: Number of loops to compute.
- `spectrum` [string, **required**]: spectrum to compute: "powerspectrum" or "bispectrum".
- External wavenumber for which the loop correction should be computed (for bispectrum one of the wavenumbers). Either provide `k_a` or both `k_a_grid` and `k_a_idx` (the latter can also be set with a command line option):
    - `k_a` [float]: The value of the external wavenumber given in the same units as the wavenumber grid for the input linear power spectrum.
    - `k_a_grid` [string]: File containing a wavenumber grid where each point is separated by a newline.
    - `k_a_idx` [integer]: The (zero-based) index of the wavenumber in the `k_a_grid` for which the loop corrections should be computed. The command line option `--k_a_idx` takes precedence.
- For the bispectrum, analogous settings exist for `k_b` (`k_b_grid`, `k_b_idx`) and `k_c` (`k_c_grid`, `k_c_idx`). If `k_b` is not given, `k_b = k_a` is used. In stead of `k_c`, the cosine of the angle between the wavevectors `k_a` and `k_b`, `cos_ab` [float in [-1,1]], can be given (but not both). If neither is given, `k_c = k_a` is used.
- `q_min` and `q_max` [float, **required**]: The lower and upper integration limits, respectively, with same units as external wavenumbers.
- `correlations` [list, default=([0,0]) or ([0,0,0])]: List of correlations (zero-indexed) to compute, each entry given as an array of two or three integers for power spectrum and bispectrum, respectively. Ignored when `rsd = true`. Example:
```
correlations = ([0,0], [1,1]) # power spectrum: compute <delta delta> and <theta theta>
correlations = ([0,0,0])      # bispectrum: compute <delta delta delta>
```
- `dynamics` [string, **required**]: The kernel dynamics, either "eds-spt", "evolve-asymp-ics" or "evolve-eds-ics". The first corresponds to EdS-SPT kernels, the second and third are used for dynamically evolved kernels with asymptotic or EdS initial conditions, respectively.
- `input_ps_file` [string, **required**]: The (full or relative path) filename of the input linear power spectrum. The file should contain two columns: the first with the wavenumber grid, and the second with the corresponding values for the power spectrum. Lines starting with \# are ignored.
- `input_ps_rescale` [float or string, default=1]: Rescale the input linear power spectrum by a number before interpolation (does not alter the input file). One can use the string "2pi^3" or "2pi^-3" to conventionally rescale it by $`(2\pi)^3`$ or $`(2\pi)^{-3}`$, respectively.
- `stdout_mode` [bool, default=false]: Print results only to stdout, not to file (can also be set with the `--stdout_mode` command line option).
- Unless `stdout_mode = true`, the output file must be set via either of the following settings (missing directories are created):
    - `output_file` [string]: The (full or relative path) filename to write to.
    - `output_path` [string]: The (full or relative) directory in which the output file should be put. The output file is named after `k_a` or `k_a_idx`, depending on which was set. For bispectrum, `k_b` and `k_c` (or the corresponding indices) are included in the name of the file.
- `description` [string]: A description of the run which is written to the output file.
- `write_header` [bool, default=true]: Write a header with information about the run to the output file.
- Redshift-space distortions (RSD):
    - `rsd` [bool, default=false]: Compute the power spectrum multipoles in redshift space.
    - `rsd_growth_f` [float, default=0]: The growth rate $`f`$ used for RSD.
- Biased tracers (only implemented for `rsd = true` and `loops <= 1`):
    - `biased_tracers` [bool, default=false]: Compute the power spectrum of biased tracers.
    - `bias_parameters` [array of 4 floats, default=[1,0,0,0]]: The bias parameters.
- IR resummation:
    - `ir_resum` [bool, default=false]: Perform IR resummation.
    - `k_s` [float, default=0.2]: Separation scale for IR resummation.
    - `k_osc` [float, default=1/110]: BAO oscillation scale.
    - `pt_order` [integer, **required** if `ir_resum = true`]: The perturbation theory order one is working at (avoids overcounting of IR contributions).
- Single hard limit:
    - `single_hard_limit` [bool, default=false]: Compute the single hard limit, where one loop momentum is fixed.
    - `single_hard_limit_q` [float, default=10*q_max]: The fixed loop momentum.
- EFT displacement dispersion:
    - `compute_eft_displacement_dispersion` [bool, default=false]: Compute the small scale displacement dispersion $`\sigma_d^2 = \int_{q_\mathrm{max}}^\infty dq P_\mathrm{lin}(q)`$ (relevant for EFT).
    - `eft_displacement_dispersion_infty` [float, default=10]: Upper limit used in place of $`\infty`$ in the integral above.
- `cuba_settings` [group]: The program uses the Suave routine from the CUBA library for Monte Carlo integration. A subset of the Suave integration settings can be set in the cuba_settings group (the rest is hardcoded and can only be changed in the source file):
    - `abs_tolerance` [float, default=1e-12]: Integration absolute tolerance; sampling stops when absolute or relative tolerance is reached.
    - `rel_tolerance` [float, default=1e-4]: Integration relative tolerance; sampling stops when absolute or relative tolerance is reached.
    - `max_evaluations` [integer, default=1e6]: Maximum number of integrand evaluations, if absolute or relative tolerance is not reached.
    - `verbosity_level` [integer, default=0]: Verbosity level, given as an integer 0-3:
        - 0: Do not print any output
        - 1: Print reasonable information on the integration progress
        - 2: As 1, but also print input parameters.
        - 3: As 2, but also print subregion results.
    - `n_cores` [integer, default=0]: Number of threads to spawn for integration (in addition to master).
    - `statefile` [string]: Store the internal state of the integration at every iteration in a file. Useful for resuming an interrupted integration or splitting long-running integrations into several sequential jobs. To resume an integration, just provide an existing state file. Obviously, restarting integration from an existing state file using a different integrand yields wrong results.
    - `statefile_path` [string]: Same as above, but only setting the path in which the state file is saved/read. The state file is named after `k_a` or `k_a_idx`, depending on which was set. For bispectrum, `k_b` and `k_c` (or the corresponding indices) are also included in name.
    - `retain_statefile` [bool, default=false]: Do not delete the state file when the integration terminated successfully (i.e. reaching either absolute or relative tolerance).
- `ode_settings` [group]: Settings for the GSL ODE solver (only used when `dynamics` is "evolve-asymp-ics" or "evolve-eds-ics"):
    - `abs_tolerance` [float, default=1e-6]: Absolute tolerance.
    - `rel_tolerance` [float, default=1e-4]: Relative tolerance.
    - `start_step` [float, default=1e-3]: Size of first step in solving the ODE numerically.
- For kernel `dynamics` equal to "evolve-asymp-ics" or "evolve-eds-ics", the value of the kernels are stored for every time step in a grid. The following settings define this grid. The time variable used is typically $`\eta = \ln(D)`$, where $`D`$ is the growth factor.
    - `eta_ini` [float, **required**]: Initial time.
    - `eta_fin` [float, **required**]: Final time.
    - `time_steps` [integer, **required**]: Number of time steps.
    - `eta_asymp` [float, **required** for "evolve-asymp-ics"]: Asymptotic time in the past.
    - `pre_time_steps` [integer, **required** for "evolve-asymp-ics"]: Number of time steps between `eta_asymp` and `eta_ini`.
- For kernel `dynamics` equal to "evolve-asymp-ics" or "evolve-eds-ics", the following quantities entering the equations of motion are **required** (see [2008.00013](http://www.arxiv.org/abs/2008.00013)):
    - `kappa_values` [array of floats]: Constant parameters $`\kappa`$.
    - `zeta_files` [array of strings]: Files containing time-dependent functions $`\zeta(\eta)`$.
    - `xi_files` [array of strings]: Files containing time- and scale-dependent functions $`\xi(\eta, k)`$.
- `omega_eigenspace_settings` [group]: For "evolve-asymp-ics", the initial conditions are set by an eigenmode of the linear evolution matrix $`\Omega`$ at `eta_ini`, which is computed on a grid in $`k`$ and interpolated:
    - `eigenmode` [integer, default=0]: Which (zero-indexed) eigenmode to use.
    - `k_min` [float, default=q_min]: Lower limit of the $`k`$-grid.
    - `k_max` [float, default=q_max]: Upper limit of the $`k`$-grid.
    - `N` [integer, default=100]: Number of grid points.
    - `imag_threshold` [float, default=1e-3]: Maximal allowed imaginary part of the eigenvalue/eigenvector components.

See [tests/ini/quijote_Mnu_0p1eV_L1.cfg](tests/ini/quijote_Mnu_0p1eV_L1.cfg) for an example configuration with evolved kernels for a cosmology with massive neutrinos.

## Debugging

To run the program in debug mode, where numerous checks and assertions are
performed at runtime, build the program with `make debug` and run
`debug`. As is to be expected, running the program in debug mode
increases the computation time greatly.

## Benchmarking

The program can be benchmarked by building the benchmark binary with `make
benchmark` and running it as `bench`. This requires the Google benchmark
library being installed.
