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

## Dependencies

- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) used to solve ODEs numerically
- [The CUBA library](http://www.feynarts.de/cuba/) library for Monte Carlo integration
- [libconfig](https://hyperrealm.github.io/libconfig/) for reading configuration files

The path to include- and library-files of the above dependencies must be
provided to the `-I` and `-L` options to the compiler in the Makefile.

## Building and running the program

The binary can be compiled by running `make` (with the `-j` option to
build in parallel). It is placed in the build directory and can be run as
```
cosPT configuration_file [--option argument]
```
where a configuration file is given as argument (see below). See
`cosPT --help` for a list of command line options that can be given.

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

##### Options that can be set in the configuration file:

- `loops` [integer, **required**]: Number of loops to compute.
- External wavenumber for which the loop correction should be computed (for bispectrum one of the wavenumbers). Either provide `k_a` or both `k_a_grid` and `k_a_idx` (the latter can however be set with a command line option):
    - `k_a` [float]: The value of the external wavenumber given in the same units as the wavenumber grid for the input linear power spectrum.
    - `k_a_grid` [string]: File containing a wavenumber grid where each point is separated by a newline.
    - `k_a_idx` [integer]: The (zero-based) index of the wavenumber in the `k_a_grid` for which the loop corrections should be computed (can also be given as command line option).
- Analogous settings exist for `k_b` and `k_c` for the bispectrum. In stead of `k_c`, the cosine of the angle between the wavevectors `k_a` and `k_b`, `cos_ab`, can be given.
- `q_min` and `q_max` [float, **required**]: The lower and upper integration limits, respectively, with same units as external wavenumbers.
- `spectrum` [string, **required**]: spectrum to compute: "powerspectrum" or "bispectrum".
- `correlations` [list, **required**]: List of correlations (zero-indexed) to compute, each entry given as an array of two or three integers for power spectrum and bispectrum, respectively. Example:
```
correlations = ([0,0], [1,1]) # power spectrum: compute <delta delta> and <theta theta>
correlations = ([0,0,0])      # bispectrum: compute <delta delta delta>
```
- `dynamics` [string, **required**]: The kernel dynamics, either "eds-spt", "evolve-ic-asymp" or "evolve-ic-eds". The first corresponds to EdS-SPT kernels, the second and third are used for dynamically evolved kernels with different initial conditions.
- `input_ps_file` [string, **required**]: The (full or relative path) filename of the input linear power spectrum. The file should contain two columns: the first with the wavenumber grid, and the second with the corresponding values for the power spectrum. Lines starting with \# are ignored.
- `input_ps_rescale` [float or string]: Rescale the input linear power spectrum by a number before interpolation (does not alter the input file). One can use the string "twopi^3" or "twopi^-3" to conventionally rescale it by $`(2\pi)^3`$ or $`(2\pi)^{-3}`$, respectively.
- Output file must be set via either of the following settings:
    - `output_file` [string]: The (full or relative path) filename to write to (either full or relative path).
    - `output_path` [string]: The (full or relative) directory in which the output file should be put. The output file is named after `k_a` or `k_a_idx`, depending on which was set. For bispectrum, `k_b` and `k_c` (or the corresponding indices) are included in the name of the file.
- `description` [string]: A description of the run which is written to the output file.
- `cuba_settings` [group]: The program uses the Suave routine from the CUBA library for Monte Carlo integration. A subset of the Suave integration settings can be set in the cuba_settings group (the rest is hardcoded and can only be changed in the source file):
    - `abs_tolerance` [float, default=1e-12]: Integration absolute tolerance; sampling stops when absolute or relative tolerance is reached.
    - `rel_tolerance` [float, default=1e-4]: Integration relative tolerance; sampling stops when absolute or relative tolerance is reached.
    - `max_evaluations` [integer, default=1e6]: Maximum number of integrand evaluations, if absolute or relative tolerance is not reached.
    - `verbosity_level` [integer, default=1]: Verbosity level, given as an integer 0-3:
        - 0: Do not print any output
        - 1: Print reasonable information on the integration progress
        - 2: As 1, but also print input parameters.
        - 3: As 2, but also print subregion results.
    - `n_cores` [integer, default=0]: Number of threads to spawn for integration (in addition to master).
    - `statefile` [string]: Store the internal state of the integration at every iteration in a file. Useful for resuming an interrupted integration or splitting long-running integrations into several sequential jobs. To resume an integration, just provide an existing state file. Obviously, restarting integration from an existing state file using a different integrand yields wrong results.
    - `statefile_path` [string]: Same as above, but only setting the path in which the state file is saved/read. The state file is named after `k_a` or `k_a_idx`, depending on which was set. For bispectrum, `k_b` and `k_c` (or the corresponding indices) are also included in name.
    - `retain_statefile` [bool, default=false]: Do not delete the state file when the integration terminated successfully (i.e. reaching either absolute or relative tolerance).
- `ode_settings` [group]: Settings for the GSL ODE solver (only used when `dynamics` is "evolve-ic-asymp" or "evolve-ic-eds"):
    - `abs_tolerance` [float, default=1e-6]: Absolute tolerance.
    - `rel_tolerance` [float, default=1e-4]: Relative tolerance.
    - `start_step` [float, default=1e-3]: Size of first step in solving the ODE numerically.
- For kernel `dynamics` equal to "evolve-ic-asymp" or "evolve-ic-eds", the value of the kernels are stored for every time step in a grid. The following settings define this grid. The time variable used is typically $`\eta = \ln(D)`$, where $`D`$ is the growth factor.
    - `eta_ini` [float]: Initial time.
    - `eta_fin` [float]: Final time.
    - `time_steps` [integer]: Number of time steps.
    - `eta_asymp` [float]: For "evolve-ic-asymp", asymptotic time time in the past.
    - `pre_time_steps` [integer]: For "evolve-ic-asymp", number of time steps between `eta_asymp` and `eta_ini`.
- In addition, various files containing certain physical quantities can be provided to be used for dynamical evolution of the kernels. See [src/parameters.cpp](src/parameters.cpp) for info and e.g. [ini/m_nu_0.02eV.cfg](ini/m_nu_0.02eV.cfg) for examples.

## Debugging

To run the program in debug mode, where numerous checks and assertions are
performed at runtime, build the program with `make debug` and run
`debug`. As is to be expected, running the program in debug mode
increases the computation time greatly.

## Benchmarking

The program can be benchmarked by building the benchmark binary with `make
benchmark` and running it as `bench`. This requires the Google benchmark
library being installed.
