# Test directory

This directory contains integration tests which can be run with ```run_tests.sh```.

Directory tree:
```
├── data                  # Data to test against; "truth"
│   ├── eds_spt_bs        # EdS-SPT bispectrum
│   └── eds_spt_ps        # EdS-STP powerspectrum
│   └── quijote_Mnu_0p1eV # Quijote Mnu = 0.1eV power spectrum
├── ini                   # ini-files for tests
├── input                 # input for tests: power spectrum, dynamics files to interpolate
│   └── quijote_Mnu_0p1eV
├── isapprox.jl           # julia-file for reading and comparing columns with uncertainties
├── run_tests.sh          # test script
└── README.md
```
