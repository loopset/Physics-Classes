# Physics Classes

This repository provides several utilities to perform several tasks demanded by a general experiment:
* Calibrate any spectra using the `Calibration` namespace classes
* Store experiment information
* Provides a set of color palettes well-suited for CVD
* Utils to fit an excitation energy spectrum using least-squares method (log-likelihood to be implemented)
* Classes to perform differential cross-sections calculations and comparisons with DWBA calculations
* Set of utilities to implement efficiencies, sigma interpolators, etc.

## Installation

Just clone and compile the project using CMake and `make -jN install`. Then, **source the thisPhysicsClasses.sh** in your `.bashrc` to get this working with ROOT.
