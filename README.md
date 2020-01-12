# SimSpin-Julia

This is repository for the pure Julia implementation of Katherine Harborne's astronomy package, SimSpin.

To see the original R implementation see Kate's GitHub repo [here](https://github.com/kateharborne/SimSpin).

[![Build Status](https://travis-ci.com/kateharborne/SimSpin-Julia.svg?branch=master)](https://travis-ci.com/kateharborne/SimSpin-Julia) |

## Installation

This package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add https://github.com/kateharborne/SimSpin-Julia
```

Once the package has been added leave the REPL's Pkg mode by pressing `backspace` and run:

```
> using SimSpin
```

All user exported functions should now be available in your Julia REPL.

## SimSpin

SimSpin - A package for the kinematic analysis of galaxy simulations

The purpose of the Simspin package is to take a galaxy simulation and measure an the kinematics of that model as if it had been observed using an IFU. A kinematic data cube can be produced using the functions in this package; from this cube, "observables" can be measured. Specifically, the observable spin parameter, ``\lambda_r``. This package, once installed, is fully documented and tested.

Simulation data needs to be laid out in HDF5 format for processing with SimSpin. If using a Gadget binary or HDF5 simulation output, see [create_SimSpinFile](https://github.com/kateharborne/create_SimSpinFile) for automatic generation of SimSpin compatible files.  If you would rather generate the SimSpin file independently, the expected format is outlined below. This SimSpin file can then be read into SimSpin using the function:

```
galaxy_data = SimSpin::sim_data(filename = SimSpin_example.hdf5)
```

This function produces a list that can be accessed by each of the basic SimSpin analysis functions listed below. While it is possible to use each function in the package in order and examine the output at each stage, there are three basic analysis functions designed to give the data in a user friendly format. We suggest first using these functions:

1. `sim_analysis()` - This function is designed to output the kinematic properties of the galaxy model to be observed. This provides the comparison to the kinematic observables produced in the following functions.

2. `build_datacube()` - This function produces the kinematic data cube prior to kinematic analysis. This allows the user to take the cubes to use in some other form of analysis without having to calculate ``\lambda_r``.

3. `find_kinematics()` - This function produces a kinematic data cube and calculates the observed spin parameter, ellipticity, inclination and the corresponding flux, line-of-sight velocity and line-of-sight velocity dispersion images.

All user-exported functions are explained in greater detail below.

By varying the effects of observational seeing, the measurement radius, projected inclination and distance, and the telescope parameters within the `find_lambda()` function, we can begin to understand how inherent limitations of observing galaxies can effect the measurement of ``\lambda_r`` by comparing to the true spin parameter than is measured in the sim_analysis() function.

## Documentation

- [**STABLE**](https://kateharborne.github.io/SimSpin-Julia/stable) &mdash; **documentation in progress.**
- [**DEVEL**](https://kateharborne.github.io/SimSpin-Julia/dev) &mdash; **Development doumentation.**

### References

A. Pontzen, R Roskar, G. Stinson and R. Woods, (2013), "pynbody: N-Body/SPH analysis for python",  Astrophysics Source Code Library, record ascl:1305.002
