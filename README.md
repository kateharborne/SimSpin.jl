# SimSpin.jl

[![Build Status](https://travis-ci.com/kateharborne/SimSpin.jl.svg?branch=master)](https://travis-ci.com/kateharborne/SimSpin.jl)

This is repository for the pure Julia implementation of [Katherine Harborne's](https://github.com/kateharborne) astronomy package, SimSpin.

To see the original R implementation see Kate's GitHub repo [here](https://github.com/kateharborne/SimSpin).

## Installation

If you are new to Julia please see the [Julia `Getting Started`](https://docs.julialang.org/en/v1/manual/getting-started/) page for instructions on how to install Julia.

Once Julia is running successfully this package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add https://github.com/kateharborne/SimSpin.jl
```

Once the package has been added leave the REPL's Pkg mode by pressing `backspace` and run:

```
> using SimSpin
```

After pre-compilation all user exported functions should be available in your Julia REPL.

## SimSpin

SimSpin - A package for the kinematic analysis of galaxy simulations

The purpose of the Simspin package is to take a galaxy simulation and measure an the kinematics of that model as if it had been observed using an IFU. A kinematic data cube can be produced using the functions in this package; from this cube, "observables" can be measured. Specifically, the observable spin parameter, &#955;<sub>R</sub>. This package, once installed, is fully documented and tested.

By varying the effects of observational seeing, the measurement radius, projected inclination, projected distance and other telescope parameters we can begin to understand some of the inherent effects and limitations of real-world galaxy observations.

## Documentation

- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kateharborne.github.io/SimSpin.jl/stable) &mdash; **Documentation in progress.**
- [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kateharborne.github.io/SimSpin.jl/dev) &mdash; **Development documentation.**
