# SimSpin.jl

[![Build Status](https://github.com/kateharborne/SimSpin.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/kateharborne/SimSpin.jl/actions) [![Coverage Status](https://coveralls.io/repos/github/kateharborne/SimSpin.jl/badge.svg?branch=master)](https://coveralls.io/github/kateharborne/SimSpin.jl?branch=master)

This is repository for the pure Julia implementation of [Katherine Harborne's](https://github.com/kateharborne) astronomy package, SimSpin.

The GitHub repo for Kate's original R implementation can be found [here](https://github.com/kateharborne/SimSpin).

## Installation

If you are new to Julia please see the [Julia `Getting Started`](https://docs.julialang.org/en/v1/manual/getting-started/) page for instructions on how to install Julia.

Once Julia is running successfully this package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add SimSpin
```

Once the package has been added leave the REPL's Pkg mode by pressing `backspace` and run:

```
> using SimSpin
```

After pre-compilation all user exported functions should be available in your Julia REPL.

## SimSpin
#### A package for the kinematic analysis of galaxy simulations

The purpose of the SimSpin package is to take a galaxy simulation and measure the kinematics of that model as if it had been observed using an IFU. A kinematic data cube can be produced using the functions in this package from which "observables" can be measured. This package is fully documented, published with the Astronomical Society of Australia and registered on the Julia Registry.

By varying the effects of observational seeing, the measurement radius, projected inclination, projected distance and other telescope parameters we can begin to understand some of the inherent effects and limitations of real-world galaxy observations.

## Documentation

- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kateharborne.github.io/SimSpin.jl/stable) &mdash; **Documentation for the most recent release.**
- [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kateharborne.github.io/SimSpin.jl/dev) &mdash; **Development documentation.**

## Citation
If you use this code in any of your own published research, please make sure to include the following citation in your bibliography:

K.E. Harborne, C.Power and A.S.G. Robotham, (2020), ["SIMSPIN - Constructing mock IFS kinematic data cubes"](https://ui.adsabs.harvard.edu/abs/2020PASA...37...16H/abstract), Publications of the Astronomical Society of Australia, Volume 37, article id. e016

K.E. Harborne, (2019), ["SimSpin: Kinematic analysis of galaxy simulations"](https://ui.adsabs.harvard.edu/abs/2019ascl.soft03006H/abstract), Astrophysics Source Code Library, record ascl:1903.006
