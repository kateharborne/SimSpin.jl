# SimSpin.jl Documentation

## Usage

For the installation procedure of the SimSpin package please follow the installation instructions on the package's [*README*](https://github.com/kateharborne/SimSpin.jl).

### Mock observables of an N-body simulation
Once installed, only five steps are required to take an observation, generate a mock-IFU datacube and then export it to a FITS file:
1.  Construct a `Telescope` object. This specifies the field of view to be used, the aperture shape, etc. In this example we will use the default [`SAMI`](@ref) telescope constructor. See [Telescope Constructors](@ref) for other default and customisable constructors.

    ```
        > telescope = SAMI()
    ```

2.  Construct an environment in which the observation is taken. This specifies the redshift of the galaxy, its inclination, the virial radius, the mass to light ratio and the seeing conditions respectively. We will set redshift to be 0.05, inclination to be 70 degrees, virial radius to be 200 kpc, to mass to light ratio to be 1 and use no blurring. See [Environment Constructor](@ref) and [Blur Constructors](@ref) for more details.

    ```
        > environment = Environment(0.05, 70, 200, 1.)
    ```

3.  Read in a simulation's particle data. In this example we will use the example file embedded in SimSpin. Usually you will need to pass [`sim_data`](@ref) the file path to your HDF5 data file. See [Data Import](@ref) for more details.

    ```
        > data = sim_data()
    ```

4.  Build the datacube as a combination of the galaxy particle data, a telescope and an environment. It returns a mock-IFU datacube and a summary of the observational properties used. See [`build_datacube`](@ref) for more details.

    ```
        > datacube, observe = build_datacube(data, telescope, environment)
    ```

5.  Export the datacube to a FITS file for viewing. See [Data Export](@ref) for more details

    ```
        > sim_FITS(datacube, observe, "SimSpin_Example_Observation.fits")
    ```

!!! note
    The most efficient way to take multiple observations is to pass [`build_datacube`](@ref) an array of Environments. These can be created individually or by passing the [Environment Constructor](@ref) an array of the values that you wish to observe. For example, to take an observation at redshifts 0.05, 0.1 and 0.15 each at an inclination of 70 degrees and 80 degrees using a constant virial radius of 200 kpc there are two options:

    ```
        > environments = Environment([0.05; 0.1; 015], [70; 80], 200)
        > environments = Environment(0.05:0.05:0.15, 70:10:80, 200)
    ```

    Both options will return an array of environments which can then be used as the environment parameter in Step 4 above. This will then return an array of tuples, each one consisting of the datacube and the observational properties used. See [`build_datacube`](@ref) and [Environment Constructor](@ref) for more details.

### Mock observables of a hydro-dynamical simulation
To use SimSpin with a hydro-dynamical simulation a similar process is followed with notable, minor differences:

1.  Construct a `Telescope` object. We must specify a filter when using a hydro-dynamical simulation. See [Telescope Constructors](@ref) for more options.

    ```
        > telescope = SAMI(filter="g")
    ```

2.  Construct an environment in which the observation is taken. This is the same as for an n-body simulation but does not require a mass to light ratio. See [Environment Constructor](@ref) and [Blur Constructors](@ref) for more details.

    ```
        > environment = Environment(0.05, 70, 200)
    ```

3.  Read in a simulation's particle data. For hydro-dynamical simulations we must specify that we want SSP to be used. The default is to not use SSP data so if SSP is not actively set to true a hydro-dynamical simulation will still be observed but a mass to light ratio will be used instead of spectra. See [Data Import](@ref) for more details.

    ```
        > data = sim_data("your/filename.hdf5", ssp=true)
    ```

4.  Build the datacube as a combination of the galaxy particle data, a telescope and an environment. It returns a mock-IFU datacube and a summary of the observational properties used. See [`build_datacube`](@ref) for more details.

    ```
        > datacube, observe = build_datacube(data, telescope, environment)
    ```

5.  Export the datacube to a FITS file for viewing. See [Data Export](@ref) for more details

    ```
        > sim_FITS(datacube, observe, "SimSpin_Example_Observation.fits")
    ```

!!! warning
    Do not multithread any SimSpin functions as a user. Each observation is already multithreaded and will run on as many threads as available. To see how to allow SimSpin to use more threads see [Multi-Threading](@ref).

## General Functions
```@docs
build_datacube
flux_grid
ifu_cube
obs_data_prep
sim_to_galaxy
```

## Constructors
### Telescope Constructors

An `IFU` object denotes all the parameters required to make a mock observation.
Any generic IFU can be made using the `IFU()` function.
Default constructors can also be used to emulate famous IFU survey instruments.
Currently `SAMI`, `MaNGA`, `CALIFA` and `Hector` are supported.
```@docs
IFU
SAMI
MaNGA
MUSE
CALIFA
Hector
```

### Environment Constructor
```@docs
Environment
```

### Blur Constructors
```@docs
Gaussian_blur
Moffat_blur
```

## Multi-Threading

The SimSpin package has multi-threading enabled in some critical functions.
To use SimSpin with more than one thread the only thing you need to do is start Julia with the desired number of threads. Instructions on how to do this can be found [here](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads-1).

!!! warning
    Do not manually multithread any SimSpin functions as a user. Each observation is already multithreaded and will run on as many threads as available.

## Data Import

To read in a SimSpin HDF5 file as defined below the [`sim_data`](@ref) function must be used.

```@docs
sim_data
```

The expected file format accepted by SimSpin is outlined below.  If you would like to generate this file automatically, a short Python function has been written that uses the [pynbody](https://github.com/pynbody/pynbody) package to read in various simulation data types and generate a SimSpin compatible HDF5 file. See [create_SimSpinFile](https://github.com/kateharborne/create_SimSpinFile).

If you would rather generate the SimSpin file independently, the expected file format is outlined below.

```
> SimSpin_example.hdf5

>> /PartType0           # Each particle type included in the simulation has its own group.
>>> /PartType0/Mass     # Each group then has a series of data sets assocaited,
>>> /PartType0/vx       #   including the position, velocity and Mass of each particle.
>>> /PartType0/vy
>>> /PartType0/vz
>>> /PartType0/x
>>> /PartType0/y
>>> /PartType0/z

>> /PartType1
>>> ...
```
We use the same PartType definition as Gadget: PartTypeX where 0 - gas, 1 - dark matter, 2 - disc, 3 - bulge, 4 - stars. For PartType0-3, each PartType group contains the same data sets as above. If the simulation contains stars, the Age and Metallicity information for each particle is also included:

```
> SimSpin_example.hdf5
>> /PartType4
>>> /PartType4/Age
>>> /PartType4/Mass
>>> /PartType4/Metallicity
>>> /PartType4/vx        
>>> /PartType4/vy
>>> /PartType4/vz
>>> /PartType4/x
>>> /PartType4/y
>>> /PartType4/z
```
Optionally SSP Star particles can also be defined with initial mass and stellar formation time:

```
> SimSpin_example.hdf5
>> /PartType4
>>> /PartType4/InitialMass
>>> /PartType4/Mass
>>> /PartType4/Metallicity
>>> /PartType4/StellarFormationTime
>>> /PartType4/vx        
>>> /PartType4/vy
>>> /PartType4/vz
>>> /PartType4/x
>>> /PartType4/y
>>> /PartType4/z
```
If the file is set up in this way, the simulation data can easily be read into the SimSpin package.

### References

A. Pontzen, R Roskar, G. Stinson and R. Woods, (2013), "pynbody: N-Body/SPH analysis for python",  Astrophysics Source Code Library, record ascl:1305.002

## Data Export

Currently only FITS file exports are supported.

### FITS Export
```@docs
sim_FITS
```
