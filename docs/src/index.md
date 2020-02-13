# SimSpin.jl Documentation

## Usage

For the installation procedure of the SimSpin package please follow the installation instructions on the package's [*README*](https://github.com/kateharborne/SimSpin.jl).

Once installed, a simple procedure of four steps is required to take an observation and generate a datacube:
1.  Construct a `Telescope` object. This specifies the field of view to be used, the aperture shape, etc. In this example we will use the default [`SAMI`](@ref) telescope constructor. See [Telescope Constructors](@ref) for other default and customisable constructors.

    ```
        > telescope = SimSpin.SAMI()
    ```

2.  Read in the simulation's particle data

    ```
        > sim_data = SimSpin.sim_data("path/to/SimSpin/example/SimSpin_example.hdf5")
    ```

3.  Construct an environment in which the observation is taken. This specifies the redshift of the galaxy, its inclination, the virial radius, the mass to light ratio and the seeing conditions respectively. We will set redshift to be 0.05, inclination to be 70 degrees, virial radius to be 200 kpc, to mass to light ratio to be 1 and use no blurring. See [Environment Constructor](@ref) and [Blur Constructors](@ref) for more details.

    ```
        > environment = SimSpin.Environment(0.05, 70, 200, 1.)
    ```

4.  Build the datacube as a combination of the galaxy particle data, a telescope and an environment. This function also returns a summary of the observational properties used.

    ```
        > datacube, observe = SimSpin.build_datacube(sim_data, telescope, environment)
    ```

5.  Export the datacube to a FITS file for viewing. See [Data Export](@ref) for more details

    ```
        > SimSpin.sim_FITS(data_cube, observe, "SimSpin_Example_Observation.fits")
    ```

## Functions
```@docs
build_datacube
flux_grid
ifu_cube
obs_data_prep
sim_data
sim_FITS
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

### Dark Matter Constructors
```@docs
Hernquist
NFW
```

## Multi-Threading

The SimSpin package has multi-threading enabled in some critical functions.
To use SimSpin with `x` threads (where x is the integer number of threads desired) you must call

```
export JULIA_NUM_THREADS=x
```
in a Terminal before the Julia REPL is started. This environment variable defaults to 1 if not set before the session has begun.

## Data Input Format

Here we outline the expected file format accepted by SimSpin.  If you would like to generate this file automatically, a short Python function has been written that uses the [pynbody](https://github.com/pynbody/pynbody) package to read in various simulation data types and generate a SimSpin compatible HDF5 file. See [create_SimSpinFile](https://github.com/kateharborne/create_SimSpinFile).

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
If the file is set up in this way, the simulation data can easily be read into the SimSpin package.

### References

A. Pontzen, R Roskar, G. Stinson and R. Woods, (2013), "pynbody: N-Body/SPH analysis for python",  Astrophysics Source Code Library, record ascl:1305.002

## Data Export

Currently only FITS file exports are supported.

### FITS Export
```@docs
sim_FITS
```
