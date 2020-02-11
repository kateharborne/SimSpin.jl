# SimSpin.jl Documentation

## Usage

For the installation procedure of the SimSpin package please follow the installation instructions on the package's [*README*](https://github.com/kateharborne/SimSpin-Julia).

Once installed, a simple procedure of four steps is required to take an observation and generate a datacube:
1.  Create a telescope object. This specifies the field of view to be used, the aperture shape, etc. In this example we will use the default SAMI telescope constructor. Alternatively, the IFU() constructor can be used to create any generic IFU.

    ```
        > telescope = SimSpin.SAMI()
    ```

2.  Read in the simulation's particle data

    ```
        > sim_data = SimSpin.sim_data("path/to/SimSpin/example/SimSpin_example.hdf5")
    ```

3.  Create an observation object. This specifies the observational redshift, inclination and virial radius of the galaxy.

    ```
        > observation = SimSpin.Observation(0.05, 70, 200)
    ```

4.  Build the datacube as a combination of a telescope, an observation and the galaxy particle data.

    ```
        > datacube = SimSpin.build_datacube(sim_data, observation, telescope)
    ```

## Functions
```@docs
build_datacube
flux_grid
ifu_cube
obs_data_prep
sim_data
```
## Constructors
### Telescope Constructors
```@docs
IFU
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
