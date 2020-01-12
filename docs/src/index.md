# SimSpin.jl Documentation

## Usage

For the installation procedure of the SimSpin package please follow the installation instructions on the package's [*README*](https://github.com/kateharborne/SimSpin-Julia).

Once installed, a kinematic datacube for the example galaxy observed from redshift 0.05 and inclined at 30 degrees from face on can be built using the following procedure:

```
> sim_data = SimSpin.sim_data("path/to/SimSpin/example/SimSpin_example.hdf5")
> datacube = SimSpin.build_datacube(sim_data, z = 0.05, inc_deg = 30)
```

## Functions
```@docs
sim_data
obs_data_prep
build_datacube
flux_grid
```

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
