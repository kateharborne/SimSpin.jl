# SimSpin.jl Documentation

## Usage

For the installation procedure of the SimSpin package please follow the installation instructions on the package's [*README*](https://github.com/kateharborne/SimSpin-Julia/tree/master/README.md).

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
