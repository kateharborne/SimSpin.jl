var documenterSearchIndex = {"docs":
[{"location":"#SimSpin.jl-Documentation-1","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"","category":"section"},{"location":"#Usage-1","page":"SimSpin.jl Documentation","title":"Usage","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"For the installation procedure of the SimSpin package please follow the installation instructions on the package's README.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Once installed, a kinematic datacube for the example galaxy observed from redshift 0.05 and inclined at 30 degrees from face on can be built using the following procedure:","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"> sim_data = SimSpin.sim_data(\"path/to/SimSpin/example/SimSpin_example.hdf5\")\n> datacube = SimSpin.build_datacube(sim_data, z = 0.05, inc_deg = 30)","category":"page"},{"location":"#Functions-1","page":"SimSpin.jl Documentation","title":"Functions","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Blur\nbuild_datacube\nflux_grid\nifu_cube\nobs_data_prep\nsim_data","category":"page"},{"location":"#SimSpin.Blur","page":"SimSpin.jl Documentation","title":"SimSpin.Blur","text":"Blur(psf; sigma, fwhm)\n\nCreate a struct containing the blurring to be used.\nMust specify blurring point spread function, `psf`,\ncurrently only \"Gaussian\" is supported.\n\nKeyword arguments (one must be specified, sigma is prioritised):\n\n    sigma       The standard deviation of the point spread function\n    fwhm        The full width half max of the point spread function\n\n\n\n\n\n","category":"type"},{"location":"#SimSpin.build_datacube","page":"SimSpin.jl Documentation","title":"SimSpin.build_datacube","text":"build_datacube(galaxy_data;\n                r200 = 200,\n                z = 0.05,\n                fov = 15.,\n                ap_shape = \"circular\",\n                central_wvl = 4800,\n                lsf_fwhm = 2.65,\n                pixel_sscale = 0.5,\n                pixel_vscale = 1.04,\n                inc_deg = 70,\n                filter = \"r\",\n                blur = Blur(\"none\"))\n\nReturns a simulated ifu datacube for input, galaxy_data, an array of Particle. Observational parameters specified as keyword arguments.\n\nKeyword arguments (optional):\n\nr200            The virial radius specified in the simulation, kpc.\nz               The projected redshift at which the observation is made.\nfov             The field of view of the IFU, diameter in arcseconds.\nap_shape        The shape of the field of view, with options \"circular\", \"square\" or \"hexagonal\".\ncentral_wvl     The central filter wavelength used for the observation, given in angstroms.\nlsf_fwhm        The line spread function full-width half-max, given in angstroms.\npixel_sscale    The corresponding spatial pixel scale associated with a given telescope output in arcseconds.\npixel_vscale    The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.\ninc_deg         The inclination at which to observe the galaxy in degrees.\nfilter          If particles type is ssp, the filter within which the SED is generated. Options include \"r\" and \"g\"  for SDSS-r and SDSS-g bands respectively.\nblur            Specified to apply observational seeing effects to the cube. Use `Blur()` function to create blur profile. If omitted no blurring occurs.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.flux_grid","page":"SimSpin.jl Documentation","title":"SimSpin.flux_grid","text":"flux_grid(parts_in_cell,\n            ap_region,\n            sbin,\n            vbin,\n            redshift,\n            filter)\n\nComputes the fluxes for each element of the IFU data-cube.\n\nThe purpose of this function is to construct the mock flux values within each cell of the IFU cube. It accepts output parameters from obs_data_prep() and returns a 3D array containing the flux at each cell position due to contributions from the particles. If ssd particles are supplied, an SED is generated in each cell using ProSpect. Else, the luminosity in each cell is converted to flux.\n\nParameters:\n\nparts_in_cell       1D array of the particles corresponding to each element in the IFU data-cube.\nap_region           The aperture region mask used to remove flux outside of the specified aperture.\nsbin                The number of spatial bins in the aperture.\nvbin                The number of velocity bins in the flux grid.\nredshift            The projected redshift at which the observation is made.\nfilter              If ssp particles are supplied, the filter within which the SED is generated.\n                    Options include \"r\" and \"g\"  for SDSS-r and SDSS-g bands respectively.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.ifu_cube","page":"SimSpin.jl Documentation","title":"SimSpin.ifu_cube","text":"ifu_cube(flux_grid,\n            parts_in_cell,\n            sbin,\n            vbin,\n            vseq,\n            lsf_size)\n\nThe purpose of this function is to construct an IFU data cube. It accepts the fluxgrid output by the `fluxgrid()` function and returns a similar, IFU-like, 3D array where each particle's flux contributes a Gaussian distribution in the velocity axis.\n\nParameters:\n\nflux_grid       Flux grid output by `flux_grid()`\nparts_in_cell   1D array of the particles corresponding to each element in the IFU data-cube.\nsbin            The number of spatial bins in the aperture.\nvbin            The number of velocity bins in the flux grid.\nvseq            The bounds of each velocity bin in the flux grid.\nlsf_size        The Gaussian standard deviation of the line spread function in km/s.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.obs_data_prep","page":"SimSpin.jl Documentation","title":"SimSpin.obs_data_prep","text":"obs_data_prep(galaxy_data;\n                r200 = 200,\n                z = 0.05,\n                fov = 15.,\n                ap_shape = \"circular\",\n                central_wvl = 4800,\n                lsf_fwhm = 2.65,\n                pixel_sscale = 0.5,\n                pixel_vscale = 1.04,\n                inc_deg = 70)\n\nThis function prepares the particle data for a given observation.\n\nReturns:\n\ngalaxy_data         Array of particle data formatted for the observation specified by the parameters.\nparts_in_cell       3D array of the particles corresponding to each element in the IFU data-cube.\nap_region           The aperture region mask used to remove flux outside of the specified aperture.\nsbin                The number of spatial bins in the aperture.\nvbin                The number of velocity bins in the flux grid.\nvseq                The bounds of each velocity bin in the flux grid.\nlsf_size            The Gaussian standard deviation of the line spread function in km/s.\nang_size            The angular size given redshift z in kpc\nsbinsize            The spatial bin size in kpc per bin\n\nKeyword arguments (optional):\n\nr200            The virial radius specified in the simulation, kpc.\nz               The projected redshift at which the observation is made.\nfov             The field of view of the IFU, diameter in arcseconds.\nap_shape        The shape of the field of view, with options \"circular\", \"square\" or \"hexagonal\".\ncentral_wvl     The central filter wavelength used for the observation, given in angstroms.\nlsf_fwhm        The line spread function full-width half-max, given in angstroms.\npixel_sscale    The corresponding spatial pixel scale associated with a given telescope output in arcseconds.\npixel_vscale    The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.\ninc_deg         The inclination at which to observe the galaxy in degrees.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.sim_data","page":"SimSpin.jl Documentation","title":"SimSpin.sim_data","text":"sim_data(filename;\n        pytpe = [],\n        ssp = false)\n\nReads in a SimSpin format HDF5 file at location, filename. Returns array of Sim_particles.\n\nKeyword arguments (optional):\n\nptype       A vector of the particles types to be read in e.g. ptype = [1,3].\n            If omitted all particles types will be read.\nssp         Boolean value to use ssp particle information.\n\n\n\n\n\n","category":"function"},{"location":"#Multi-Threading-1","page":"SimSpin.jl Documentation","title":"Multi-Threading","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"The SimSpin package has multi-threading enabled in some critical functions. To use SimSpin with x threads (where x is the integer number of threads desired) you must call","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"export JULIA_NUM_THREADS=x","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"in a Terminal before the Julia REPL is started. This environment variable defaults to 1 if not set before the session has begun.","category":"page"},{"location":"#Data-Input-Format-1","page":"SimSpin.jl Documentation","title":"Data Input Format","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Here we outline the expected file format accepted by SimSpin.  If you would like to generate this file automatically, a short Python function has been written that uses the pynbody package to read in various simulation data types and generate a SimSpin compatible HDF5 file. See create_SimSpinFile.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"If you would rather generate the SimSpin file independently, the expected file format is outlined below.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"> SimSpin_example.hdf5\n\n>> /PartType0           # Each particle type included in the simulation has its own group.\n>>> /PartType0/Mass     # Each group then has a series of data sets assocaited,\n>>> /PartType0/vx       #   including the position, velocity and Mass of each particle.\n>>> /PartType0/vy\n>>> /PartType0/vz\n>>> /PartType0/x\n>>> /PartType0/y\n>>> /PartType0/z\n\n>> /PartType1\n>>> ...","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"We use the same PartType definition as Gadget: PartTypeX where 0 - gas, 1 - dark matter, 2 - disc, 3 - bulge, 4 - stars. For PartType0-3, each PartType group contains the same data sets as above. If the simulation contains stars, the Age and Metallicity information for each particle is also included:","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"> SimSpin_example.hdf5\n>> /PartType4\n>>> /PartType4/Age\n>>> /PartType4/Mass\n>>> /PartType4/Metallicity\n>>> /PartType4/vx        \n>>> /PartType4/vy\n>>> /PartType4/vz\n>>> /PartType4/x\n>>> /PartType4/y\n>>> /PartType4/z","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"If the file is set up in this way, the simulation data can easily be read into the SimSpin package.","category":"page"},{"location":"#References-1","page":"SimSpin.jl Documentation","title":"References","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"A. Pontzen, R Roskar, G. Stinson and R. Woods, (2013), \"pynbody: N-Body/SPH analysis for python\",  Astrophysics Source Code Library, record ascl:1305.002","category":"page"}]
}
