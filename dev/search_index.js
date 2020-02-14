var documenterSearchIndex = {"docs":
[{"location":"#SimSpin.jl-Documentation-1","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"","category":"section"},{"location":"#Usage-1","page":"SimSpin.jl Documentation","title":"Usage","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"For the installation procedure of the SimSpin package please follow the installation instructions on the package's README.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Once installed, a simple procedure of four steps is required to take an observation and generate a datacube:","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Construct a Telescope object. This specifies the field of view to be used, the aperture shape, etc. In this example we will use the default SAMI telescope constructor. See Telescope Constructors for other default and customisable constructors.\n> telescope = SimSpin.SAMI()\nConstruct an environment in which the observation is taken. This specifies the redshift of the galaxy, its inclination, the virial radius, the mass to light ratio and the seeing conditions respectively. We will set redshift to be 0.05, inclination to be 70 degrees, virial radius to be 200 kpc, to mass to light ratio to be 1 and use no blurring. See Environment Constructor and Blur Constructors for more details.\n> environment = SimSpin.Environment(0.05, 70, 200, 1.)\nRead in the simulation's particle data\n> sim_data = SimSpin.sim_data(\"path/to/SimSpin/example/SimSpin_example.hdf5\")\nBuild the datacube as a combination of the galaxy particle data, a telescope and an environment. This function also returns a summary of the observational properties used.\n> datacube, observe = SimSpin.build_datacube(sim_data, telescope, environment)\nExport the datacube to a FITS file for viewing. See Data Export for more details\n> SimSpin.sim_FITS(data_cube, observe, \"SimSpin_Example_Observation.fits\")","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"note: Note\nIn the above procedure we read in our simulation data from file (Step 3) and immediately use it to build a kinematic datacube (Step 4). Although this is perfectly functional for a single observation, if taking more than one observation the build_datacube function must convert the particle data from the simulation reference frame to the galaxy's reference frame every time. The function sim_to_galaxy is provided to do this conversion so that preconverted data can be supplied to build_datacube instead. Run this command before passing galaxy_data to build_datacube for massively increased performance:> galaxy_data = SimSpin.sim_to_galaxy(sim_data)","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"warning: Warning\nDo not multithread any SimSpin functions as a user. Each observation is already multithreaded and will run on as many cores as available. To see how to allow SimSpin to use more cores see Multi-Threading.","category":"page"},{"location":"#General-Functions-1","page":"SimSpin.jl Documentation","title":"General Functions","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"build_datacube\nflux_grid\nifu_cube\nobs_data_prep\nsim_to_galaxy","category":"page"},{"location":"#SimSpin.build_datacube","page":"SimSpin.jl Documentation","title":"SimSpin.build_datacube","text":"build_datacube(galaxy_data, ifu, envir)\n\nReturns a simulated ifu datacube for input, galaxy_data, an array of Particle and a summary of observation properties used. The instrument to be used for observation is given in the Telescope class parameter and the environmental variables (redshift, inclination, seeing, etc) are given by the Environment class parameter.\n\nParameters:\n\ngalaxy_data         Array of `Particle` describing galaxy\nifu                 Struct of type `Telescope`\nenvir               Struct of type `Environment`\n\nReturns:\n\ndata_cube           3D array of a simulated IFU datacube with spatial and velocity binned fluxes\nobserve             A struct of type `Observation` summarising the observational properties used.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.flux_grid","page":"SimSpin.jl Documentation","title":"SimSpin.flux_grid","text":"flux_grid(parts_in_cell,\n            observe,\n            filter)\n\nComputes the fluxes for each element of the IFU data-cube.\n\nThe purpose of this function is to construct the mock flux values within each cell of the IFU cube. It accepts output parameters from obs_data_prep() and returns a 3D array containing the flux at each cell position due to contributions from the particles. If ssd particles are supplied, an SED is generated in each cell using ProSpect. Else, the luminosity in each cell is converted to flux.\n\nParameters:\n\nparts_in_cell       1D array of the particles corresponding to each element in the IFU data-cube.\nobserve             Struct of type `Observation` containing all observation parameters.\nfilter              If ssp particles are supplied, the filter within which the SED is generated.\n                    Options include \"r\" and \"g\"  for SDSS-r and SDSS-g bands respectively.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.ifu_cube","page":"SimSpin.jl Documentation","title":"SimSpin.ifu_cube","text":"ifu_cube(flux_grid,\n            parts_in_cell,\n            observe)\n\nThe purpose of this function is to construct an IFU data cube. It accepts a flux grid in the format output by the flux_grid() function and returns a similar, IFU-like, 3D array where each particle's flux contributes a Gaussian distribution in the velocity axis.\n\nParameters:\n\nflux_grid       Flux grid output by `flux_grid()`\nparts_in_cell   1D array of the particles corresponding to each element in the IFU data-cube.\nobserve         Struct of type `Observation` containing all observation parameters.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.obs_data_prep","page":"SimSpin.jl Documentation","title":"SimSpin.obs_data_prep","text":"obs_data_prep(galaxy_data, ifu, envir)\n\nThis function prepares the particle data for a given observation with the given telescope.\n\nParameters:\n\ngalaxy_data         Array of `Particle` describing galaxy\nifu                 Struct of type `Telescope`\nenvir               Struct of type `Environment`\n\nReturns:\n\ngalaxy_data         Array of particle data formatted for the observation specified by the parameters.\nparts_in_cell       3D array of the particles corresponding to each element in the IFU data-cube.\nobserve             Struct of type `Observation` containing all observation parameters.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.sim_to_galaxy","page":"SimSpin.jl Documentation","title":"SimSpin.sim_to_galaxy","text":"sim_to_galaxy(sim_data)\n\nReturns particle data wrt centre of galaxy in both spatial and velocity space instead of wrt arbitrary point in simulation.\n\n\n\n\n\n","category":"function"},{"location":"#Constructors-1","page":"SimSpin.jl Documentation","title":"Constructors","text":"","category":"section"},{"location":"#Telescope-Constructors-1","page":"SimSpin.jl Documentation","title":"Telescope Constructors","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"An IFU object denotes all the parameters required to make a mock observation. Any generic IFU can be made using the IFU() function. Default constructors can also be used to emulate famous IFU survey instruments. Currently SAMI, MaNGA, CALIFA and Hector are supported.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"IFU\nSAMI\nMaNGA\nCALIFA\nHector","category":"page"},{"location":"#SimSpin.IFU","page":"SimSpin.jl Documentation","title":"SimSpin.IFU","text":"IFU(fov, ap_shape, pixel_sscale, pixel_vscale, central_wvl, lsf_fwhm, filter)\n\nCreates a customisable, mock IFU telescope which can be used to \"observe\" simulations.\n\nParameters:\n\nfov             The field of view of the IFU, diameter in arcseconds.\nap_shape        The shape of the field of view, with options \"circular\", \"square\" or \"hexagonal\".\ncentral_wvl     The central filter wavelength used for the observation, given in angstroms.\nlsf_fwhm        The line spread function full-width half-max, given in angstroms.\npixel_sscale    The corresponding spatial pixel scale associated with a given telescope output in arcseconds.\npixel_vscale    The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.\nfilter          Optional. If particles type is ssp, the filter within which the SED is generated. Options include \"r\" and \"g\"  for SDSS-r and SDSS-g bands respectively.\n\n\n\n\n\n","category":"type"},{"location":"#SimSpin.SAMI","page":"SimSpin.jl Documentation","title":"SimSpin.SAMI","text":"SAMI(;filter)\n\nCreates an IFU using parameters of the SAMI survey. Optional filters \"r\" or \"g\" can also be specified.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.MaNGA","page":"SimSpin.jl Documentation","title":"SimSpin.MaNGA","text":"MaNGA()\n\nCreates an IFU using parameters of the MaNGA survey. Optional filters \"r\" or \"g\" can also be specified.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.CALIFA","page":"SimSpin.jl Documentation","title":"SimSpin.CALIFA","text":"CALIFA()\n\nCreates an IFU using parameters of the CALIFA survey. Optional filters \"r\" or \"g\" can also be specified.\n\n\n\n\n\n","category":"function"},{"location":"#SimSpin.Hector","page":"SimSpin.jl Documentation","title":"SimSpin.Hector","text":"Hector()\n\nCreates an IFU using parameters of the CALIFA survey. Optional filters \"r\" or \"g\" can also be specified.\n\n\n\n\n\n","category":"function"},{"location":"#Environment-Constructor-1","page":"SimSpin.jl Documentation","title":"Environment Constructor","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Environment","category":"page"},{"location":"#SimSpin.Environment","page":"SimSpin.jl Documentation","title":"SimSpin.Environment","text":"Environment(z, inc_deg, r200, mass2light, blur)\n\nCreates a struct containing environmental parameters required for a mock observation of a simulated galaxy.\n\nParameters:\n\nz               The projected redshift at which the observation is made.\ninc_deg         The inclination at which to observe the galaxy in degrees. Relative to face on, rotated around semi-major axis.\nr200            The virial radius specified in the simulation, kpc.\nmass2light      Optional. The mass to light ratio for non-ssp, luminous particles. Defaults to 1.\nblur            Optional. Struct of type `Blur` containing seeing information. If ommitted no blurring is used.\n\n\n\n\n\n","category":"type"},{"location":"#Blur-Constructors-1","page":"SimSpin.jl Documentation","title":"Blur Constructors","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Gaussian_blur\nMoffat_blur","category":"page"},{"location":"#SimSpin.Gaussian_blur","page":"SimSpin.jl Documentation","title":"SimSpin.Gaussian_blur","text":"Gaussian_blur(;sigma, fwhm)\n\nCreate a struct containing seeing information.\n\nKeyword arguments (at least one must be specified, sigma is prioritised):\n\nsigma       The standard deviation of the point spread function\nfwhm        The full width half max of the point spread function\n\n\n\n\n\n","category":"type"},{"location":"#SimSpin.Moffat_blur","page":"SimSpin.jl Documentation","title":"SimSpin.Moffat_blur","text":"Moffat_blur(β;\n            α,\n            fwhm)\n\nCreate a struct containing seeing information. β and either α or fwhm must be specified. If both α and fwhm are specified, α is prioritised.\n\nArguments:\n\nβ           The power component in the Moffat distribution\nα           The core width of the Moffat distribution (optional)\nfwhm        The full width half max of the Moffat distribution (optional)\n\n\n\n\n\n","category":"type"},{"location":"#Dark-Matter-Constructors-1","page":"SimSpin.jl Documentation","title":"Dark Matter Constructors","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Hernquist\nNFW","category":"page"},{"location":"#SimSpin.Hernquist","page":"SimSpin.jl Documentation","title":"SimSpin.Hernquist","text":"Hernquist(;dm_mass = 186.9,\n            dm_a = 34.5)\n\nCreates a Hernquist analytic dark matter mass potential. Default values are shown above.\n\n\n\n\n\n","category":"type"},{"location":"#SimSpin.NFW","page":"SimSpin.jl Documentation","title":"SimSpin.NFW","text":"NFW(;dm_vm = 186.9,\n        dm_a = 34.5,\n        dm_rhof = 0.035)\n\nCreates an NFW analytic dark matter mass potential. Default values are shown above.\n\n\n\n\n\n","category":"type"},{"location":"#Multi-Threading-1","page":"SimSpin.jl Documentation","title":"Multi-Threading","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"The SimSpin package has multi-threading enabled in some critical functions. To use SimSpin with x threads (where x is the integer number of threads desired) you must call","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"export JULIA_NUM_THREADS=x","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"in a Terminal before the Julia REPL is started. This environment variable defaults to 1 if not set before the session has begun.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"warning: Warning\nDo not multithread any SimSpin functions as a user. Each observation is already multithreaded and will run on as many cores as available.","category":"page"},{"location":"#Data-Import-1","page":"SimSpin.jl Documentation","title":"Data Import","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"To read in a SimSpin HDF5 file as defined below the sim_data function must be used.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"sim_data","category":"page"},{"location":"#SimSpin.sim_data","page":"SimSpin.jl Documentation","title":"SimSpin.sim_data","text":"sim_data(filename;\n        pytpe = [],\n        ssp = false)\n\nReads in a SimSpin format HDF5 file at location, filename. Returns array of Sim_particles.\n\nKeyword arguments (optional):\n\nptype       A vector of the particles types to be read in e.g. ptype = [1,3].\n            If omitted all particles types will be read.\nssp         Boolean value to use ssp particle information.\n\n\n\n\n\n","category":"function"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"The expected file format accepted by SimSpin is outlined below.  If you would like to generate this file automatically, a short Python function has been written that uses the pynbody package to read in various simulation data types and generate a SimSpin compatible HDF5 file. See create_SimSpinFile.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"If you would rather generate the SimSpin file independently, the expected file format is outlined below.","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"> SimSpin_example.hdf5\n\n>> /PartType0           # Each particle type included in the simulation has its own group.\n>>> /PartType0/Mass     # Each group then has a series of data sets assocaited,\n>>> /PartType0/vx       #   including the position, velocity and Mass of each particle.\n>>> /PartType0/vy\n>>> /PartType0/vz\n>>> /PartType0/x\n>>> /PartType0/y\n>>> /PartType0/z\n\n>> /PartType1\n>>> ...","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"We use the same PartType definition as Gadget: PartTypeX where 0 - gas, 1 - dark matter, 2 - disc, 3 - bulge, 4 - stars. For PartType0-3, each PartType group contains the same data sets as above. If the simulation contains stars, the Age and Metallicity information for each particle is also included:","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"> SimSpin_example.hdf5\n>> /PartType4\n>>> /PartType4/Age\n>>> /PartType4/Mass\n>>> /PartType4/Metallicity\n>>> /PartType4/vx        \n>>> /PartType4/vy\n>>> /PartType4/vz\n>>> /PartType4/x\n>>> /PartType4/y\n>>> /PartType4/z","category":"page"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"If the file is set up in this way, the simulation data can easily be read into the SimSpin package.","category":"page"},{"location":"#References-1","page":"SimSpin.jl Documentation","title":"References","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"A. Pontzen, R Roskar, G. Stinson and R. Woods, (2013), \"pynbody: N-Body/SPH analysis for python\",  Astrophysics Source Code Library, record ascl:1305.002","category":"page"},{"location":"#Data-Export-1","page":"SimSpin.jl Documentation","title":"Data Export","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"Currently only FITS file exports are supported.","category":"page"},{"location":"#FITS-Export-1","page":"SimSpin.jl Documentation","title":"FITS Export","text":"","category":"section"},{"location":"#","page":"SimSpin.jl Documentation","title":"SimSpin.jl Documentation","text":"sim_FITS","category":"page"},{"location":"#SimSpin.sim_FITS","page":"SimSpin.jl Documentation","title":"SimSpin.sim_FITS","text":"sim_FITS(data_cube,\n            observe,\n            out_file)\n\nOutputs FITS file containing 3D, velocity binned datacube.\n\nParameters:\n\ndata_cube   3D array of a simulated IFU datacube with spatial and velocity binned fluxes\nobserve     A struct of type `Observation` summarising the observational properties used.\nout_file    String denoting path and name of FITS file to be output.\n\n\n\n\n\n","category":"function"}]
}
