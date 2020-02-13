# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

module SimSpin

    export sim_data,
            sim_analysis,
            obs_data_prep,
            build_datacube,
            flux_grid,
            ifu_cube,
            Gaussian_blur,
            Moffat_blur,
            Hernquist,
            NFW,
            IFU,
            SAMI,
            MaNGA,
            CALIFA,
            Hector,
            Environment

    include("types/blur_types.jl")
    include("types/dark_matter_types.jl")
    include("types/environment_type.jl")
    include("types/particle_types.jl")
    include("types/Galaxy_particle.jl")
    include("types/observation_types.jl")
    include("types/Sim_particle.jl")
    include("types/shell_types.jl")
    include("types/telescope_types.jl")

    include("io/sim_data.jl")

    include("celestial/cosdistAngScale.jl")
    include("celestial/cosdistTravelTime.jl")
    include("celestial/cosdistLumDist.jl")
    include("celestial/cosgrowRhoDE.jl")
    include("celestial/e_inv.jl")
    include("celestial/getcos.jl")

    include("ProSpect/data/bc_data.jl")
    include("ProSpect/data/constants.jl")
    include("ProSpect/data/get_filter.jl")
    include("ProSpect/Lum2FluxFactor.jl")
    include("ProSpect/photom_lum.jl")
    include("ProSpect/utilities.jl")

    include("functions/build_datacube.jl")
    include("functions/blur_cube.jl")
    include("functions/flux_grid.jl")
    include("functions/ifu_cube.jl")
    include("functions/obs_data_prep.jl")
    include("functions/sim_analysis.jl")

    include("utilities/ap_shapes.jl")
    include("utilities/assign_flux.jl")
    include("utilities/galaxy_centre.jl")
    include("utilities/galaxy_orient.jl")
    include("utilities/particle_kinematics.jl")
    include("utilities/r_functions.jl")
    include("utilities/shell_analysis.jl")
    include("utilities/sim_to_galaxy.jl")
end
