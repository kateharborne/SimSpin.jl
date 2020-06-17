# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

module SimSpin

    export sim_data,
            sim_FITS,
            sim_to_galaxy,
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
    include("types/telescope_types.jl")

    include("io/sim_data.jl")
    include("io/sim_FITS.jl")

    include("celestial/celestial.jl")
    using .celestial

    include("ProFit/ProFit.jl")
    using .ProFit

    include("ProSpect/ProSpect.jl")
    using .ProSpect

    include("functions/build_datacube.jl")
    include("functions/blur_cube.jl")
    include("functions/flux_grid.jl")
    include("functions/ifu_cube.jl")
    include("functions/obs_data_prep.jl")

    include("utilities/ap_shapes.jl")
    include("utilities/assign_flux.jl")
    include("utilities/galaxy_centre.jl")
    include("utilities/galaxy_orient.jl")
    include("utilities/particle_kinematics.jl")
    include("utilities/sim_to_galaxy.jl")
end
