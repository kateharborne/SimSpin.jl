# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    sim_analysis(sim_data;
                    bin_type="r",
                    rmax=200,
                    rbin=200,
                    dm_profile=nothing)

The purpose of this function is to calculate the kinematic properties of a simulated galaxy,
specifically the velocity and dispersion of particles within certain user defined bins. The user
must specify which direction they wish to study the kinematics using `bin_type` (where
`r` specifies out radial 3D spherical bins, `cr` specifies radial 2D circular
bins, and `z` specifies directly in 1D out of the plane of the galaxy). The default,
`bin_type = "r"` will return a comprehensive list of the simulation's kinematic properties
(i.e. contained radial mass and densities, velocities, velocity dispersions, anisotropy,
rotational and circular velocities and spin parameter). Other options for `bin_type` will
not contain the anisotropy, rotational and circular velocities or spin parameter as these are
only physical when defined using 3D spherical shells.

    Keyword arguments (optional):\n
        bin_type    The direction in which to bin the simulation model - `r` (default) bins
radially in 3D spherical shells, `cr` bins radially in 2D circular rings, `z` bins
in 1D off the plane of the galaxy.
        rmax        The maximum radial coordinate considered within the simulated galaxy in kpc.
        rbin        The number of radial bins considered.
        dm_profile  If dark matter particles are not included in the analysis, this option allows
you to use the DM profile for the mass distribution such that the circular velocity can be
correctly determined. Options include nothing (default),
`Dict("profile"=>"NFW", "dm_vm"=>186.9, "dm_a"=>34.5, "dm_rhof"=>0.035)` - where dm_vm is the virial mass, dm_a is the
scale radius, and dm_rhof is the density evaluated at the flattening radius - and
`Dict("profile"=>"Hernquist", "dm_mass"=>184.9, "dm_a"=>34.5)` - where dm_mass is the total mass of the dark
matter component and dm_a is the scale radius of the halo.
"""
function sim_analysis(galaxy_data::Array{Galaxy_particle, 1};
                      bin_type::String="r",
                      rmax::Int64=200,
                      rbin::Int64=200,
                      dm_profile::Union{Dark_matter, Nothing} = nothing)

    if(!any(getfield.(galaxy_data, :type) == "Dark Matter") && isnothing(dm_profile))
        error("No dark matter component is included in this analysis. Describe an analytic potential to calculate the total kinematic profile correctly.")
    end

    if bin_type == "r"
        if !isnothing(dm_profile)
            return r_shell(galaxy_data, rmax, rbin, dm_profile=dm_profile)
        else
            return r_shell(galaxy_data, rmax, rbin)
        end
    elseif bin_type == "cr"
        if !isnothing(dm_profile)
            return cr_shell(galaxy_data, rmax, rbin, dm_profile=dm_profile)
        else
            return cr_shell(galaxy_data, rmax, rbin)
        end
    elseif bin_type == "z"
        if !isnothing(dm_profile)
            return z_shell(galaxy_data, rmax, rbin, dm_profile=dm_profile)
        else
            return z_shell(galaxy_data, rmax, rbin)
        end
    else
        error("The specified bin type: ", bin_type, "is not supported. Please use 'r', 'cr' or 'z'.")
    end
end

function sim_analysis(sim_data::Array{Sim_particle, 1};
                      bin_type::String="r",
                      rmax::Int64=200,
                      rbin::Int64=200,
                      dm_profile::Union{Dark_matter, Nothing} = nothing)

    galaxy_data = sim_to_galaxy(sim_data)
    sim_data = nothing

    return sim_analysis(galaxy_data,
                        bin_type = bin_type,
                        rmax = rmax,
                        rbin = rbin,
                        dm_profile = dm_profile)
end
