# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

function sim_analysis(sim_data::Array{Galaxy_particle, 1};
                      bin_type::String="r",
                      rmax::Int64=200,
                      rbin::Int64=200,
                      dm_profile=nothing)

    if(!any(getfield.(sim_data, :type) == "Dark Matter") && isnothing(dm_profile))
        error("No dark matter component is included in this analysis. Describe an analytic potential to calculate the total kinematic profile correctly.")
    end
end

function sim_analysis(sim_data::Array{Sim_particle, 1};
                      bin_type::String="r",
                      rmax::Int64=200,
                      rbin::Int64=200,
                      dm_profile=nothing)

    galaxy_data = sim_to_galaxy(sim_data)
    sim_data = nothing

    return sim_analysis(galaxy_data,
                        bin_type = bin_type,
                        rmax = rmax,
                        rbin = rbin,
                        dm_profile = dm_profile)
end
