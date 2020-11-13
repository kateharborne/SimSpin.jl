# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using HDF5

"""
    sim_data(filename;
            pytpe = [],
            ssp = false)

Reads in a [SimSpin format](https://github.com/kateharborne/create_SimSpinFile) HDF5 file at location, `filename`.
Returns array of Sim_particles.

Keyword arguments (optional):\n
    ptype       A vector of the particles types to be read in e.g. ptype = [1,3].
                If omitted all particles types will be read.
    ssp         Boolean value to use ssp particle information.
"""
function sim_data(filename::String;
                    ptype::Vector{} = [],
                    ssp::Bool = false)

    sim_data = Sim_particle[]
    count = 1

    translation = Dict(
                "PartType0" => "Gas",
                "PartType1" => "Dark Matter",
                "PartType2" => "Disc",
                "PartType3" => "Bulge",
                "PartType4" => "Star")

    h5open(filename, "r") do galaxy_file
        ppart = keys(galaxy_file)

        if(isempty(ptype)) ptype = ppart
        else ptype = string.("PartType", ptype)
        end

        has_no_data = setdiff(ptype, ppart)
        if(!isempty(has_no_data))
            error("The requested particle type, ", has_no_data, ", are not present in this model.")
        end

        for type in ptype

            x_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/x")])
            y_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/y")])
            z_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/z")])
            vx_arr::Array{Float64, 1}= read(galaxy_file[string(type, "/vx")])
            vy_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/vy")])
            vz_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/vz")])
            mass_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/Mass")])

            if(type == "PartType0" || type == "PartType1")
                for i = 1:length(x_arr)
                    particle = Sim_dark(
                        count, translation[type],
                        x_arr[i], y_arr[i], z_arr[i],
                        vx_arr[i], vy_arr[i], vz_arr[i],
                        mass_arr[i])
                    count += 1
                    push!(sim_data, particle)
                end
            else()
                if(!ssp)
                    for i = 1:length(x_arr)
                        particle = Sim_lum(
                            count, translation[type],
                            x_arr[i], y_arr[i], z_arr[i],
                            vx_arr[i], vy_arr[i], vz_arr[i],
                            mass_arr[i])
                        count += 1
                        push!(sim_data, particle)
                    end

                    if type == "PartType4" && length(keys(galaxy_file[type])) > 7
                        @warn("SSP data is available for Star particles in this simulation file but has not been read in. If spectra is desired set ssp=true.")
                    end

                elseif(ssp && length(keys(galaxy_file[type])) > 7)
                    if any(isequal.(keys(galaxy_file[type]), "Age"))
                        age_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/Age")])
                    elseif any(isequal.(keys(galaxy_file[type]), "StellarFormationTime"))
                        sft::Array{Float64, 1} = read(galaxy_file[string(type, "/StellarFormationTime")])
                        z_sft = ((1 ./ sft) .- 1)
                        age_arr = celestial.cosdistTravelTime.(z_sft)     #convert stellar formation time to particle age
                    else
                        error("SSP particle data must include either particle ages or particle stellar formation time.")
                    end

                    met_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/Metallicity")])

                    if any(isequal.(keys(galaxy_file[type]), "InitialMass"))
                        init_mass_arr::Array{Float64, 1} = read(galaxy_file[string(type, "/InitialMass")])
                    else
                        init_mass_arr = zeros(Float64, length(x_arr))
                    end

                    for i = 1:length(x_arr)
                        particle = Sim_ssp(
                            count, translation[type],
                            x_arr[i], y_arr[i], z_arr[i],
                            vx_arr[i], vy_arr[i], vz_arr[i],
                            mass_arr[i],
                            age_arr[i],
                            met_arr[i],
                            init_mass_arr[i])
                        count += 1
                        push!(sim_data, particle)
                    end

                elseif(ssp && length(keys(galaxy_file[type])) <= 7)
                    error("SSP requested, but no Age/Metallicity information contained within supplied file, ", filename, ". \n",
                            "Please set SSP=false, or provide additional particle information.")

                end
            end
        end
    end
    return sim_data
end

"""
    sim_data(;ptype::Vector{} = [],
                ssp::Bool = false)

if name filename is specified `sim_data` returns particle data from SimSpin's example file.

Keyword arguments (optional):\n
    ptype       A vector of the particles types to be read in e.g. ptype = [1,3].
                If omitted all particles types will be read.
    ssp         Boolean value to use ssp particle information.
"""
function sim_data(;ptype::Vector{} = [],
                    ssp::Bool = false)
    filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_example.hdf5")

    return sim_data(filename, ptype=ptype, ssp=ssp)
end
