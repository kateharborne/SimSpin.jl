# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using HDF5

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
        ppart = names(galaxy_file)

        if(isempty(ptype)) ptype = ppart
        else ptype = string.("PartType", ptype)
        end

        has_no_data = setdiff(ptype, ppart)
        if(!isempty(has_no_data))
            error("The requested particle type, ", has_no_data, ", are not present in this model.")
        end

        for type in ptype

            x_arr = read(galaxy_file[string(type, "/x")])
            y_arr = read(galaxy_file[string(type, "/y")])
            z_arr = read(galaxy_file[string(type, "/z")])
            vx_arr = read(galaxy_file[string(type, "/vx")])
            vy_arr = read(galaxy_file[string(type, "/vy")])
            vz_arr = read(galaxy_file[string(type, "/vz")])
            mass_arr = read(galaxy_file[string(type, "/Mass")])

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

                elseif(ssp && length(names(galaxy_file[type])) == 10)
                    if(any(isequal.(names(galaxy_file[type]), "Age")))
                        age_arr = read(galaxy_file[string(type, "/Age")])
                    else
                        sft = read(galaxy_file[string(type, "/StellarFormationTime")])
                        age_arr = sftToAge(sft)
                    end

                    met_arr = read(galaxy_file[string(type, "/Metallicity")])
                    init_mass_arr = read(galaxy_file[string(type, "/InitialMass")])

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

                elseif(ssp && length(names(galaxy_file[type])) != 10)
                    error("SSP requested, but no Age/Metallicity information contained within supplied file, ", filename, ". \n",
                            "Please set SSP=false, or provide additional particle information.")

                end
            end
        end
    end

    return(sim_data)
end
