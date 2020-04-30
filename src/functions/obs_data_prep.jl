# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    obs_data_prep(galaxy_data, ifu, envir)

This function prepares the particle data for a given observation with the given telescope.

Parameters:\n
    galaxy_data         Array of `Particle` describing galaxy
    ifu                 Struct of type `Telescope`
    envir               Struct of type `Environment`

Returns:\n
    galaxy_data         Array of particle data formatted for the observation specified by the parameters.
    parts_in_cell       3D array of the particles corresponding to each element in the IFU data-cube.
    observe             Struct of type `Observation` containing all observation parameters.
"""
function obs_data_prep(galaxy_data::Array{Galaxy_particle, 1},
                        ifu::Telescope,
                        envir::Environment)

    ap_size     = envir.ang_size * ifu.fov                        # diameter size of the telescope, kpc
    sbinsize    = ap_size / ifu.sbin                        # spatial bin size (kpc per bin)

    set_observables!.(galaxy_data, envir.inc_deg)     #set each particles mutable struct `Observables` for the given observation inclination

    deleteat!(galaxy_data, findall(part -> typeof(part) == Galaxy_dark, galaxy_data)) #remove all non-luminous/dark particles
    if(length(galaxy_data) == 0)
        error("There are no particles representing luminous matter in this simulation (i.e. no stars, bulge or disc particles).")
    end

    deleteat!(galaxy_data, findall(part -> part.r >= envir.r200, galaxy_data))  # remove particles beyond r200

    if (ifu.ap_shape == "circular")                    # remove particles outside aperture
      galaxy_data  = circular_ap_cut(galaxy_data, ap_size)

    elseif (ifu.ap_shape == "square")
        galaxy_data = square_ap_cut(galaxy_data, ifu.sbin, sbinsize)

    elseif (ifu.ap_shape == "hexagonal")
        galaxy_data = hexagonal_ap_cut(galaxy_data, ifu.sbin, sbinsize)
    end

    x = getfield.(galaxy_data, :x)
    z_obs = getfield.(getfield.(galaxy_data, :obs), :z_obs)
    vy_obs = getfield.(getfield.(galaxy_data, :obs), :vy_obs)

    max_vy_obs = maximum(vy_obs)
    min_vy_obs = minimum(vy_obs)

    vbin::Int64 = ceil((max_vy_obs - min_vy_obs) / ifu.vbinsize)    # number of velocity bins
    vbin_range = vbin * ifu.vbinsize / 2
    vseq = -vbin_range : ifu.vbinsize : vbin_range

    sbin_range = ifu.sbin * sbinsize / 2
    sseq = -sbin_range : sbinsize : sbin_range

    x_coord = searchsortedlast.(Ref(sseq), x)
    y_coord = searchsortedlast.(Ref(sseq), z_obs)
    z_coord = searchsortedlast.(Ref(vseq), vy_obs)

    x_invalid = findall(x-> x == 0 || x == length(sseq), x_coord)   #find indexs for particles outside of seen range
    y_invalid = findall(y->y==0 || y == length(sseq), y_coord)
    z_invalid = findall(z->z==0 || z == length(vseq), z_coord)
    invalid = sort(unique(vcat(x_invalid, y_invalid, z_invalid)))

    bins = x_coord + ifu.sbin * y_coord + ifu.sbin^2 * z_coord .- (ifu.sbin^2 + ifu.sbin)

    deleteat!(bins, invalid)
    deleteat!(galaxy_data, invalid)

    used_cells = unique(bins)
    parts_in_cell = [Galaxy_particle[] for i = 1:(ifu.sbin^2)*vbin]

    [parts_in_cell[cell] = galaxy_data[findall(x -> x == cell, bins)] for cell in used_cells]

    observe = Observation(envir.z,
                            envir.inc_deg,
                            envir.r200,
                            envir.blur,
                            ifu.ap_region,
                            ifu.sbin,
                            sseq,
                            vbin,
                            vseq,
                            ifu.lsf_size,
                            envir.ang_size,
                            sbinsize,
                            envir.mass2light,
                            envir.lum_dist,
                            envir.redshift_coef)

    return  galaxy_data, parts_in_cell, observe
end

function obs_data_prep(sim_data::Array{Sim_particle, 1},
                        ifu::Telescope,
                        envir::Environment)

    galaxy_data = sim_to_galaxy(sim_data)

    return obs_data_prep(galaxy_data, ifu, envir)
end
