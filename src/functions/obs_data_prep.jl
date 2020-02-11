# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    obs_data_prep(galaxy_data, observation, ifu)

This function prepares the particle data for a given observation with the given telescope.

Returns:\n
    galaxy_data         Array of particle data formatted for the observation specified by the parameters.
    parts_in_cell       3D array of the particles corresponding to each element in the IFU data-cube.
    observed            Struct of type `Observation` containing all observation parameters.

Parameters:\n
    galaxy_data         Array of `Particle` describing galaxy
    ifu                 Struct of type `Telescope`
    z                   The projected redshift at which the observation is made.
    inc_deg             The inclination at which to observe the galaxy in degrees.
    r200                The virial radius specified in the simulation, kpc.
"""
function obs_data_prep(galaxy_data::Array{Galaxy_particle, 1},
                        ifu::Telescope,
                        z::Float64,
                        inc_deg::Real,
                        r200::Real)

    ang_size    = angleSize(z)                      # angular size given z, kpc
    ap_size     = ang_size * ifu.fov                # diameter size of the telescope, kpc
    sbinsize    = ap_size / ifu.sbin                # spatial bin size (kpc per bin)

    set_observables!.(galaxy_data, inc_deg)     #set each particles mutable struct `Observables` for the given observation inclination

    deleteat!(galaxy_data, findall(part -> typeof(part) == Galaxy_dark, galaxy_data)) #remove all non-luminous/dark particles
    if(length(galaxy_data) == 0)
        error("There are no particles representing luminous matter in this simulation (i.e. no stars, bulge or disc particles).")
    end

    deleteat!(galaxy_data, findall(part -> part.r >= r200, galaxy_data))  # remove particles beyond r200

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

    bins = searchsortedlast.(Ref(sseq), x) + ifu.sbin .* searchsortedlast.(Ref(sseq), z_obs) +
            ifu.sbin^2 .* searchsortedlast.(Ref(vseq), vy_obs) .- (ifu.sbin^2 + ifu.sbin)

    parts_in_cell = [Galaxy_particle[] for i = 1:(ifu.sbin^2)*vbin]

    for cell in 1:length(parts_in_cell)
        this = galaxy_data[findall(x->x==cell, bins)]
        parts_in_cell[cell] = this
    end

    observe = Observation(z, inc_deg, r200, ifu.ap_region, ifu.sbin, vbin, vseq, ifu.lsf_size, ang_size, sbinsize)

    return  galaxy_data, parts_in_cell, observe
end

function obs_data_prep(sim_data::Array{Sim_particle, 1},
                        ifu::Telescope,
                        z::Float64,
                        inc_deg::Real,
                        r200::Real)

    galaxy_data = sim_to_galaxy(sim_data)

    return obs_data_prep(galaxy_data, ifu, z, inc_deg, r200)
end