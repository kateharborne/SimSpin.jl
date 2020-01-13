# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    obs_data_prep(galaxy_data;
                    r200 = 200,
                    z = 0.05,
                    fov = 15.,
                    ap_shape = "circular",
                    central_wvl = 4800,
                    lsf_fwhm = 2.65,
                    pixel_sscale = 0.5,
                    pixel_vscale = 1.04,
                    inc_deg = 70)

This function prepares the particle data for a given observation.

Returns:\n
    galaxy_data         Array of particle data formatted for the observation specified by the parameters.
    parts_in_cell       3D array of the particles corresponding to each element in the IFU data-cube.
    ap_region           The aperture region mask used to remove flux outside of the specified aperture.
    sbin                The number of spatial bins in the aperture.
    vbin                The number of velocity bins in the flux grid.
    vseq                The bounds of each velocity bin in the flux grid.
    lsf_size            The Gaussian standard deviation of the line spread function in km/s.
    ang_size            The angular size given redshift z in kpc
    sbinsize            The spatial bin size in kpc per bin

Keyword arguments (optional):\n
    r200            The virial radius specified in the simulation, kpc.
    z               The projected redshift at which the observation is made.
    fov             The field of view of the IFU, diameter in arcseconds.
    ap_shape        The shape of the field of view, with options "circular", "square" or "hexagonal".
    central_wvl     The central filter wavelength used for the observation, given in angstroms.
    lsf_fwhm        The line spread function full-width half-max, given in angstroms.
    pixel_sscale    The corresponding spatial pixel scale associated with a given telescope output in arcseconds.
    pixel_vscale    The corresponding velocity pixel scale associated with a given telescope filter output in angstroms.
    inc_deg         The inclination at which to observe the galaxy in degrees.
"""
function obs_data_prep(galaxy_data::Array{Galaxy_particle, 1};
                        r200::Int64 = 200,
                        z::Float64 = 0.05,
                        fov::Float64 = 15.,
                        ap_shape::String = "circular",
                        central_wvl::Int64=4800,
                        lsf_fwhm::Float64=2.65,
                        pixel_sscale::Float64=0.5,
                        pixel_vscale::Float64=1.04,
                        inc_deg::Int64=70)

    ang_size    = angleSize(z)                      # angular size given z, kpc
    ap_size     = ang_size * fov                    # diameter size of the telescope, kpc
    sbin::Int64 = floor(fov / pixel_sscale)         # number of spatial bins
    sbinsize    = ap_size / sbin                    # spatial bin size (kpc per bin)
    vbinsize    = (pixel_vscale / central_wvl) * (3e8 / 1e3)  # km/s per velocity bin
    lsf_size    = ((lsf_fwhm / central_wvl) * (3e8 / 1e3)) / (2 * sqrt(2*log(2))) # velocity uncertainty (sd)

    set_observables!.(galaxy_data, inc_deg)

    if (ap_shape == "circular")       # circular apperture mask
      ap_region = circular_ap(sbin)
    elseif (ap_shape == "square")     # square apperture mask
      ap_region = square_ap(sbin)
    elseif (ap_shape == "hexagonal")  # hexagonal apperture mask
      ap_region = hexagonal_ap(sbin)
    else
       error("Unsupported aperture shape specified.")
    end

    deleteat!(galaxy_data, findall(part -> typeof(part) == Galaxy_dark, galaxy_data)) #remove all non-luminous/dark particles
    if(length(galaxy_data) == 0)
        error("There are no particles representing luminous matter in this simulation (i.e. no stars, bulge or disc particles).")
    end

    deleteat!(galaxy_data, findall(part -> part.r >= r200, galaxy_data))  # remove particles beyond r200

    if (ap_shape == "circular")                    # remove particles outside aperture
      galaxy_data  = circular_ap_cut(galaxy_data, ap_size)

    elseif (ap_shape == "square")
      galaxy_data = square_ap_cut(galaxy_data, sbin, sbinsize)

    elseif (ap_shape == "hexagonal")
      galaxy_data = hexagonal_ap_cut(galaxy_data, sbin, sbinsize)
    end

    x = getfield.(galaxy_data, :x)
    z_obs = getfield.(getfield.(galaxy_data, :obs), :z_obs)
    vy_obs = getfield.(getfield.(galaxy_data, :obs), :vy_obs)

    max_vy_obs = maximum(vy_obs)
    min_vy_obs = minimum(vy_obs)

    vbin::Int64 = ceil((max_vy_obs - min_vy_obs) / vbinsize)    # number of velocity bins
    vbin_range = vbin * vbinsize / 2
    vseq = -vbin_range : vbinsize : vbin_range

    sbin_range = sbin * sbinsize / 2
    sseq = -sbin_range : sbinsize : sbin_range

    bins = searchsortedlast.(Ref(sseq), x) + sbin .* searchsortedlast.(Ref(sseq), z_obs) +
            sbin^2 .* searchsortedlast.(Ref(vseq), vy_obs) .- (sbin^2 + sbin)

    parts_in_cell = [Galaxy_particle[] for i = 1:sbin*sbin*vbin]

    for cell in 1:length(parts_in_cell)
        this = galaxy_data[findall(x->x==cell, bins)]
        parts_in_cell[cell] = this
    end

    return  galaxy_data,
            parts_in_cell,
            ap_region,
            sbin,
            vbin,
            vseq,
            lsf_size,
            ang_size,
            sbinsize
end

function obs_data_prep(sim_data::Array{Sim_particle, 1};
                        r200::Int64 = 200,
                        z::Float64 = 0.05,
                        fov::Float64 = 15.,
                        ap_shape::String = "circular",
                        central_wvl::Int64=4800,
                        lsf_fwhm::Float64=2.65,
                        pixel_sscale::Float64=0.5,
                        pixel_vscale::Float64=1.04,
                        inc_deg::Int64=70)

    galaxy_data = sim_to_galaxy(sim_data)

    return obs_data_prep(galaxy_data, r200=r200, z=z, fov=fov,
                            ap_shape=ap_shape, central_wvl=central_wvl,
                            lsf_fwhm=lsf_fwhm,
                            pixel_sscale=pixel_sscale, pixel_vscale=pixel_vscale,
                            inc_deg=inc_deg)
end
