# Date created: 18/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    r_shell(galaxy_data, rmax, rbin; dm_profile)

3D spherical shell analysis
"""
function r_shell(galaxy_data::Array{Galaxy_particle, 1},
                        rmax::Int64,
                        rbin::Int64;
                        dm_profile::Union{Dark_matter, Nothing} = nothing)


    G = 4.516e-29   # gravitational constant in units of kpc^3/([1e10 Msolar]s^2)

    galaxy_data = galaxy_data[findall(x -> x.r < rmax, galaxy_data)] # remove particles further than rmax (spherical)

    z = getfield.(galaxy_data, :z)
    sigma_z = sqrt(mean(z.^2) - (mean(z)^2))

    r = getfield.(galaxy_data, :r)
    bin_bounds = 0:rmax/rbin:rmax
    part_bin = searchsortedlast.(Ref(bin_bounds), r)    # assigns each particle into an rbin

    parts_in_bin = [Galaxy_particle[] for i = 1:rbin]

    for bin in 1:rbin
        this = galaxy_data[findall(x->x==bin, part_bin)]
        parts_in_bin[bin] = this
    end

    if isnothing(dm_profile) #if no DM profile is specified
        dm_mass_profile = zeros(Float64, rbin)
    elseif typeof(dm_profile) == Hernquist
        dm_mass_profile = dm_profile.dm_mass .* bin_bounds[2:end] .^ 2 ./ (dm_profile.dm_a .+ bin_bounds[2:end]).^2
    elseif typeof(dm_profile) == NFW
        dm_mass_profile = 4 * pi * dm_profile.dm_rhof .* dm_profile.dm_a .^ 3 * (log.(1 .+ (bin_bounds[2:end] ./ dm_profile.dm_a)) - ((bin_bounds[2:end] ./ dm_profile.dm_a) ./ (1 .+ (bin_bounds[2:end] ./ dm_profile.dm_a))))
    else
        error("Dark matter analytic mass profile must be a dict containing 'profile' => 'Hernquist' or 'NFW'.")
    end

    profile = Spherical_3d[]

    for (index, bin) in enumerate(parts_in_bin)

        part_num = length(bin)              #number of particles in bin

        if part_num == 0
            mass = 0.0

            Jx = 0.0  # angular momentum components of shell
            Jy = 0.0
            Jz = 0.0

            each_vr = 0.0
            each_vtheta = 0.0
            each_vphi = 0.0

            each_vx = 0.0
            each_vz = 0.0
        else
            mass = sum(getfield.(bin, :mass))  # sums the mass of each particle in this radial bin

            Jx = sum(getfield.(bin, :jx))   # angular momentum components of shell
            Jy = sum(getfield.(bin, :jy))
            Jz = sum(getfield.(bin, :jz))

            each_vr = getfield.(bin, :vr)
            each_vtheta = getfield.(bin, :vtheta)
            each_vphi = getfield.(bin, :vphi)

            each_vx = getfield.(bin, :vx)
            each_vz = getfield.(bin, :vz)
        end

        rad = bin_bounds[index + 1]         # the outer edge of each radial bin

        mass_cum = mass + sum(getfield.(profile, :mass))    #cumulative mass from centre

        logp = log10(mass / ((4/3) * pi * (rad^3 - bin_bounds[index]^3)))   # log10 of shell density

        vc = pc_to_m * sqrt(G * (mass_cum + dm_mass_profile[index]) / rad) # the circular velocity of particles at this radius, km/s

        J = sqrt(Jx^2 + Jy^2 + Jz^2)    # magnitude of angular momentum

        Jx_cum = Jx + sum(getfield.(profile, :Jx))      # enclosed angular momentum components
        Jy_cum = Jy + sum(getfield.(profile, :Jy))
        Jz_cum = Jz + sum(getfield.(profile, :Jz))
        J_cum = sqrt(Jx_cum^2 + Jy_cum^2 + Jz_cum^2)    # enclosed magnitude of angular momentum

        mean_vr = mean(each_vr)             # mean radial velocity in the shell, km/s
        mean_vtheta = mean(each_vtheta)     # mean velocity along theta in the shell, km/s
        mean_vphi = mean(each_vphi)         # mean velocity along phi in the shell, km/s

        sigma_vr = sqrt(mean(each_vr.^2) - mean_vr^2)       # radial velocity dispersion, km/s
        sigma_vt = sqrt(mean(each_vtheta.^2) - mean_vtheta^2) + (mean(each_vphi.^2) - mean_vphi^2)  # tangential velocity dispersion, km/s

        mean_vx = mean(each_vx)
        sigma2_vx = sum((each_vx .- mean_vx).^2) / part_num

        mean_vz = mean(each_vz)
        sigma2_vz = sum((each_vz .- mean_vz).^2) / part_num

        b = 1 - ((sigma_vt^2) / (2 * sigma_vr^2))   # velocity anisotropy, unitless
        vrot = J / (mass * rad * pc_to_m)          # rotational velocity, km/s
        lamda = J_cum / (1.414214 * mass_cum * vc * rad * pc_to_m)   # spin parameter, unitless

        shell = Spherical_3d(rad, mass_cum, logp, vc,
                                Jx_cum, Jy_cum, Jz_cum, J_cum,
                                mean_vr, sigma_vr, sigma_vt,
                                sigma2_vx, sigma2_vz,
                                b, vrot, lamda)
        push!(profile, shell)
    end

    return profile, sigma_z
end


"""
    cr_shell(galaxy_data, rmax, rbin; dm_profile)

2D circular radial shell analysis
"""
function cr_shell(galaxy_data::Array{Galaxy_particle, 1},
                        rmax::Int64,
                        rbin::Int64;
                        dm_profile::Union{Dark_matter, Nothing} = nothing)

    G = 4.516e-29   # gravitational constant in units of kpc^3/([1e10 Msolar]s^2)

    galaxy_data = galaxy_data[findall(x -> (x.cr < rmax) && (abs(x.z) < rmax), galaxy_data)] # remove particles further than rmax (spherical)

    z = getfield.(galaxy_data, :z)
    sigma_z = sqrt(mean(z.^2) - mean(z)^2)

    r = getfield.(galaxy_data, :r)
    bin_bounds = 0:rmax/rbin:rmax
    part_bin = searchsortedlast.(Ref(bin_bounds), r)    # assigns each particle into an rbin

    parts_in_bin = [Galaxy_particle[] for i = 1:rbin]

    for bin in 1:rbin
        this = galaxy_data[findall(x->x==bin, part_bin)]
        parts_in_bin[bin] = this
    end

    profile = Circ_radial_2d[]

    for (index, bin) in enumerate(parts_in_bin)

        part_num = length(bin)              #number of particles in bin

        if part_num == 0
            mass = 0.0

            Jx = 0.0  # angular momentum components of shell
            Jy = 0.0
            Jz = 0.0

            each_vcr = 0.0
        else
            mass = sum(getfield.(bin, :mass))  # sums the mass of each particle in this radial bin

            Jx = sum(getfield.(bin, :jx))   # angular momentum components of shell
            Jy = sum(getfield.(bin, :jy))
            Jz = sum(getfield.(bin, :jz))

            each_vcr = getfield.(bin, :vcr)
        end


        cr = bin_bounds[index + 1]          # the outer edge of each radial bin

        mass_cum = mass + sum(getfield.(profile, :mass))    #cumulative mass from centre

        logp = log10(mass / (pi * (rmax * 2) * (cr^2 - bin_bounds[index]^2)))   # log10 of shell density

        Jx_cum = Jx + sum(getfield.(profile, :Jx))      # enclosed angular momentum components
        Jy_cum = Jy + sum(getfield.(profile, :Jy))
        Jz_cum = Jz + sum(getfield.(profile, :Jz))
        J_cum = sqrt(Jx_cum^2 + Jy_cum^2 + Jz_cum^2)    # enclosed magnitude of angular momentum

        mean_vcr = mean(each_vcr)
        sigma_vcr = sqrt(mean(each_vcr.^2) - mean_vcr^2)

        shell = Circ_radial_2d(cr, mass_cum, logp,
                                Jx_cum, Jy_cum, Jz_cum, J_cum,
                                mean_vcr, sigma_vcr)
        push!(profile, shell)
    end

    return profile, sigma_z
end

"""
    z_shell(galaxy_data, rmax, rbin; dm_profile)

2D circular planar shell analysis
"""
function z_shell(galaxy_data::Array{Galaxy_particle, 1},
                        rmax::Int64,
                        rbin::Int64;
                        dm_profile::Union{Dark_matter, Nothing} = nothing)

    G = 4.516e-29   # gravitational constant in units of kpc^3/([1e10 Msolar]s^2)

    galaxy_data = galaxy_data[findall(x -> (x.cr < rmax) && (x.z < rmax) && (x.z > 0), galaxy_data)] # remove particles further than rmax (spherical)

    z = getfield.(galaxy_data, :z)
    sigma_z = sqrt(mean(z.^2) - mean(z)^2)

    bin_bounds = 0:rmax/rbin:rmax
    part_bin = searchsortedlast.(Ref(bin_bounds), z)    # assigns each particle an rbin

    parts_in_bin = [Galaxy_particle[] for i = 1:rbin]

    for bin in 1:rbin
        this = galaxy_data[findall(x->x==bin, part_bin)]    #Assigns all rbins a list of particles
        parts_in_bin[bin] = this
    end

    profile = Circ_planar_2d[]

    for (index, bin) in enumerate(parts_in_bin)

        part_num = length(bin)              #number of particles in bin

        if part_num == 0
            mass = 0.0

            Jx = 0.0  # angular momentum components of shell
            Jy = 0.0
            Jz = 0.0

            each_vz = 0.0
        else
            mass = sum(getfield.(bin, :mass))  # sums the mass of each particle in this radial bin

            Jx = sum(getfield.(bin, :jx))   # angular momentum components of shell
            Jy = sum(getfield.(bin, :jy))
            Jz = sum(getfield.(bin, :jz))

            each_vz = getfield.(bin, :vz)
        end

        zbin = bin_bounds[index + 1]           # the outer edge of each radial bin

        mass_cum = mass + sum(getfield.(profile, :mass))    #cumulative mass from centre

        logp = log10(mass / (pi * (rmax^2) * (zbin - bin_bounds[index])))   # log10 of shell density

        Jx_cum = Jx + sum(getfield.(profile, :Jx))      # enclosed angular momentum components
        Jy_cum = Jy + sum(getfield.(profile, :Jy))
        Jz_cum = Jz + sum(getfield.(profile, :Jz))
        J_cum = sqrt(Jx_cum^2 + Jy_cum^2 + Jz_cum^2)    # enclosed magnitude of angular momentum

        mean_vz = mean(each_vz)
        sigma_vz = sqrt(mean(each_vz.^2) - mean_vz^2)

        shell = Circ_planar_2d(zbin, mass_cum, logp,
                                Jx_cum, Jy_cum, Jz_cum, J_cum,
                                mean_vz, sigma_vz)
        push!(profile, shell)
    end

    return profile, sigma_z
end
