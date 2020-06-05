# Date created: 12/02/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Environment(z, inc_deg, r200, mass2light, blur)

Creates a `struct` containing environmental parameters required for a mock observation of a simulated galaxy. Each argument can either be an array
or a single value. If an array is passed then a corresponding array of Environments will be returned containing every permutation of the given parameters.

Parameters:\n
    z               The projected redshift at which the observation is made.
    inc_deg         The inclination at which to observe the galaxy in degrees. Relative to face on, rotated around semi-major axis.
    r200            The virial radius specified in the simulation, kpc.
    mass2light      Optional. The mass to light ratio for non-ssp, luminous particles (disk and bulge). If not specified both particle types default to a mass to light ratio of 1. If a single value is specified both use that value and if a tuple is specified, bulge and disk mass-to-light values are taken respectively.
    blur            Optional. Struct of type `Blur` containing seeing information. If ommitted no blurring is used.
"""
struct Environment

    z::Float64
    inc_deg::Real
    r200::Real
    disc_mass2light::Real
    bulge_mass2light::Real
    blur::Union{Blur, Nothing}
    lum_dist::Real              #luminosity distance in Mpc
    ang_size::Real              # angular size in kpc
    redshift_coef::Real

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Tuple{Real, Real},
                        blur::Blur)

        lum_dist = celestial.cosdistLumDist(z, ref="Planck")
        ang_size = celestial.cosdistAngScale(z, ref="Planck")
        redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)

        new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Tuple{Real, Real},
                        blur::Nothing)

        lum_dist = celestial.cosdistLumDist(z, ref="Planck")
        ang_size = celestial.cosdistAngScale(z, ref="Planck")
        redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)

        new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Real,
                        blur::Blur)

        mass2light = (mass2light, mass2light)
        Environment(z, inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Real,
                        blur::Nothing)

        mass2light = (mass2light, mass2light)
        Environment(z, inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        blur::Blur)
        mass2light = (1., 1.)
        Environment(z, inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Tuple{Real, Real})
        blur = nothing
        Environment(z,inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Real)

        mass2light = (mass2light, mass2light)
        blur = nothing
        Environment(z,inc_deg, r200, mass2light, blur)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real)
        mass2light = (1., 1.)
        blur = nothing
        Environment(z, inc_deg, r200, mass2light, blur)
    end

    function Environment(z_array::Union{Float64, Array{Float64, 1}, AbstractRange{Float64}},
                        inc_deg_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        r200_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        mass2light_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        blur_array::Union{Blur, Array{Blur, 1}, AbstractRange{<:Real}})

        envir_array = Environment[]

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    for mass2light in mass2light_array
                        for blur in blur_array
                            envir = new(z, inc_deg, r200, mass2light, mass2light, blur, lum_dist, ang_size, redshift_coef)
                            push!(envir_array, envir)
                        end
                    end
                end
            end
        end

        return envir_array
    end

    function Environment(z_array::Union{Float64, Array{Float64, 1}, AbstractRange{Float64}},
                        inc_deg_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        r200_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        mass2light_array::Union{Tuple{Real, Real}, Array{<:Tuple, 1}},
                        blur_array::Union{Blur, Array{Blur, 1}, AbstractRange{<:Real}})

        envir_array = Environment[]
        m2l_is_array = isa(mass2light_array, Array{<:Tuple, 1})

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    for blur in blur_array
                        if m2l_is_array
                            for mass2light in mass2light_array
                                envir = new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
                                push!(envir_array, envir)
                            end
                        else
                            mass2light = mass2light_array
                            envir = new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
                            push!(envir_array, envir)
                        end
                    end
                end
            end
        end

        return envir_array
    end

    function Environment(z_array::Union{Float64, Array{Float64, 1}, AbstractRange{Float64}},
                        inc_deg_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        r200_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        blur_array::Union{Blur, Array{Blur, 1}})

        envir_array = Environment[]

        mass2light = (1., 1.)

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    for blur in blur_array
                        envir = new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
                        push!(envir_array, envir)
                    end
                end
            end
        end

        return envir_array
    end

    function Environment(z_array::Union{Float64, Array{Float64, 1}, AbstractRange{Float64}},
                        inc_deg_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        r200_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        mass2light_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}})

        envir_array = Environment[]

        blur = nothing

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    for mass2light in mass2light_array
                        envir = new(z, inc_deg, r200, mass2light, mass2light, blur, lum_dist, ang_size, redshift_coef)
                        push!(envir_array, envir)
                    end
                end
            end
        end

        return envir_array
    end

    function Environment(z_array::Union{Float64, Array{Float64, 1}, AbstractRange{Float64}},
                        inc_deg_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        r200_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        mass2light_array::Union{Tuple{Real, Real}, Array{<:Tuple, 1}})

        envir_array = Environment[]

        blur = nothing
        m2l_is_array = isa(mass2light_array, Array{<:Tuple, 1})

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    if m2l_is_array
                        for mass2light in mass2light_array
                            envir = new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
                            push!(envir_array, envir)
                        end
                    else
                        mass2light = mass2light_array
                        envir = new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
                        push!(envir_array, envir)
                    end
                end
            end
        end

        return envir_array
    end

    function Environment(z_array::Union{Float64, Array{Float64, 1}, AbstractRange{Float64}},
                        inc_deg_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}},
                        r200_array::Union{Real, Array{<:Real, 1}, AbstractRange{<:Real}})

        envir_array = Environment[]

        mass2light = (1., 1.)
        blur = nothing

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    envir = new(z, inc_deg, r200, mass2light[1], mass2light[2], blur, lum_dist, ang_size, redshift_coef)
                    push!(envir_array, envir)
                end
            end
        end

        return envir_array
    end
end
