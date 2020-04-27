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
    mass2light      Optional. The mass to light ratio for non-ssp, luminous particles. Defaults to 1.
    blur            Optional. Struct of type `Blur` containing seeing information. If ommitted no blurring is used.
"""
struct Environment

    z::Float64
    inc_deg::Real
    r200::Real
    mass2light::Real
    blur::Union{Blur, Nothing}
    lum_dist::Real              #luminosity distance in Mpc
    ang_size::Real              # angular size in kpc
    redshift_coef::Real

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Real,
                        blur::Blur)

        lum_dist = celestial.cosdistLumDist(z, ref="Planck")
        ang_size = celestial.cosdistAngScale(z, ref="Planck")
        redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)

        new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        blur::Blur)
        mass2light = 1.
        lum_dist = celestial.cosdistLumDist(z, ref="Planck")
        ang_size = celestial.cosdistAngScale(z, ref="Planck")
        redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)

        new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real,
                        mass2light::Real)
        blur = nothing
        lum_dist = celestial.cosdistLumDist(z, ref="Planck")
        ang_size = celestial.cosdistAngScale(z, ref="Planck")
        redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)

        new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
    end

    function Environment(z::Float64,
                        inc_deg::Real,
                        r200::Real)
        mass2light = 1.
        blur = nothing
        lum_dist = celestial.cosdistLumDist(z, ref="Planck")
        ang_size = celestial.cosdistAngScale(z, ref="Planck")
        redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)

        new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
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
                            envir = new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
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

        mass2light = 1.

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    for blur in blur_array
                        envir = new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
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
                        envir = new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
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

        mass2light = 1.
        blur = nothing

        for z in z_array
            lum_dist = celestial.cosdistLumDist(z, ref="Planck")
            ang_size = celestial.cosdistAngScale(z, ref="Planck")
            redshift_coef = ProSpect.Lum2FluxFactor(z, lum_dist)
            for inc_deg in inc_deg_array
                for r200 in r200_array
                    envir = new(z, inc_deg, r200, mass2light, blur, lum_dist, ang_size, redshift_coef)
                    push!(envir_array, envir)
                end
            end
        end

        return envir_array
    end
end
