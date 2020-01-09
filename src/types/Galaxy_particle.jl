# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

mutable struct Observable
    z_obs::Float64
    r_obs::Float64
    vy_obs::Float64
end

struct Galaxy_dark <: Galaxy_particle
    id::Int64
    type::String
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    r::Float64
    cr::Float64
    theta::Float64
    phi::Float64
    vr::Float64
    vt::Float64
    vcr::Float64
    vtheta::Float64
    vphi::Float64
    mass::Float64
    jx::Float64
    jy::Float64
    jz::Float64
    obs::Observable

        function Galaxy_dark(id::Int64, type::String,
                            x::Number, y::Number, z::Number,
                            vx::Number, vy::Number, vz::Number,
                            r::Number, cr::Number, theta::Number, phi::Number,
                            vr::Number, vt::Number, vcr::Number,
                            vtheta::Number, vphi::Number,
                            mass::Number,
                            jx::Number, jy::Number, jz::Number)
            obs = Observable(0., 0., 0.)

            new(id, type,
                x, y, z,
                vx, vy, vz,
                r, cr, theta, phi,
                vr, vt, vcr,
                vtheta, vphi,
                mass, jx, jy, jz,
                obs)
        end
end

struct Galaxy_lum <: Galaxy_particle
    id::Int64
    type::String
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    r::Float64
    cr::Float64
    theta::Float64
    phi::Float64
    vr::Float64
    vt::Float64
    vcr::Float64
    vtheta::Float64
    vphi::Float64
    mass::Float64
    jx::Float64
    jy::Float64
    jz::Float64
    obs::Observable

    function Galaxy_lum(id::Int64, type::String,
                        x::Number, y::Number, z::Number,
                        vx::Number, vy::Number, vz::Number,
                        r::Number, cr::Number, theta::Number, phi::Number,
                        vr::Number, vt::Number, vcr::Number,
                        vtheta::Number, vphi::Number,
                        mass::Number,
                        jx::Number, jy::Number, jz::Number)
        obs = Observable(0., 0., 0.)

        new(id, type,
            x, y, z,
            vx, vy, vz,
            r, cr, theta, phi,
            vr, vt, vcr,
            vtheta, vphi,
            mass, jx, jy, jz,
            obs)
    end
end

struct Galaxy_ssp <: Galaxy_particle
    id::Int64
    type::String
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    r::Float64
    cr::Float64
    theta::Float64
    phi::Float64
    vr::Float64
    vt::Float64
    vcr::Float64
    vtheta::Float64
    vphi::Float64
    mass::Float64
    jx::Float64
    jy::Float64
    jz::Float64
    age::Float64
    metallicity::Float64
    initial_mass::Float64
    obs::Observable

    function Galaxy_ssp(id::Int64, type::String,
                        x::Number, y::Number, z::Number,
                        vx::Number, vy::Number, vz::Number,
                        r::Number, cr::Number, theta::Number, phi::Number,
                        vr::Number, vt::Number, vcr::Number,
                        vtheta::Number, vphi::Number,
                        mass::Number,
                        jx::Number, jy::Number, jz::Number,
                        age::Number, metallicity::Number, initial_mass::Number)

        obs = Observable(0., 0., 0.)

        new(id, type,
            x, y, z,
            vx, vy, vz,
            r, cr, theta, phi,
            vr, vt, vcr,
            vtheta, vphi,
            mass, jx, jy, jz,
            age, metallicity, initial_mass, obs)
    end
end

"""
function galaxy_particle(particle, centre, rot_mat)

    Converts a given Sim_particle with values from the simulation reference frame
    to a Galaxy_particle in the galaxy reference frame. Some kinematic properties also calculated.
"""
function galaxy_particle(particle::Sim_particle, centre::Array{Float64, 1}, rot_mat::Array{Float64, 2})
    id = particle.id
    type = particle.type

    space_coords =  rot_mat * [(particle.x - centre[1]), (particle.y - centre[2]), (particle.z - centre[3])]
    x = space_coords[1]
    y = space_coords[2]
    z = space_coords[3]

    vel_coords = rot_mat * [(particle.vx - centre[4]), (particle.vy - centre[5]), (particle.vz - centre[6])]
    vx = vel_coords[1]
    vy = vel_coords[2]
    vz = vel_coords[3]

    mass = particle.mass

    r = get_r(x, y, z)
    cr = get_cr(x, y)
    theta, phi = get_spherical(x, y, z, r)
    vr = get_vr(x, y, z, vx, vy, vz)
    vt = get_vt(x, y, z, vx, vy, vz, r)
    vcr = get_vcr(x, y, vx, vy, cr)
    vtheta = get_vtheta(z, vz, r, vr)
    vphi = get_vphi(x, y, vx, vy, r, theta)
    jx, jy, jz = get_j(x, y, z, vx, vy, vz, mass)

    if(typeof(particle) == Sim_dark)
        return Galaxy_dark(id, type,
                            x, y, z,
                            vx, vy, vz,
                            r, cr, theta, phi,
                            vr, vt, vcr,
                            vtheta, vphi,
                            mass, jx, jy, jz)

    elseif(typeof(particle) == Sim_lum)
        return Galaxy_lum(id, type,
                            x, y, z,
                            vx, vy, vz,
                            r, cr, theta, phi,
                            vr, vt, vcr,
                            vtheta, vphi,
                            mass, jx, jy, jz)

    elseif(typeof(particle) == Sim_ssp)
        age = particle.age
        metallicity=particle.metallicity
        initial_mass=particle.initial_mass
        return Galaxy_ssp(id, type,
                            x, y, z,
                            vx, vy, vz,
                            r, cr, theta, phi,
                            vr, vt, vcr,
                            vtheta, vphi,
                            mass, jx, jy, jz,
                            age, metallicity, initial_mass)
    end
end

function set_observables!(particle::Galaxy_particle, inc_deg::Int64)
    particle.obs.z_obs = sind(inc_deg) * particle.z + cosd(inc_deg) * particle.y
    particle.obs.vy_obs = cosd(inc_deg) * particle.vz - sind(inc_deg) * particle.vy
    particle.obs.r_obs = sqrt(particle.x^2 + particle.obs.z_obs^2)
end
