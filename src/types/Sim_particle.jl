# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Sim_dark

A Sim_particle type with no luminosity.
Can be gas (PartType0) or dark matter (PartType1).
"""
struct Sim_dark <: Sim_particle
    id::Int64
    type::String
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    mass::Float64

    function Sim_dark(id::Int64, type::String,
        x::Number, y::Number, z::Number,
        vx::Number, vy::Number, vz::Number,
        mass::Number)

        new(id, type,
        x, y, z,
        vx, vy, vz,
        mass)
    end
end

"""
    Sim_lum

    A Sim_particle type with luminosity but no associated spectra.
    Can be of type bulge (PartType2), disc (PartType3) or star (PartType4).
"""
struct Sim_lum <: Sim_particle #possibly not necessary lum assigned at cell stage
    id::Int64
    type::String
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    mass::Float64

    function Sim_lum(id::Int64, type::String,
        x::Number, y::Number, z::Number,
        vx::Number, vy::Number, vz::Number,
        mass::Number)

        new(id, type,
        x, y, z,
        vx, vy, vz,
        mass)
    end
end

"""
    Sim_ssp

    A Sim_particle type with luminosity and associated spectra.
    Can be of type bulge (PartType2), disc (PartType3) or star (PartType4).
"""
struct Sim_ssp <: Sim_particle
    id::Int64
    type::String
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    mass::Float64
    age::Float64
    metallicity::Float64
    initial_mass::Float64

    function Sim_ssp(id::Int64, type::String,
        x::Number, y::Number, z::Number,
        vx::Number, vy::Number, vz::Number,
        mass::Number, age::Number, metallicity::Number,
        initial_mass::Number)

        new(id, type,
        x, y, z,
        vx, vy, vz,
        mass, age, metallicity, initial_mass)
    end
end
