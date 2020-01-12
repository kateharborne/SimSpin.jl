# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    Particle

Struct containing information about an individual particle.
A particle always has a spatial location, velocity and mass
but may also have a number of other attributes eg observed radius. 
"""
abstract type Particle end

"""
    Sim_particle

A type of Particle. All values are in simulation reference frame, ie taken with respect
to an arbitrary location in the simulation.
"""
abstract type Sim_particle <: Particle end

"""
    Galaxy_particle

A type of Particle. All values are in galaxy reference frame, ie taken with respect
to the median spatial and velocity values of all particles in the galaxy.
"""
abstract type Galaxy_particle <: Particle end
