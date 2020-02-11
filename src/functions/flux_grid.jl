# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    flux_grid(parts_in_cell,
                observe,
                filter)

Computes the fluxes for each element of the IFU data-cube.

The purpose of this function is to construct the mock flux values within each cell of the IFU
cube. It accepts output parameters from `obs_data_prep()` and returns a 3D array containing
the flux at each cell position due to contributions from the particles. If ssd particles are supplied,
an SED is generated in each cell using ProSpect. Else, the luminosity in each cell is converted to flux.

Parameters:\n
    parts_in_cell       1D array of the particles corresponding to each element in the IFU data-cube.
    observe             Struct of type `Observation` containing all observation parameters.
    filter              If ssp particles are supplied, the filter within which the SED is generated.
                        Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
"""
function flux_grid(parts_in_cell::Array{Array{Galaxy_particle, 1}, 1},
                    observe::Observation,
                    filter::Union{String, Nothing})

    flux = zeros(Float64, length(parts_in_cell))
    redshiftCoef = Lum2FluxFactor(z=observe.z)

    if !isnothing(filter)
        filter = get_filter(filter)
    end

    Threads.@threads for index in axes(parts_in_cell, 1)
        cell_flux = 0
        for particle in parts_in_cell[index]
            cell_flux += assign_flux(particle, filter, redshiftCoef)
        end
        flux[index] = cell_flux
    end

    flux_grid = reshape(flux, observe.sbin, observe.sbin, observe.vbin)

    for bin = 1:observe.vbin #Set all cells outside aperture to zero
        flux_grid[:,:,bin] = flux_grid[:,:,bin] .* observe.ap_region
    end

    return flux_grid
end
