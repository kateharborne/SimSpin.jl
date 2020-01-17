# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using StaticArrays

"""
    flux_grid(parts_in_cell,
                ap_region,
                sbin,
                vbin,
                redshift,
                filter)

Computes the fluxes for each element of the IFU data-cube.

The purpose of this function is to construct the mock flux values within each cell of the IFU
cube. It accepts output parameters from `obs_data_prep()` and returns a 3D array containing
the flux at each cell position due to contributions from the particles. If ssd particles are supplied,
an SED is generated in each cell using ProSpect. Else, the luminosity in each cell is converted to flux.

Parameters:\n
    parts_in_cell       1D array of the particles corresponding to each element in the IFU data-cube.
    ap_region           The aperture region mask used to remove flux outside of the specified aperture.
    sbin                The number of spatial bins in the aperture.
    vbin                The number of velocity bins in the flux grid.
    redshift            The projected redshift at which the observation is made.
    filter              If ssp particles are supplied, the filter within which the SED is generated.
                        Options include "r" and "g"  for SDSS-r and SDSS-g bands respectively.
"""
function flux_grid(parts_in_cell::Array{Array{Galaxy_particle, 1}, 1},
                    ap_region::MArray,
                    sbin::Int64,
                    vbin::Int64,
                    redshift::Float64,
                    filter::String)

    flux = zeros(Float64, length(parts_in_cell))
    redshiftCoef = Lum2FluxFactor(z = redshift)

    filter = get_filter(filter)

    for index in axes(parts_in_cell, 1)
        cell_flux = 0
        for particle in parts_in_cell[index]
            cell_flux += assign_flux(particle, filter, redshiftCoef)
        end
        flux[index] = cell_flux
    end

    flux_grid = reshape(flux, sbin, sbin, vbin)

    for bin = 1:vbin #Set all cells outside aperture to zero
        flux_grid[:,:,bin] = flux_grid[:,:,bin] .* ap_region
    end

    return flux_grid
end
