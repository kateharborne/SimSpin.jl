# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using StaticArrays

function flux_grid(parts_in_cell::Array{Array{Galaxy_particle, 1}, 1},
                    ap_region::MArray,
                    sbin::Int64,
                    vbin::Int64,
                    redshift::Float64,
                    filter::String)

    flux = zeros(Float64, length(parts_in_cell))
    redshiftCoef = Lum2FluxFactor(z = redshift)

    for index in axes(parts_in_cell, 1)
        cell_flux = 0
        for particle in parts_in_cell[index]
            cell_flux += assign_flux(particle, filter, redshiftCoef)
        end
        flux[index] = cell_flux
    end

    flux_grid = reshape(flux, sbin, sbin, vbin)

    outside_ap = findall(x->x==0, ap_region)
    if(length(outside_ap) != 0)
        flux_grid[outside_ap[1], outside_ap[2], :] .= 0
    end

    return flux_grid
end
