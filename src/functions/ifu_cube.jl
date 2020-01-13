# Date created: 12/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using Distributions

"""
    ifu_cube(flux_grid,
                parts_in_cell,
                sbin,
                vbin,
                vseq,
                lsf_size)

The purpose of this function is to construct an IFU data cube. It accepts the flux_grid
output by the `flux_grid()` function and returns a similar, IFU-like, 3D array
where each particle's flux contributes a Gaussian distribution in the velocity axis.

Parameters:\n
    flux_grid       Flux grid output by `flux_grid()`
    parts_in_cell   1D array of the particles corresponding to each element in the IFU data-cube.
    sbin            The number of spatial bins in the aperture.
    vbin            The number of velocity bins in the flux grid.
    vseq            The bounds of each velocity bin in the flux grid.
    lsf_size        The Gaussian standard deviation of the line spread function in km/s.
"""
function ifu_cube(flux_grid::Array{Float64, 3},
                    parts_in_cell::Array{Array{SimSpin.Galaxy_particle,1},1},
                    sbin::Int64,
                    vbin::Int64,
                    vseq::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},
                    lsf_size::Float64)

    cube = zeros(Float64, (sbin, sbin, vbin))

    for (index, cell) in enumerate(parts_in_cell) # for each cell in cube
        coord = [index % sbin, index % (sbin^2) ÷ sbin + 1, index ÷ (sbin^2) + 1]
        if(index % sbin == 0 && index % (sbin^2) == 0 )
            coord = [sbin, sbin, index ÷ (sbin^2)]
        elseif(index % sbin == 0 && index % (sbin^2) != 0 )
            coord = [sbin, index % (sbin^2) ÷ sbin, index ÷ (sbin^2) + 1]
        end

        num_of_parts = length(cell) # number of particles in that cell

        if num_of_parts > 0  # if particles in that cell
            for particle in cell
                cell_mass = sum(getfield.(cell, :mass))
                cell_flux = flux_grid[coord[1], coord[2], coord[3]] * particle.mass / cell_mass

                distribution = Normal(particle.obs.vy_obs, lsf_size^2)  #Normal distribution of particle's vy_obs
                # adding the "gaussians" of each particle to the velocity bins
                cube[coord[1], coord[2], :] += diff(cell_flux .* cdf.(distribution, vseq))
            end
        end
    end
    
    return cube
end
