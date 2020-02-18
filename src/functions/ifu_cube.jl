# Date created: 12/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using Distributions

"""
    ifu_cube(flux_grid,
                parts_in_cell,
                observe)

The purpose of this function is to construct an IFU data cube. It accepts a flux grid
in the format output by the `flux_grid()` function and returns a similar, IFU-like, 3D array
where each particle's flux contributes a Gaussian distribution in the velocity axis.

Parameters:\n
    flux_grid       Flux grid output by `flux_grid()`
    parts_in_cell   1D array of the particles corresponding to each element in the IFU data-cube.
    observe         Struct of type `Observation` containing all observation parameters.
"""
function ifu_cube(flux_grid::Array{Float64, 3},
                    parts_in_cell::Array{Array{SimSpin.Galaxy_particle,1},1},
                    observe::Observation)

    cube = zeros(Float64, (observe.sbin, observe.sbin, observe.vbin))

    sbin2 = observe.sbin^2

    Threads.@threads for index in 1:length(parts_in_cell) # for each cell in cube

        cell = parts_in_cell[index]
        
        coord = [index % observe.sbin, index % (sbin2) ÷ observe.sbin + 1, index ÷ (sbin2) + 1]
        if(index % observe.sbin == 0 && index % (sbin2) == 0 )
            coord = [observe.sbin, observe.sbin, index ÷ (sbin2)]
        elseif(index % observe.sbin == 0 && index % (sbin2) != 0 )
            coord = [observe.sbin, index % (sbin2) ÷ observe.sbin, index ÷ (sbin2) + 1]
        end

        num_of_parts = length(cell) # number of particles in that cell

        if num_of_parts > 0  # if particles in that cell
            cell_mass = sum(getfield.(cell, :mass))

            for particle in cell
                cell_flux = flux_grid[coord[1], coord[2], coord[3]] * particle.mass / cell_mass

                distribution = Normal(particle.obs.vy_obs, observe.lsf_size)  #Normal distribution of particle's vy_obs
                # adding the "gaussians" of each particle to the velocity bins
                cube[coord[1], coord[2], :] .+= diff(cell_flux .* cdf.(distribution, observe.vseq))
            end
        end
    end

    return cube
end
