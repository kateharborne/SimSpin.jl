using Distributions

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

                distribution = Normal(particle.obs.vy_obs, lsf_size^2)
                # adding the "gaussians" of each particle to the velocity bins
                cube[coord[1], coord[2], :] += cell_flux .* cdf.(distribution, vseq)[2:end]
            end
        end
    end
    return cube
end
