# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using Interpolations

function interp_param(x::Array{Float64,1},
                        params::Array{Float64,1};
                        log::Bool=false,
                        method::String="linear")

    flag = fill(2, length(x))

    minp = minimum(params)
    min_index = findall(z->z < minp, x)
    maxp = maximum(params)
    max_index = findall(z->z > maxp, x)

    flag[min_index] .= 1
    flag[max_index] .= 3

    x[min_index] .= minp
    x[max_index] .= maxp

    if(log && method=="linear")
        itp = interpolate((log10.(params),),
                        1:length(params),
                        Gridded(Linear()))
        temp = itp.(log10.(x))

    elseif(!log && method=="linear")
        itp = interpolate((params,),
                        1:length(params),
                        Gridded(Linear()))
        temp = itp.(x)
    else
        error("Unsupported interpolation type. Use linear interpolation.")
    end

    ind_zero = findall(x->mod(x, 1) == 0, temp)
    ind_two = findall(x-> x == 2, temp)
    flag[intersect(ind_zero, ind_two)] .= 0

    ID_lo::Array{Int64, 1} = floor.(temp)
    ID_hi::Array{Int64, 1} = ceil.(temp)

    param_lo = params[ID_lo]
    param_hi = params[ID_hi]

    ID_mode = copy(ID_lo)
    ind_half = findall(x-> mod.(x, 1) > 0.5, temp)
    ID_mode[ind_half] .= ID_hi[ind_half]

    weight_hi = mod.(temp, 1)
    weight_lo = 1 .- weight_hi

    return Dict(
        "x"=>x,
        "param_lo"=>param_lo,
        "param_hi"=>param_hi,
        "ID_lo"=>ID_lo,
        "ID_hi"=>ID_hi,
        "ID_mode"=>ID_mode,
        "weight_lo"=>weight_lo,
        "weight_hi"=>weight_hi,
        "flag"=>flag
    )
end

function qdiff(vect::Array{Float64, 1}, pad0=true)
    if pad0
        return vcat(0, vect[2:end] .- vect[1:(end-1)])
    else
        return vect[2:end].-vect[1:(end-1)]
    end
end
