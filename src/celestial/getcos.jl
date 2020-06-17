# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using DelimitedFiles

"""
function getcos(ref::String)

    Returns cosref data for corresponding ref parameter
    Output is  [Ref, H0, omegaM, omegaL, omegaR, sigma8]
"""
function getcos(ref::String)

    ref = lowercase(ref)            #convert to all lower case letters
    int_ref = tryparse(Int64, ref)  #if possible, convert to integer

    if !isnothing(int_ref)
        ref = int_ref
    end

    data = readdlm(joinpath(@__DIR__, "data/cosref.tab"), header=true)
    row = findall(x -> x == ref, data[1][:,1])

    if length(row) == 0
        println(ref)
        error("Reference is not supported")
    end

    return data[1][row,:]  #Ref, H0, omegaM, omegaL, omegaR, sigma8
end
