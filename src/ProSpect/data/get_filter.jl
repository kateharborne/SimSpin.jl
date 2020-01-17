using Interpolations
using DelimitedFiles

"""
    get_filter(filter)

Reads in ProSpect filter specified by parameter `filter`. Converts to
object of type `Interpolations.FilledExtrapolation`.

Filter arguments currently supported are `r` and `g` corresponding to
filt_r_SDSS.tab and filt_g_SDSS.tab respectively.
"""
function get_filter(filter::String)

    if(filter == "g")   tempfilt = readdlm(joinpath(@__DIR__, "filt_g_SDSS.tab"), header=true)
    elseif(filter == "r")   tempfilt = readdlm(joinpath(@__DIR__, "filt_r_SDSS.tab"), header=true)
    else error("An unsupported filter type was specified.")
    end

    filt_wave = tempfilt[1][:,1]
    filt_response = tempfilt[1][:,2]

    itp = interpolate((filt_wave,),
                        abs.(filt_response),
                        Gridded(Linear()))
    extrap = extrapolate(itp, 0) #return zero for all of wave outside of filter bounds

    return extrap
end
