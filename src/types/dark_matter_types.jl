abstract type Dark_matter end

struct Hernquist <: Dark_matter
    dm_mass::Float64
    dm_a::Float64

    """
        Hernquist(;dm_mass = 186.9,
                    dm_a = 34.5)

        Creates a Hernquist analytic dark matter mass potential.
        Default values are shown above.
    """
    function Hernquist(dm_mass::Float64 = 184.9,
                        dm_a::Float64 = 34.5)
        new(dm_vm, dm_a)
    end
end

struct NFW <: Dark_matter
    dm_vm::Float64
    dm_a::Float64
    dm_rhof::Float64

    """
        NFW(;dm_vm = 186.9,
                dm_a = 34.5,
                dm_rhof = 0.035)

        Creates an NFW analytic dark matter mass potential.
        Default values are shown above.
    """
    function NFW(dm_vm::Float64 = 186.9,
                    dm_a::Float64 = 34.5,
                    dm_rhof::Float64 = 0.035)
        new(dm_vm, dm_a, dm_rhof)
    end
end


function Dark_matter(profile::String;
                        dm_mass::Float64 = 184.9,
                        dm_vm::Float64=186.9,
                        dm_a::Float64=34.5,
                        dm_rhof::Float64=0.035)
    if lowercase(profile) == "hernquist"
        return Hernquist(dm_vm, dm_a)
    elseif lowercase(profile) == "nfw"
        return NFW(dm_vm, dm_a, dm_rhof)
    else
        error("Dark matter profile: ", profile, "is unsupported. Please specify 'Hernquist' or 'NFW'.")
    end
end
