# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

function bc_data()

h5open(joinpath(@__DIR__, "BC03lr.h5"), "r") do galaxy_file

    data = galaxy_file["BC03lr"]

    bcDict = Dict( "Age" => HDF5.readmmap(data["Age"]),
                    "AgeBins" => HDF5.readmmap(data["AgeBins"]),
                    "AgeWeights" => HDF5.readmmap(data["AgeWeights"]),
                    "Wave" => HDF5.readmmap(data["Wave"]),
                    "Z" => HDF5.readmmap(data["Z"]),
                    "Zevo" => Dict(
                        "1" => Dict("SFR" => HDF5.readmmap(data["Zevo"]["1"]["SFR"]),
                                    "SMgas" => HDF5.readmmap(data["Zevo"]["1"]["SMgas"]),
                                    "SMrem" => HDF5.readmmap(data["Zevo"]["1"]["SMrem"]),
                                    "SMstar" => HDF5.readmmap(data["Zevo"]["1"]["SMstar"]),
                                    "SMtot" => HDF5.readmmap(data["Zevo"]["1"]["SMtot"])),
                        "2" => Dict("SFR" => HDF5.readmmap(data["Zevo"]["2"]["SFR"]),
                                    "SMgas" => HDF5.readmmap(data["Zevo"]["2"]["SMgas"]),
                                    "SMrem" => HDF5.readmmap(data["Zevo"]["2"]["SMrem"]),
                                    "SMstar" => HDF5.readmmap(data["Zevo"]["2"]["SMstar"]),
                                    "SMtot" => HDF5.readmmap(data["Zevo"]["2"]["SMtot"])),
                        "3" => Dict("SFR" => HDF5.readmmap(data["Zevo"]["3"]["SFR"]),
                                    "SMgas" => HDF5.readmmap(data["Zevo"]["3"]["SMgas"]),
                                    "SMrem" => HDF5.readmmap(data["Zevo"]["3"]["SMrem"]),
                                    "SMstar" => HDF5.readmmap(data["Zevo"]["3"]["SMstar"]),
                                    "SMtot" => HDF5.readmmap(data["Zevo"]["3"]["SMtot"])),
                        "4" => Dict("SFR" => HDF5.readmmap(data["Zevo"]["4"]["SFR"]),
                                    "SMgas" => HDF5.readmmap(data["Zevo"]["4"]["SMgas"]),
                                    "SMrem" => HDF5.readmmap(data["Zevo"]["4"]["SMrem"]),
                                    "SMstar" => HDF5.readmmap(data["Zevo"]["4"]["SMstar"]),
                                    "SMtot" => HDF5.readmmap(data["Zevo"]["4"]["SMtot"])),
                        "5" => Dict("SFR" => HDF5.readmmap(data["Zevo"]["5"]["SFR"]),
                                    "SMgas" => HDF5.readmmap(data["Zevo"]["5"]["SMgas"]),
                                    "SMrem" => HDF5.readmmap(data["Zevo"]["5"]["SMrem"]),
                                    "SMstar" => HDF5.readmmap(data["Zevo"]["5"]["SMstar"]),
                                    "SMtot" => HDF5.readmmap(data["Zevo"]["5"]["SMtot"])),
                        "6" => Dict("SFR" => HDF5.readmmap(data["Zevo"]["6"]["SFR"]),
                                    "SMgas" => HDF5.readmmap(data["Zevo"]["6"]["SMgas"]),
                                    "SMrem" => HDF5.readmmap(data["Zevo"]["6"]["SMrem"]),
                                    "SMstar" => HDF5.readmmap(data["Zevo"]["6"]["SMstar"]),
                                    "SMtot" => HDF5.readmmap(data["Zevo"]["6"]["SMtot"]))),
                    "ZSpec" => Dict(
                                    1 => HDF5.readmmap(data["Zspec"]["1"]),
                                    2 => HDF5.readmmap(data["Zspec"]["2"]),
                                    3 => HDF5.readmmap(data["Zspec"]["3"]),
                                    4 => HDF5.readmmap(data["Zspec"]["4"]),
                                    5 => HDF5.readmmap(data["Zspec"]["5"]),
                                    6 => HDF5.readmmap(data["Zspec"]["6"])))
    return bcDict
    end
end
