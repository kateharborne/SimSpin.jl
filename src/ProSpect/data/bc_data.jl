# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

function bc_data()

h5open(joinpath(@__DIR__, "BC03lr.h5"), "r") do galaxy_file

    data = galaxy_file["BC03lr"]

    bcDict = Dict( "Age" => readmmap(data["Age"]),
                    "AgeBins" => readmmap(data["AgeBins"]),
                    "AgeWeights" => readmmap(data["AgeWeights"]),
                    "Wave" => readmmap(data["Wave"]),
                    "Z" => readmmap(data["Z"]),
                    "Zevo" => Dict(
                        "1" => Dict("SFR" => readmmap(data["Zevo"]["1"]["SFR"]),
                                    "SMgas" => readmmap(data["Zevo"]["1"]["SMgas"]),
                                    "SMrem" => readmmap(data["Zevo"]["1"]["SMrem"]),
                                    "SMstar" => readmmap(data["Zevo"]["1"]["SMstar"]),
                                    "SMtot" => readmmap(data["Zevo"]["1"]["SMtot"])),
                        "2" => Dict("SFR" => readmmap(data["Zevo"]["2"]["SFR"]),
                                    "SMgas" => readmmap(data["Zevo"]["2"]["SMgas"]),
                                    "SMrem" => readmmap(data["Zevo"]["2"]["SMrem"]),
                                    "SMstar" => readmmap(data["Zevo"]["2"]["SMstar"]),
                                    "SMtot" => readmmap(data["Zevo"]["2"]["SMtot"])),
                        "3" => Dict("SFR" => readmmap(data["Zevo"]["3"]["SFR"]),
                                    "SMgas" => readmmap(data["Zevo"]["3"]["SMgas"]),
                                    "SMrem" => readmmap(data["Zevo"]["3"]["SMrem"]),
                                    "SMstar" => readmmap(data["Zevo"]["3"]["SMstar"]),
                                    "SMtot" => readmmap(data["Zevo"]["3"]["SMtot"])),
                        "4" => Dict("SFR" => readmmap(data["Zevo"]["4"]["SFR"]),
                                    "SMgas" => readmmap(data["Zevo"]["4"]["SMgas"]),
                                    "SMrem" => readmmap(data["Zevo"]["4"]["SMrem"]),
                                    "SMstar" => readmmap(data["Zevo"]["4"]["SMstar"]),
                                    "SMtot" => readmmap(data["Zevo"]["4"]["SMtot"])),
                        "5" => Dict("SFR" => readmmap(data["Zevo"]["5"]["SFR"]),
                                    "SMgas" => readmmap(data["Zevo"]["5"]["SMgas"]),
                                    "SMrem" => readmmap(data["Zevo"]["5"]["SMrem"]),
                                    "SMstar" => readmmap(data["Zevo"]["5"]["SMstar"]),
                                    "SMtot" => readmmap(data["Zevo"]["5"]["SMtot"])),
                        "6" => Dict("SFR" => readmmap(data["Zevo"]["6"]["SFR"]),
                                    "SMgas" => readmmap(data["Zevo"]["6"]["SMgas"]),
                                    "SMrem" => readmmap(data["Zevo"]["6"]["SMrem"]),
                                    "SMstar" => readmmap(data["Zevo"]["6"]["SMstar"]),
                                    "SMtot" => readmmap(data["Zevo"]["6"]["SMtot"]))),
                    "Zspec" => Dict(
                        "1" => readmmap(data["Zspec"]["1"]),
                        "2" => readmmap(data["Zspec"]["2"]),
                        "3" => readmmap(data["Zspec"]["3"]),
                        "4" => readmmap(data["Zspec"]["4"]),
                        "5" => readmmap(data["Zspec"]["5"]),
                        "6" => readmmap(data["Zspec"]["6"])))
    return bcDict
    end
end
