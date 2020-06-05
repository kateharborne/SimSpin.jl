# Date created: 25/02/2020
# Author: Gerry Gralton

using SimSpin, Test, Random

@testset "functions_with_SSP" begin
    @testset "sim_data" begin
        filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_SSP.hdf5")
        particles = sim_data(filename, ssp=true)
        particles_dm = sim_data(filename, ptype=[1])

        @test length(particles) == 311179
        @test length(particles_dm) == 190516

        @test_throws ErrorException sim_data(filename, ptype=[3])

        warning = (:warn, "SSP data is available for Star particles in this simulation file but has not been read in. If spectra is desired set ssp=true.")
        @test_logs warning sim_data(filename, ssp=false)
    end

    @testset "obs_data_prep" begin
        filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_SSP.hdf5")
        particles = sim_data(filename, ssp=true)

        filt = rand(["r"; "g"])
        tele = SAMI(filter=filt)

        z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
        envir = Environment(z, inc_deg, r200)

        galaxy_data, parts_in_cell, observe = obs_data_prep(particles, tele, envir)

        @test size(galaxy_data)[1]  == sum(first.(size.(parts_in_cell)))     #output particle array has same number of particles as datacube cells
        @test size(parts_in_cell)[1]== tele.sbin * tele.sbin * observe.vbin   #number of cells in datacube

        @test observe.z             == envir.z
        @test observe.inc_deg       == envir.inc_deg
        @test observe.r200          == envir.r200
        @test observe.blur          == envir.blur
        @test observe.mass2light    == envir.mass2light

        @test observe.ap_region     == tele.ap_region
        @test observe.sbin          == tele.sbin
        @test observe.lsf_size      == tele.lsf_size

        @test size(observe.vseq)[1] == observe.vbin + 1        #number of velocity bins matches bounds
    end

    @testset "build_datacube" begin
        filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_SSP.hdf5")
        particles = sim_data(filename, ssp=true)

        filt = rand(["r"; "g"])
        tele = SAMI(filter=filt)

        z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
        envir = Environment(z, inc_deg, r200)

        datacube, observe = build_datacube(particles, tele, envir)

        @test size(datacube) == (tele.sbin, tele.sbin, observe.vbin)

        aperture = repeat(tele.ap_region, outer=[1,1,observe.vbin])
        zero_indices = findall(x-> x == 0, aperture)
        @test all(datacube[zero_indices] .== 0.0)      #no data outside telescope aperture

        @test observe.z             == envir.z
        @test observe.inc_deg       == envir.inc_deg
        @test observe.r200          == envir.r200
        @test observe.blur          == envir.blur
        @test observe.mass2light    == envir.mass2light

        @test observe.ap_region     == tele.ap_region
        @test observe.sbin          == tele.sbin
        @test observe.lsf_size      == tele.lsf_size

        @test size(observe.vseq)[1] == observe.vbin + 1        #number of velocity bins matches bounds
    end
end
