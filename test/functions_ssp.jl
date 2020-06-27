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
            @test observe.sbin      == tele.sbin
        @test observe.lsf_size      == tele.lsf_size

        @test size(observe.vseq)[1]== observe.vbin + 1        #number of velocity bins matches bounds
    end

    @testset "obs_data_prep_with_blur" begin
        filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_SSP.hdf5")
        particles = sim_data(filename, ssp=true)

        filter = rand(["r"; "g"])
        tele = IFU(8,"square", 7000, 2.63, 0.2, 1.25, filter)
        blur = Gaussian_blur(σ = rand() + 2)

        z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
        envir = Environment(z, inc_deg, r200, blur)

        galaxy_data, parts_in_cell, observe = obs_data_prep(particles, tele, envir)

        @test size(galaxy_data)[1]  == sum(first.(size.(parts_in_cell)))     #output particle array has same number of particles as datacube cells
        @test size(parts_in_cell)[1]== tele.sbin * tele.sbin * observe.vbin   #number of cells in datacube

        @test observe.z             == envir.z
        @test observe.inc_deg       == envir.inc_deg
        @test observe.r200          == envir.r200
        @test observe.blur          == blur
        @test observe.mass2light    == envir.mass2light

        @test observe.ap_region     == tele.ap_region
        @test observe.sbin          == tele.sbin
        @test observe.lsf_size      == tele.lsf_size

        @test size(observe.vseq)[1] == observe.vbin + 1        #number of velocity bins matches bounds
    end

    @testset "build_datacube" begin
        filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_SSP.hdf5")
        particles = sim_data(filename, ssp=true)

        z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
        envir = Environment(z, inc_deg, r200)

        tele_no_filt = SAMI()
        @test_throws Exception build_datacube(particles, tele_no_filt, envir)

        filt = rand(["r"; "g"])
        tele = SAMI(filter=filt)
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

    @testset "build_datacube_with_blur" begin
        filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_SSP.hdf5")
        particles = sim_data(filename, ssp=true)

        filter = rand(["r"; "g"])
        tele = IFU(8,"hexagonal", 7000, 2.63, 0.2, 1.25, filter)

        z = rand(); inc_deg = rand(0:90); r200 = rand(100:200); m2l = (rand() + 1, rand() * 0.5 + 0.5);
        blur = Moffat_blur(rand() * 3 + 2, fwhm=rand() * 3 + 1);
        envir = Environment(z, inc_deg, r200, m2l, blur)

        datacube, observe = build_datacube(particles, tele, envir)

        @test size(datacube) == (tele.sbin, tele.sbin, observe.vbin)

        aperture = repeat(tele.ap_region, outer=[1,1,observe.vbin])
        zero_indices = findall(x-> x == 0, aperture)
        @test all(datacube[zero_indices] .== 0.0)      #no data outside telescope aperture

        @test observe.z             == envir.z
        @test observe.inc_deg       == envir.inc_deg
        @test observe.r200          == envir.r200
        @test observe.blur          == blur
        @test observe.mass2light    == m2l

        @test observe.ap_region     == tele.ap_region
        @test observe.sbin          == tele.sbin
        @test observe.lsf_size      == tele.lsf_size

        @test size(observe.vseq)[1]== observe.vbin + 1        #number of velocity bins matches bounds
    end

    @testset "build_datacube_with_environment_array" begin
        filename = joinpath(dirname(pathof(SimSpin)), "..", "data", "SimSpin_SSP.hdf5")
        particles = sim_data(filename, ssp=true)

        filter = rand(["r"; "g"])
        tele = IFU(12,"hexagonal", 2300, 2.63, 0.2, 1.25, filter)

        z = [rand(), rand(), rand()]; inc_deg = rand(0:90); r200 = rand(100:200); m2l = (rand() * 0.2 + 0.4, rand() + 1);
        blur = [Moffat_blur(rand() *2 + 3, fwhm=rand() + 1),
                Moffat_blur(rand() *2 + 3, α=rand() + 1),
                Gaussian_blur(σ = rand() * 3)];

        envirs = Environment(z, inc_deg, r200, m2l, blur)
        @test length(envirs) == length(z) * length(blur)

        data_array = build_datacube(particles, tele, envirs)
        @test length(data_array) == length(envirs)
    end
end
