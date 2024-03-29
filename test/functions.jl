# Date created: 25/02/2020
# Author: Gerry Gralton

using SimSpin, Test, Random


@testset "functions_without_SSP" begin
    @testset "sim_data" begin
        particles = sim_data()
        particles3 = sim_data(ptype=[3])

        @test length(particles) == 25000
        @test length(particles3) == 15000

        @test_throws ErrorException sim_data(ptype=[4])
        @test_throws ErrorException sim_data(ssp=true)
    end

    @testset "sim_to_galaxy" begin
        particles = sim_data()
        particles_galaxy = sim_to_galaxy(particles)
        @test particles_galaxy != particles

        warning = (:warn, "This galaxy was already in the galaxy reference frame.")
        @test_logs warning sim_to_galaxy(particles_galaxy)

        centre = [(rand() + 1) * rand(-100.:10.:100),
                    (rand() + 1) * rand(-100.:10.:100),
                    (rand() + 1) * rand(-100.:10.:100)]
        vel_centre = [(rand() + 1) * rand(-100.:10.:100),
                        (rand() + 1) * rand(-100.:10.:100),
                        (rand() + 1) * rand(-100.:10.:100)]
        tot_centre = vcat(centre, vel_centre)
        @test_logs warning sim_to_galaxy(particles_galaxy, centre=centre)
        @test_logs warning sim_to_galaxy(particles_galaxy, centre=tot_centre)

        particles_centre = sim_to_galaxy(particles, centre=centre)
        @test particles_centre != particles_galaxy

        particles_tot_centre = sim_to_galaxy(particles, centre=vel_centre)
        @test particles_tot_centre != particles_galaxy
        @test particles_tot_centre != particles_centre

        @test_throws ErrorException sim_to_galaxy(particles, centre=vel_centre[1:end-1])
    end

    @testset "obs_data_prep" begin
        particles = sim_data()
        tele = SAMI()

        z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
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

        @test size(observe.vseq)[1] == observe.vbin + 1        #number of velocity bins matches bounds
    end

    @testset "obs_data_prep_with_blur" begin
        particles = sim_data()
        tele = IFU(8,"square", 7000, 2.63, 0.2, 1.25)
        blur = Gaussian_blur(fwhm=0.5)

        z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
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
        particles = sim_data()
        tele = SAMI()

        z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
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

        @test size(observe.vseq)[1]== observe.vbin + 1        #number of velocity bins matches bounds
    end

    @testset "build_datacube_with_blur" begin
        particles = sim_data()
        tele = IFU(8,"hexagonal", 7000, 2.63, 0.2, 1.25)

        z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
        m2l = (1.2, 0.8); blur = Moffat_blur(5., fwhm=1.5);
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
        particles = sim_data()
        tele = IFU(12,"hexagonal", 2300, 2.63, 0.2, 1.25)

        z = [rand(), rand(), rand()] .+ 0.1; inc_deg = rand(0:90); r200 = rand(100:200); m2l = (1.2, 0.8);
        blur = [Moffat_blur(rand() *2 + 3, fwhm=rand() + 1),
                Moffat_blur(rand() *2 + 3, α=rand() + 1),
                Gaussian_blur(σ = rand() * 3)];
        envirs = Environment(z, inc_deg, r200, m2l, blur)
        @test length(envirs) == length(z) * length(blur)

        data_array = build_datacube(particles, tele, envirs)

        @test length(data_array) == length(envirs)
    end
end
