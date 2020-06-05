# Date created: 25/02/2020
# Author: Gerry Gralton

using SimSpin, Test, Random

@testset "telescopes" begin
    @testset "Circular aperture" begin
        tele = IFU(8,"circular", 7000, 2.63, 0.2, 1.25, "r")

        @test tele.sbin == 40
        @test size(tele.ap_region) == (tele.sbin, tele.sbin)    #aperture is same size as sbin
        @test sum(tele.ap_region) == 1264.
    end

    @testset "Square aperture" begin
        tele = IFU(8,"square", 7000, 2.63, 0.2, 1.25, "r")

        @test tele.sbin == 40
        @test size(tele.ap_region) == (tele.sbin, tele.sbin)    #aperture is same size as sbin
        @test sum(tele.ap_region) == tele.sbin^2
    end

    @testset "Hexagonal aperture" begin
        tele = IFU(8,"hexagonal", 7000, 2.63, 0.2, 1.25, "r")

        @test tele.sbin == 40
        @test size(tele.ap_region) == (tele.sbin, tele.sbin)    #aperture is same size as sbin
        @test sum(tele.ap_region) == 1140.
    end

    @testset "SAMI" begin
        tele = SAMI()

        @test tele.sbin == 30                                   #number of spatial bins
        @test size(tele.ap_region) == (tele.sbin, tele.sbin)    #aperture is same size as sbin
        @test sum(tele.ap_region) == 716.

        @test tele.vbinsize == 65.0
        @test tele.lsf_size ≈ 70.33 atol=0.01
    end

    @testset "MaNGA" begin
        tele = MaNGA()

        @test tele.sbin == 88                                   #number of spatial bins
        @test size(tele.ap_region) == (tele.sbin, tele.sbin)    #aperture is same size as sbin
        @test sum(tele.ap_region) == 5512.

        @test tele.vbinsize ≈ 72.73 atol=0.01
        @test tele.lsf_size ≈ 72.06 atol=0.01
    end

    @testset "CALIFA" begin
        tele = CALIFA()

        @test tele.sbin == 30                                   #number of spatial bins
        @test size(tele.ap_region) == (tele.sbin, tele.sbin)    #aperture is same size as sbin
        @test sum(tele.ap_region) == 644.

        @test tele.vbinsize ≈ 142.86 atol=0.01
        @test tele.lsf_size ≈ 171.38 atol=0.01
    end

    @testset "Hector" begin
        tele = Hector()

        @test tele.sbin == 18                                   #number of spatial bins
        @test size(tele.ap_region) == (tele.sbin, tele.sbin)    #aperture is same size as sbin
        @test sum(tele.ap_region) == 232.

        @test tele.vbinsize ≈ 3.125 atol=0.001
        @test tele.lsf_size ≈ 34.50 atol=0.01
    end
end
