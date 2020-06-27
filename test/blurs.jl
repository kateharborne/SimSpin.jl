# Date created: 25/02/2020
# Author: Gerry Gralton

using SimSpin, Test, Random

@testset "blurs" begin
    @testset "gaussian_blur" begin

        @test_throws ErrorException Gaussian_blur()
        @test_throws ErrorException Gaussian_blur(σ=0.1, fwhm=10.)

        σ = rand() + 1.
        gblur = Gaussian_blur(σ=σ)
        @test gblur.σ == σ
        @test isapprox(gblur.fwhm, σ * 2.355, rtol = 0.01)

        fwhm = 1 - rand() / 2
        gblur = Gaussian_blur(fwhm=fwhm)
        @test gblur.fwhm == fwhm
        @test isapprox(gblur.σ, fwhm / 2.355, rtol = 0.01)
    end

    @testset "moffat_blur" begin

        β = rand() * 2 + 3
        @test_throws ErrorException Moffat_blur(β)
        @test_throws ErrorException Moffat_blur(β, α=0.1, fwhm=100.)

        α = rand() + 1.
        mblur = Moffat_blur(β, α=α)
        @test mblur.α == α
        @test isapprox(mblur.fwhm, α * 2 * sqrt(2^(1/β) - 1), rtol = 0.01)

        fwhm = 1.5 - rand()
        mblur = Moffat_blur(β, fwhm=fwhm)
        @test mblur.fwhm == fwhm
        @test isapprox(mblur.α, fwhm / (2 * sqrt(2^(1/β) - 1)), rtol = 0.01)
    end
end
