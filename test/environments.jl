# Date created: 25/02/2020
# Author: Gerry Gralton

using SimSpin, Test, Random

@testset "environments" begin
    gauss_blur=Gaussian_blur(sigma=1.)
    moffat_blur=Moffat_blur(1.5, α=1.)

    z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
    m2l = 1.
    environment = Environment(z, inc_deg, r200)

    @test environment.mass2light == m2l
    @test isnothing(environment.blur)

    environment = Environment(z, inc_deg, r200, gauss_blur)
    @test !isnothing(environment.blur)

    environment = Environment(z, inc_deg, r200, m2l, gauss_blur)
    @test !isnothing(environment.blur)

    z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
    m2l = 1.3
    environment = Environment(z, inc_deg, r200, m2l)

    @test environment.mass2light == m2l

    z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
    m2l = (1.3, 2.5)
    environment = Environment(z, inc_deg, r200, m2l)

    @test environment.mass2light[1] == m2l[1]
    @test environment.mass2light[2] == m2l[2]

    z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
    m2l = [(1.3, 2.5), (1., 1.)]
    environment = Environment(z, inc_deg, r200, m2l)

    @test environment[1].mass2light == m2l[1]
    @test environment[2].mass2light == m2l[2]

    zlength = rand(1:10)
    zmax = z * zlength
    inc_deg_arr = [30; 45; 60]

    environment_arr = Environment(z:z:zmax, inc_deg_arr, r200)
    @test length(environment_arr) == zlength * length(inc_deg_arr)

    environment_arr_m2l= Environment(z:z:zmax, inc_deg_arr, r200, m2l)
    @test length(environment_arr_m2l) == 2 * length(environment_arr)
end
