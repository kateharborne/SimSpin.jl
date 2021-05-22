# Date created: 25/02/2020
# Author: Gerry Gralton

using SimSpin, Test, Random

@testset "environments" begin
    gauss_blur=Gaussian_blur(σ=1.)
    moffat_blur=Moffat_blur(1.5, α=1.)

    z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
    m2l = 1.
    m2l_tup = (1.2, 1.5)
    environment = Environment(z, inc_deg, r200)

    @test environment.mass2light == m2l
    @test isnothing(environment.blur)

    environment = Environment(z, inc_deg, r200, gauss_blur)
    @test !isnothing(environment.blur)

    environment = Environment(z, inc_deg, r200, m2l, gauss_blur)
    @test !isnothing(environment.blur)

    environment = Environment(z, inc_deg, r200, m2l_tup, gauss_blur)
    @test !isnothing(environment.blur)
    @test length(environment.mass2light) == 2

    z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
    m2l = 1.3
    environment = Environment(z, inc_deg, r200, m2l)

    @test environment.mass2light == m2l

    z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
    m2l = (1.3, 2.5)
    environment = Environment(z, inc_deg, r200, m2l)

    @test environment.mass2light[1] == m2l[1]
    @test environment.mass2light[2] == m2l[2]

    z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
    m2l_arr_tup = [(1.3, 2.5), (1., 1.)]
    environment = Environment(z, inc_deg, r200, m2l_arr_tup)

    @test environment[1].mass2light == m2l_arr_tup[1]
    @test environment[2].mass2light == m2l_arr_tup[2]

    z = rand(0.1:0.1:1.2)
    zlength = rand(2.:10.)
    zmax = round(z * zlength, digits=2)
    inc_deg_arr = [30; 45; 60]

    environment_arr = Environment(z:z:zmax, inc_deg_arr, r200)
    @test length(environment_arr) == zlength * length(inc_deg_arr)

    environment_arr_m2l= Environment(z:z:zmax, inc_deg_arr, r200, m2l_arr_tup)
    @test length(environment_arr_m2l) == 2 * length(environment_arr)

    m2l_arr = [1.5, 2.5, 1.]
    environment_arr = Environment(z:z:zmax, inc_deg_arr, r200, m2l_arr, gauss_blur)
    @test length(environment_arr) == zlength * length(inc_deg_arr) * length(m2l_arr)

    environment_arr = Environment(z:z:zmax, inc_deg_arr, r200, m2l_arr)
    @test length(environment_arr) == zlength * length(inc_deg_arr) * length(m2l_arr)

    environment_arr = Environment(z:z:zmax, inc_deg_arr, r200, m2l_arr_tup, gauss_blur)
    @test length(environment_arr) == zlength * length(inc_deg_arr) * length(m2l_arr_tup)

    blur_array = [moffat_blur, gauss_blur]
    environment_arr = Environment(z:z:zmax, inc_deg_arr, r200, blur_array)
    @test length(environment_arr) == zlength * length(inc_deg_arr) * length(blur_array)


end
