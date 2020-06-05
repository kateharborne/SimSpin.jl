# Date created: 25/02/2020
# Author: Gerry Gralton

using SimSpin, Test, Random

@testset "environments" begin
    z = rand(); inc_deg = rand(0:90); r200 = rand(100:200);
    environment = Environment(z, inc_deg, r200)

    @test environment.mass2light == 1.
    @test isnothing(environment.blur)

    zlength = rand(1:10)
    zmax = z * zlength
    inc_deg_arr = [30; 45; 60]

    environment_arr = Environment(z:z:zmax, inc_deg_arr, r200)
    @test length(environment_arr) == zlength * length(inc_deg_arr)
end
