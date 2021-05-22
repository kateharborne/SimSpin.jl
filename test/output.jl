# Date created: 03/08/2020
# Author: Gerry Gralton

using SimSpin, Test, Random
using FITSIO

@testset "output" begin
    particles = sim_data()
    tele = SAMI()

    z = rand() + 0.1; inc_deg = rand(0:90); r200 = rand(100:200);
    envir = Environment(z, inc_deg, r200)

    datacube, observe = build_datacube(particles, tele, envir)

    #Test default output function
    sim_FITS(datacube, observe, "test")
    file = FITS("test.fits")
    read_data = read(file[1])
    @test read_data == datacube

    #Test output function with build_datacube output as tuple
    sim_FITS((datacube, observe), "test_tuple")
    file = FITS("test_tuple.fits")
    read_data = read(file[1])
    @test read_data == datacube

    #Test output function with array of build_datacube outputs
    z_arr = rand(Float64, 4)  .+ 0.1
    envir_arr = Environment(z_arr, inc_deg, r200)
    datacube_arr = build_datacube(particles, tele, envir_arr)

    sim_FITS(datacube_arr, folder = "", file_prefix="test")

    for ind in 1:length(datacube_arr)
        file = FITS(string("test_z", envir_arr[ind].z,
                            "_INC", envir_arr[ind].inc_deg,
                            "_R", envir_arr[ind].r200,
                            "_M2L", envir_arr[ind].mass2light, ".fits"))
        read_data = read(file[1])
        @test read_data == datacube_arr[ind][1]
    end

end
