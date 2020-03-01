using ConvertBNG
using BenchmarkTools

function benchmark_convert_bng()
    # UK bounding box
    N = 60.83
    E =  1.77
    S = 49.93
    W = -8.57

    num_coords = 1000000
    lons = rand(W:0.01:E, num_coords, 1)
    lats = rand(S:0.01:N, num_coords, 1)

    lons_lats = hcat(lons, lats)
    convert_osgb36(lons_lats)
end

show(@benchmark benchmark_convert_bng())