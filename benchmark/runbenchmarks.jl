using ConvertBNG
using BenchmarkTools

function benchmark_convert_bng()
    # London bounding box
    N = 51.691874116909894
    E = 0.3340155643740321
    S = 51.28676016315085
    W = -0.5103750689005356

    num_coords = 1000000
    lons = rand(W:0.01:E, num_coords, 1)
    lats = rand(S:0.01:N, num_coords, 1)

    lons_lats = hcat(lons, lats)
    convert_osgb36(lons_lats)
end

show(@benchmark benchmark_convert_bng())