using BenchmarkTools
using ConvertBNG
using Printf

const NUM_COORDS = 1_000_000
const REPEAT_RUNS = 50
const LONS_LATS = Ref{Array{Float64,2}}()
const ENS = Ref{Array{Float64,2}}()

function initial_setup(num_coords)
    # London bounding box
    N = 51.691874117
    E =  0.334015564
    S = 51.286760163
    W = -0.510375069
    
    # Create random testset of lon, lats
    lons = rand(W:0.01:E, num_coords, 1)
    lats = rand(S:0.01:N, num_coords, 1)

    # Create random testset of eastings, northings
    e_max, n_max = convert_bng(E, N)
    e_min, n_min = convert_bng(W, S)
    eastings = rand(e_min:0.1:e_max, num_coords, 1)
    northings = rand(n_min:0.1:n_max, num_coords, 1)

    return hcat(lons, lats), hcat(eastings, northings)
end

function benchmark(num_coords::Integer, repeat_runs::Integer)
    # Run benchmark for convert_bng
    println("Running convert_bng benchmark for $repeat_runs repeat runs of $num_coords points...")
    starttime = time()
    for x in 1:repeat_runs
        convert_bng(LONS_LATS[])
    end
    totaltime = @sprintf("%.2f", time() - starttime)
    println("Time taken: $totaltime s")
    
    # Run benchmark for convert_lonlat
    println("Running convert_lonlat benchmark for $repeat_runs repeat runs of $num_coords points...")
    starttime = time()
    for x in 1:repeat_runs
        convert_lonlat(ENS[])
    end
    totaltime = @sprintf("%.2f", time() - starttime)
    println("Time taken: $totaltime s")
    print("Benchmarks complete!")
end

# Setup test coordinate arrays
LONS_LATS[], ENS[] = initial_setup(NUM_COORDS)

# Run benchmarks
benchmark(NUM_COORDS, REPEAT_RUNS)
