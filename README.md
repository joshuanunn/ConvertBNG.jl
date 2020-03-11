## Introduction
A Julia library for fast conversion between WGS84 longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates. Conversions use a standard 7-element Helmert transform with the addition of OSTN15 corrections for [accuracy](#accuracy).

This is a derived work of original repositories written by Stephan Hügel, which include a Rust binary [lonlat_bng](https://github.com/urschrei/lonlat_bng) and the corresponding Python utility library [convertbng](https://github.com/urschrei/convertbng). As Python is relatively slow, [convertbng](https://github.com/urschrei/convertbng) utilises a [Rust binary](https://github.com/urschrei/lonlat_bng) to do the heavy lifting.

As Julia is JIT compiled, all of the conversions can remain within the same language and still exhibit good performance - this package has been written exclusively in Julia and has only a single dependency.

## Accuracy
Conversions which solely use Helmert transforms are only accurate to approximately 5 metres and are **not suitable** for calculations or conversions where accuracy is important. To correct for these distortions in the OSGB36 terrestrial reference frame, Ordnance Survey has developed a ‘rubber-sheet’ transformation model called OSTN15. This model is available as grids of northing and easting shifts between ETRS89 and OSGB36 covering the UK at a one km resolution. Precision easting and northing shifts for each point can be obtained by a bilinear interpolation. [See here](https://www.ordnancesurvey.co.uk/documents/resources/guide-coordinate-systems-great-britain.pdf) for more information.

This package is capable of converting all forty ETRS89 lon,lat OS test coordinates to OSGB36 eastings, northings with a maximum error of less than 2 mm, and completes the reverse transformation at an accuracy of a single digit at 8 d.p (~1 mm).

## Installation
This package requires Julia v1.0.5 as a minimum and it is recommended that Julia v1.3.1 onwards is used. As this Julia package is currently unregistered, simply enter the following in Pkg mode:
```julia
Pkg> add https://github.com/joshuanunn/ConvertBNG.jl
```

## Multi-threading for high performance
Multi-threading is used automatically when passing arrays of coordinates to the ```convert_bng``` and ```convert_lonlat``` functions. However, by default Julia only utilises a single thread so calculations will be limited to a single CPU core. If high performance is desired, the environment variable ```JULIA_NUM_THREADS``` should be set to the match the number of cores available. Note that this must be set *before* launching the Julia REPL.

For example, enter ```set JULIA_NUM_THREADS=4``` into the Windows command prompt or ```export JULIA_NUM_THREADS=4``` into the Linux terminal. These settings only last the length of the terminal session, so you should launch Julia from the same terminal afterwards. If you want this setting to remain persistent, change the user environment variable in the Windows control panel (```User Accounts``` -> ```Change my environment variables```) or add the above export command to your ```.bashrc``` file or equivalent in Linux. See the **benchmarks** section below for timings.

If you already have the environment variable ```JULIA_NUM_THREADS``` set to > 1, but you only want to use this package in a single threaded mode, simply call the equivalent ```convert_bng_nothread``` or ```convert_lonlat_nothread``` functions with your array of coordinates.

## Use
Following installation, import the package to make all functions available:
```julia
julia> using ConvertBNG
```

The primary functions available in this package are:
* ```convert_bng```	        Perform Lon, Lat to OSGB36 Eastings, Northings conversion, using OSTN15 data
* ```convert_lonlat```  Convert OSGB36 Eastings, Northing coordinates to Lon, Lat using OSTN15 data

Coordinates can be passed to these functions as a single pair or as a 2-column array of any length. Transformed coordinates are always returned as a 2-column array, irrespective of the input format. Note that `lon`, `lat` coordinates outside the [UK bounding box](http://spatialreference.org/ref/epsg/27700/) will be transformed to `[NaN NaN]`, which cannot be mapped.

Can be used to transform single points:
```julia
julia> e,n = convert_bng(-2.594679, 51.450521)
1×2 Array{Float64,2}:
 358772.71  172560.997

julia> e
358772.71

julia> n
172560.997
```

Or an array of points:
```julia
julia> lats_lons = [-5.003508 56.79685; -3.211528 54.454222; -4.076231 53.068497]
3×2 Array{Float64,2}:
 -5.003508  56.79685
 -3.211528  54.454222
 -4.076231  53.068497

julia> convert_bng(lats_lons)
3×2 Array{Float64,2}:
 216676.483  771283.429
 321549.542  507211.985
 260985.952  354375.893
 
julia> convert_lonlat(convert_bng(lats_lons))
3×2 Array{Float64,2}:
 -5.003508    56.79685
 -3.211528    54.454222
 -4.076231    53.068497
```

## In Detail
The above conversions can be broken into their component steps by using the following piecewise functions:
* ```convert_bng``` is equivalent to:
  1. Perform Longitude, Latitude to ETRS89 conversion: ```easting, northing = convert_etrs89(longitude, latitude)```
  2. Perform ETRS89 to OSGB36 conversion, using OSTN15 data: ```easting_corr, northing_corr = convert_etrs89_to_osgb36(easting, northing)```

* ```convert_lonlat``` is equivalent to:
  1. Convert OSGB36 coordinates to ETRS89, using OSTN15 data: ```x, y = convert_osgb36_to_etrs89(easting, northing)```
  2. Convert ETRS89 coordinates to Lon, Lat: ```longitude, latitude = convert_etrs89_to_ll(x, y)```

Using the Caister Water Tower example in the OS Transformations and OSGM15 User Guide - Annexe A, we can describe this location using the following coordinates:
* ETRS89 longitude, latitude: (1.716073973, 52.658007833)
* ETRS89 eastings, northings: (651307.003, 313255.686)
* OSGB36 eastings, northings: (651409.804, 313177.450)

All of these coordinates can be interrelated using the following function calls, which give the corresponding return values:
```julia
julia> convert_bng(1.716073973, 52.658007833)             == [651409.804 313177.450]
julia> convert_etrs89(1.716073973, 52.658007833)          == [651307.003 313255.686]
julia> convert_etrs89_to_osgb36(651307.003, 313255.686)   == [651409.804 313177.450]

julia> convert_lonlat(651409.804, 313177.450)             == [1.71607397 52.65800783]
julia> convert_osgb36_to_etrs89(651409.804, 313177.450)   == [651307.003 313255.686]
julia> convert_etrs89_to_ll(651307.003, 313255.686)       == [1.71607397 52.65800783]
```
The difference between ETRS89 and OSGB36 is determined by extracting the raw OSTN15 shifts for the area of interest and using bilinear interpolation to derive a more accurate set of shifts. For this example the raw shifts are ```[102.787, -78.242, 44.236]```, which are subsequently refined to ```[102.801, -78.236, 44.228]```.

## Benchmarks
Benchmarks were run on two different systems - a high performance Linux system with an 8-core i7-9700K CPU and a mainstream Windows 10 system with a 4-core i5-7500 CPU. The [benchmark code](/benchmark/runbenchmarks.jl) is derived from that used in the [lonlat_bng](https://github.com/urschrei/lonlat_bng) Rust library and involves getting the total time for converting 50 lots of 1 million random lon,lat pairs to BNG and vice versa. The code can be run in the Julia REPL: ```include("/path/to/package/benchmark/runbenchmarks.jl")```.
* The benchmarks were run a number of times to ensure that the total times for each case were consistent and the best time from any run is quoted in the results tables below.
* Julia versions ```v1.0.5```, ```v1.3.1``` and the latest ```v1.4.0-rc2.0``` were tested - all quoted times are those using Julia ```v1.3.1```.
* The quoted times for Julia ```v1.3.1``` were very consistent and similar to those obtained with the latest ```v1.4.0-rc2.0```. However, Julia ```v1.0.5``` times were upto 25% slower and inconsistent - it is recommended to use later versions.
* All Julia benchmarks were run with multi-threading, using 1, 4 or 8 CPU cores for comparison.
* For context, the tables also include the equivalent benchmark times for the [convertbng](https://github.com/urschrei/convertbng) library which uses a [Rust binary](https://github.com/urschrei/lonlat_bng). This has been run in both exposed operational modes (Ctypes and Cython). These binaries are setup to fully utilise the CPU.

Results for **8-Core** i7-9700K CPU @ 3.60GHz

| Function | Threads/Cores | Julia (s) | Rust Ctypes (s) | Rust Cython (s) |
| :---: | :---: | :---: | :---: | :---: |
| convert_bng | 1 | 14.0 | - | - |
| convert_bng | 4 | 4.6 | - | - |
| convert_bng | 8 | 3.1 | 3.3 | 2.0 |
| convert_latlon | 1 | 22.8 | - | - |
| convert_latlon | 4 | 6.8 | - | - |
| convert_latlon | 8 | 4.2 | 6.3 | 3.2 |

Results for **4-Core** i5-7500 CPU @ 3.40GHz

| Function | Threads/Cores | Julia (s) | Rust Ctypes (s) | Rust Cython (s) |
| :---: | :---: | :---: | :---: | :---: |
| convert_bng | 1 | 22.3 | - | - |
| convert_bng | 4 | 8.0 | 8.5 | 5.5 |
| convert_latlon | 1 | 34.5 | - | - |
| convert_latlon | 4 | 10.8 | 16.2 | 9.8 |

## Tests
After installation, the unit tests can (and should) be run by entering the following at the Julia REPL:
```julia
julia> using Pkg; Pkg.test("ConvertBNG")
```
## License
[MIT](license.txt)  

This software is a derived work of the Rust library [lonlat_bng](https://github.com/urschrei/lonlat_bng), released by Stephan Hügel under the MIT license.

This software makes use of OSTN15 data, which is © Copyright and database rights Ordnance Survey Limited 2016. All rights reserved. Provided under the BSD 2-clause [license](OSTN15_license.txt).
