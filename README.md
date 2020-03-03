## Introduction
A Julia library for fast conversion between WGS84 longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates. Conversions use a standard 7-element Helmert transform with the addition of OSTN15 corrections for [accuracy](#accuracy).

This is a derived work of original repositories written by Stephan Hügel, which include a Rust binary [lonlat_bng](https://github.com/urschrei/lonlat_bng) and the corresponding Python utility library [convertbng](https://github.com/urschrei/convertbng). As Python is relatively slow, [convertbng](https://github.com/urschrei/convertbng) utilises a [Rust binary](https://github.com/urschrei/lonlat_bng) to do the heavy lifting.

As Julia is JIT compiled, all of the conversions can remain within the same lanuguage and still exhibit good performance - this package has been written exclusively in Julia and has only a single external dependency.

## Accuracy
Conversions which solely use Helmert transforms are only accurate to approximately 5 metres and are **not suitable** for calculations or conversions where accuracy is important. To correct for these distortions in the OSGB36 terrestrial reference frame, Ordnance Survey has developed a ‘rubber-sheet’ transformation model called OSTN15. This model is available as grids of northing and easting shifts between ETRS89 and OSGB36 covering the UK at a one km resolution. Precision easting and northing shifts for each point can be obtained by a bilinear interpolation. [See here](https://www.ordnancesurvey.co.uk/documents/resources/guide-coordinate-systems-great-britain.pdf) for more information.

This package is capable of converting all forty ETRS89 lon,lat OS test coordinates to OSGB36 eastings, northings with a maximum error of less than 2 mm, and completes the reverse transformation at an accuracy of a single digit at 8 d.p (~1 mm).


## Installation
As this Julia package is currently unregistered, simply enter the following in Pkg mode:
```julia
Pkg> add https://github.com/joshuanunn/ConvertBNG.jl
```

## Use
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
Multi-threading support to be added soon for Julia users of v1.3 onwards. Detailed performance comparisons will be added at this point.

## Tests
After installation, the unit tests can (and should) be run by entering the following at the Julia REPL:
```julia
julia> using Pkg; Pkg.test("ConvertBNG")
```
## License
[MIT](license.txt)  

This software is a derived work of the Rust library [lonlat_bng](https://github.com/urschrei/lonlat_bng), released by Stephan Hügel under the MIT license.

This software makes use of OSTN15 data, which is © Copyright and database rights Ordnance Survey Limited 2016. All rights reserved. Provided under the BSD 2-clause [license](OSTN15_license.txt).
