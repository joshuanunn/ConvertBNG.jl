
# Introduction
A Julia library for fast conversion between WGS84 longitude and latitude and British National Grid ([epsg:27700](http://spatialreference.org/ref/epsg/osgb-1936-british-national-grid/)) coordinates. Conversions use a standard 7-element Helmert transform with the addition of OSTN15 corrections for [accuracy](#accuracy).

This is a derived work of original repositories written by Stephan Hügel, which include a Rust binary [lonlat_bng](https://github.com/urschrei/lonlat_bng) and the corresponding Python utility library [convertbng](https://github.com/urschrei/convertbng).

The moviation for [convertbng](https://github.com/urschrei/convertbng) was to offer fast and accurate coordinate conversion utilising OSTN15 to Python users. As Python is relatively slow, this library utilised a Rust binary [lonlat_bng](https://github.com/urschrei/lonlat_bng) to do the main conversions. As Julia is JIT compiled, all of the conversions can remain within the same lanuguage - this package has been written exclusively in Julia and has minimal external dependancies.

# Accuracy
Conversions which solely use Helmert transforms are only accurate to within around 5 metres, and are **not suitable** for calculations or conversions used in e.g. surveying. Thus, we use the OSTN15 transform, which adjusts for local variation within the Terrestrial Reference Frame by incorporating OSTN15 data. [See here](http://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/surveying.html) for more information.

This package sucessfully converts all 40 ETRS89 lon,lat OS test coordinates to OSGB36 eastings, northings with a maximum error of < 2 mm, and complete the reverse transformation at an accuracy of a single digit at 8 d.p (~1 mm).

## Use
As this Julia package is currently unregistered, installation should be as simple as:
```julia
julia> using Pkg; Pkg.add("https://github.com/joshuanunn/ConvertBNG.jl")
```

**Note that `lon`, `lat` coordinates outside the [UK bounding box](http://spatialreference.org/ref/epsg/27700/) will be transformed to `(NaN, NaN)`, which cannot be mapped.**

## Functions
The primary functions available in this package are:
* ```convert_bng```	        Perform Lon, Lat to OSGB36 Eastings, Northings conversion, using OSTN15 data
* ```convert_lonlat```  Convert OSGB36 Eastings, Northing coordinates to Lon, Lat using OSTN15 data

If required, the above conversions can be broken into their component steps by using the following piecewise functions:
* ```convert_bng``` is equivalent to:
  1. Perform Longitude, Latitude to ETRS89 conversion: ```easting, northing = convert_etrs89(longitude, latitude)```
  2. Perform ETRS89 to OSGB36 conversion, using OSTN15 data: ```easting_corr, northing_corr = convert_etrs89_to_osgb36(easting, northing)```

* ```convert_lonlat``` is equivalent to:
  1. Convert OSGB36 coordinates to ETRS89, using OSTN15 data: ```x, y = convert_osgb36_to_etrs89(easting, northing)```
  2. Convert ETRS89 coordinates to Lon, Lat: ```longitude, latitude = convert_etrs89_to_ll(x, y)```

The following uses the Caister Water Tower example taken from OS Transformations and OSGM15 User Guide - Annexe A. This location can be described using the following coordinates:
* ETRS89 longitude, latitude: (1.716073973, 52.658007833)
* ETRS89 eastings, northings: (651307.003, 313255.686)
* OSGB36 eastings, northings: (651409.804, 313177.450)    
All of these coordinates can be interelated using the following function calls, which give the corresponding return values:
```julia
julia> convert_bng(1.716073973, 52.658007833)             == [651409.804 313177.450]
julia> convert_etrs89(1.716073973, 52.658007833)          == [651307.003 313255.686]
julia> convert_etrs89_to_osgb36(651307.003, 313255.686)   == [651409.804 313177.450]

julia> convert_lonlat(651409.804, 313177.450)             == [1.716073973 52.658007833]
julia> convert_etrs89_to_ll(651307.003, 313255.686)       == [1.716073973 52.658007833]
```
The differences between ETRS89 eastings, northings and OSGB36 eastings, northings are described by extracting the raw xyz OSTN15 shifts for the area of interest, and using bilinear interpolation to derive a more accurate set of shifts at the specific location of interest. For this example the raw shifts are ```[102.787, -78.242, 44.236]```, which are subsequently refined to ```[102.801, -78.236, 44.228]```.

## Benchmark
TBC

## Tests
The unittests can be run by entering the following at the Julia REPL:
```julia
julia> using Pkg; Pkg.test("ConvertBNG")
```
## License
[MIT](license.txt)  

This software makes use of OSTN15 data, which is © Copyright and database rights Ordnance Survey Limited 2016. All rights reserved. Provided under the BSD 2-clause [license](OSTN15_license.txt).
