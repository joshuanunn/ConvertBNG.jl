module ConvertBNG

import DelimitedFiles

export convert_osgb36, convert_etrs89

const AIRY_1830_SEMI_MAJOR = 6377563.396
const AIRY_1830_SEMI_MINOR = 6356256.909
const GRS80_SEMI_MAJOR = 6378137.000
const GRS80_SEMI_MINOR = 6356752.3141
const RAD = π / 180.0
# lon and lat of true origin
const LAM0 = RAD * -2.0
const PHI0 = RAD * 49.0
# Easting and Northing of origin
const E0 = 400000.0
const N0 = -100000.0
# convergence factor
const F0 = 0.9996012717
# Coordinate bounds for eastings, northings
const MAX_EASTING = 700000.000
const MAX_NORTHING = 1250000.000
# Coordinate bounds for lon, lat
const MIN_LONGITUDE = -7.5600
const MAX_LONGITUDE = 1.7800
const MIN_LATITUDE = 49.9600
const MAX_LATITUDE = 60.8400
# Const container for ostn15 shift data
const OSTN15_DATA = Ref{Array{Float64,2}}()

"""
Import ostn15 shifts data file into array.
Line position in the file becomes array index and so it is critical that the
raw datafile is not midified.
"""
function __init__()
    # Import ostn15 lookup array
    ostn15_file = joinpath(@__DIR__, "..", "data", "ostn15.dat")
    global OSTN15_DATA[] = DelimitedFiles.readdlm(ostn15_file, ',', Float64)
end

"""
Bounds checking for input values
Return input value if within specified bounds, else return NaN
"""
function check(to_check::T, bounds::Tuple{T,T}) where {T<:Real}
    if bounds[1] <= to_check <= bounds[2]
        return to_check
    end
    return NaN
end

"""
Round a float to nearest mm
"""
function round_to_mm(position_m::T) where {T<:AbstractFloat}
    return round(position_m * 1000.0, digits=0) / 1000.0
end

"""
Round a float to eight decimal places
"""
function round_to_eight(x::T, y::T) where {T<:AbstractFloat}
    new_x = round(x * 100000000.0, digits=0) / 100000000.0
    new_y = round(y * 100000000.0, digits=0) / 100000000.0
    return new_x, new_y
end

"""
Intermediate calculation used for lon, lat to ETRS89 and reverse conversion
"""
function compute_m(phi::T, b::T, n::T) where {T<:Real}
    p_plus = phi + PHI0
    p_minus = phi - PHI0
    result = b * F0 * ((1.0 + n * (1.0 + 5.0 / 4.0 * n * (1.0 + n))) * p_minus -
             3.0 * n * (1.0 + n * (1.0 + 7.0 / 8.0 * n)) * sin(p_minus) * cos(p_plus) +
             (15.0 / 8.0 * n * (n * (1.0 + n))) * sin(2.0 * p_minus) * cos(2.0 * p_plus) -
             35.0 / 24.0 * n ^ 3 * sin(3.0 * p_minus) * cos(3.0 * p_plus))
    return result
end

"""
Perform Longitude, Latitude to ETRS89 conversion
"""
function convert_etrs89(longitude::T, latitude::T) where {T<:Real}
    # Input is restricted to the UK bounding box
    # Convert bounds-checked input to degrees, or return an Err
    lambd = RAD * check(longitude, (MIN_LONGITUDE, MAX_LONGITUDE))
    phi = RAD * check(latitude, (MIN_LATITUDE, MAX_LATITUDE))
    # ellipsoid squared eccentricity constant
    e2 = (GRS80_SEMI_MAJOR ^ 2 - GRS80_SEMI_MINOR ^ 2) / GRS80_SEMI_MAJOR ^ 2
    n = (GRS80_SEMI_MAJOR - GRS80_SEMI_MINOR) / (GRS80_SEMI_MAJOR + GRS80_SEMI_MINOR)
    
    sp2 = sin(phi) ^ 2
    nu = GRS80_SEMI_MAJOR * F0 * (1.0 - e2 * sp2) ^ -0.5 # v
    rho = GRS80_SEMI_MAJOR * F0 * (1.0 - e2) * (1.0 - e2 * sp2) ^ -1.5
    eta2 = nu / rho - 1.0
    
    m = compute_m(phi, GRS80_SEMI_MINOR, n)
    
    cp = cos(phi)
    sp = sin(phi)
    tp = tan(phi)
    tp2 = tp ^ 2
    tp4 = tp ^ 4
    
    I = m + N0
    II = nu / 2.0 * sp * cp
    III = nu / 24.0 * sp * cp ^ 3 * (5.0 - tp2 + 9.0 * eta2)
    @fastmath IIIA = nu / 720.0 * sp * cp ^ 5 * (61.0 - 58.0 * tp2 + tp4)
    IV = nu * cp
    V = nu / 6.0 * cp ^ 3 * (nu / rho - tp2)
    @fastmath VI = nu / 120.0 * cp ^ 5 * (5.0 - 18.0 * tp2 + tp4 + 14.0 * eta2 - 58.0 * tp2 * eta2)
    l = lambd - LAM0
    @fastmath north = I + II * l ^ 2 + III * l ^ 4 + IIIA * l ^ 6
    @fastmath east = E0 + IV * l + V * l ^ 3 + VI * l ^ 5

    return round_to_mm(east), round_to_mm(north)
end

"""
Easting and Northing to Lon, Lat conversion using a Helmert transform
Note that either GRS80 or Airy 1830 ellipsoids can be passed
"""
@fastmath function convert_to_ll(eastings::T, northings::T, ell_a::T, ell_b::T) where {T<:Real}
    # ellipsoid squared eccentricity constant
    a = ell_a
    b = ell_b
    e2 = (a ^ 2 - b ^ 2) / a ^ 2
    n = (a - b) / (a + b)

    dN = check(northings, (0.0, MAX_NORTHING)) - N0
    phi = PHI0 + dN / (a * F0)
    m = compute_m(phi, b, n)
    while (dN - m) >= 0.001
        phi += (dN - m) / (a * F0)
        m = compute_m(phi, b, n)
    end
    sp2 = sin(phi) ^ 2
    nu = a * F0 * (1.0 - e2 * sp2) ^ -0.5
    rho = a * F0 * (1.0 - e2) * (1.0 - e2 * sp2) ^ -1.5
    eta2 = nu / rho - 1.0

    tp = tan(phi)
    tp2 = tp ^ 2
    tp4 = tp ^ 4

    VII = tp / (2.0 * rho * nu)
    VIII = tp / (24.0 * rho * nu ^ 3) * (5.0 + 3.0 * tp2 + eta2 - 9.0 * tp2 * eta2)
    IX = tp / (720.0 * rho * nu ^ 5) * (61.0 + 90.0 * tp2 + 45.0 * tp4)

    sp = 1.0 / cos(phi)
    tp6 = tp4 * tp2

    X = sp / nu
    XI = sp / (6.0 * nu ^ 3) * (nu / rho + 2.0 * tp2)
    XII = sp / (120.0 * nu ^ 5) * (5.0 + 28.0 * tp2 + 24.0 * tp4)
    XIIA = sp / (5040.0 * nu ^ 7) * (61.0 + 662.0 * tp2 + 1320.0 * tp4 + 720.0 * tp6)

    e = check(eastings, (0.000, MAX_EASTING)) - E0

    phi = phi - VII * e ^ 2 + VIII * e ^ 4 - IX * e ^ 6
    lambd = LAM0 + X * e - XI * e ^ 3 + XII * e ^ 5 - XIIA * e ^ 7

    phi = rad2deg(phi)
    lambd = rad2deg(lambd)
    return round_to_eight(lambd), round_to_eight(phi)
end

function convert_etrs89_to_ll(E::T, N::T) where {T<:Real}
    # ETRS89 uses the WGS84 / GRS80 ellipsoid constants
    return convert_to_ll(E, N, GRS80_SEMI_MAJOR, GRS80_SEMI_MINOR)
end

"""
Perform Longitude, Latitude to OSGB36 conversion, using [OSTN15] data
"""
function convert_osgb36(longitude::T, latitude::T) where {T<:Real}
    # Convert input to ETRS89
    eastings, northings = convert_etrs89(longitude, latitude)
    # Obtain OSTN15 corrections, and incorporate
    e_shift, n_shift, _ = ostn15_shifts(eastings, northings)
    e_corr = round_to_mm(eastings + e_shift)
    n_corr = round_to_mm(northings + n_shift)
    return e_corr, n_corr
end

function convert_osgb36(lonlat::Array{T,2}) where{T<:Real}
    # Create empty output array
    en_arr = zeros(T, size(lonlat))
    # Transform
    #Threads.@threads
    for i in 1:size(en_arr)[1]
        en = convert_osgb36(lonlat[i,1], lonlat[i,2])
        en_arr[i,1] = round(en[1], digits=3)
        en_arr[i,2] = round(en[2], digits=3)
    end
    return en_arr
end

"""
Perform ETRS89 to OSGB36 conversion, using [OSTN15] data
"""
function convert_etrs89_to_osgb36(eastings::T, northings::T) where {T<:Real}
    # Obtain OSTN15 corrections, and incorporate
    e_shift, n_shift, _ = ostn15_shifts(eastings, northings)
    return round_to_mm(eastings + e_shift), round_to_mm(northings + n_shift)
end

"""
Convert OSGB36 coordinates to Lon, Lat using OSTN15 data
"""
function convert_osgb36_to_ll(E::T, N::T) where {T<:Real}
    # Apply reverse OSTN15 adustments
    epsilon = 0.009
    dx, dy, _ = ostn15_shifts(E, N)
    x, y = E - dx, N - dy
    last_dx, last_dy = dx, dy
    while true
        dx, dy = ostn15_shifts(x, y)
        x = E - dx
        y = N - dy
        # If the difference […] is more than 0.00010m (User Guide, p15)
        # TODO: invert this logic
        if abs(dx - last_dx) < epsilon && abs(dy - last_dy) < epsilon
            break
        end
        last_dx = dx
        last_dy = dy
    end
    x = round_to_mm(E - dx)
    y = round_to_mm(N - dy)
    # We've converted to ETRS89, so we need to use the WGS84/ GRS80 ellipsoid constants
    return convert_to_ll(x, y, GRS80_SEMI_MAJOR, GRS80_SEMI_MINOR)
end

"""
Convert OSGB36 coordinates to ETRS89 using OSTN15 data
"""
function convert_osgb36_to_etrs89(E::T, N::T) where {T<:Real}
    # Apply reverse OSTN15 adustments
    epsilon = 0.00001
    dx, dy, _ = ostn15_shifts(E, N)
    x, y = E - dx, N - dy
    last_dx, last_dy = dx, dy
    while true
        dx, dy = ostn15_shifts(x, y)
        x = E - dx
        y = N - dy
        if abs(dx - last_dx) < epsilon && abs(dy - last_dy) < epsilon
            break
        end
        last_dx = dx
        last_dy = dy
    end
    x = round_to_mm(E - dx)
    y = round_to_mm(N - dy)
    return x, y
end

"""
Try to get OSTN15 shift parameters, and calculate offsets
"""
function get_ostn_ref(x::T, y::T) where {T<:Integer}
    # Calculate unique index key
    index = x + (y * 701) + 1
    # Some or None, so convert to Result, which we can try!
    se = OSTN15_DATA[][index,1]
    sn = OSTN15_DATA[][index,2]
    sg = OSTN15_DATA[][index,3]
    return se, sn, sg
end

"""
Input values must be valid ETRS89 grid references
Calculate OSTN15 shifts for a given coordinate
"""
function ostn15_shifts(x::T, y::T) where {T<:Real}
    # Ensure that we're within the boundaries
    if check(x, (0.0, MAX_EASTING)) === NaN || check(y, (0.0, MAX_NORTHING)) === NaN
        return NaN, NaN, NaN
    end
    # Calculate base ostn15 lookup indices (type must be integer)
    e_index = convert(Int, round(x / 1000.0, digits=0))
    n_index = convert(Int, round(y / 1000.0, digits=0))
    
    # eastings and northings of the south-west corner of the cell
    x0 = e_index * 1000
    y0 = n_index * 1000

    # The easting, northing and geoid shifts for the four corners of the cell
    # intersections
    # this is a 3 x 4 matrix (using column-major order)
    # bottom-left grid intersection
    s0_1, s0_2, s0_3 = get_ostn_ref(e_index, n_index)
    # bottom-right
    s1_1, s1_2, s1_3 = get_ostn_ref(e_index + 1, n_index)
    # top-left
    s2_1, s2_2, s2_3 = get_ostn_ref(e_index, n_index + 1)
    # top-right
    s3_1, s3_2, s3_3 = get_ostn_ref(e_index + 1, n_index + 1)
    
    # offset within square
    dx = x - x0
    dy = y - y0

    t = dx / 1000.0
    u = dy / 1000.0

    # Calculation of the weights for each intersection (W)
    # this is a 4 x 1 matrix
    f0 = (1.0 - t) * (1.0 - u)
    f1 = t * (1.0 - u)
    f2 = (1.0 - t) * u
    f3 = t * u
    
    # Bilinear interpolation, to obtain the actual shifts
    # We could also do a dot product: weights.dot(offsets)
    se = f0 * s0_1 + f1 * s1_1 + f2 * s2_1 + f3 * s3_1
    sn = f0 * s0_2 + f1 * s1_2 + f2 * s2_2 + f3 * s3_2
    # this isn't needed for this library, since it's a height offset
    sg = f0 * s0_3 + f1 * s1_3 + f2 * s2_3 + f3 * s3_3
    
    return round_to_mm(se), round_to_mm(sn), round_to_mm(sg)
end

end # module
