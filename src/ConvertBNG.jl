module ConvertBNG

import DelimitedFiles

export convert_osgb36
export convert_etrs89
export convert_osgb36_to_ll
export convert_to_ll
export get_ostn_ref
export ostn15_shifts
export check
export convert_etrs89_to_ll
export convert_etrs89_to_osgb36

const AIRY_1830_SEMI_MAJOR = 6377563.396
const AIRY_1830_SEMI_MINOR = 6356256.909
const GRS80_SEMI_MAJOR = 6378137.000
const GRS80_SEMI_MINOR = 6356752.3141
const RAD = π / 180.0
# lon and lat of true origin
const λ₀ = RAD * -2.0
const ϕ₀ = RAD * 49.0
# Easting and Northing of origin
const E₀ = 400000.0
const N₀ = -100000.0
# convergence factor
const F₀ = 0.9996012717
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
function check(to_check::T1, bounds::Tuple{T2,T2}) where {T1<:Real, T2<:Real}
    if bounds[1] <= to_check <= bounds[2]
        return to_check
    else
        return NaN
    end
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
function compute_m(ϕ::T, b::T, n::T) where {T<:Real}
    n² = n ^ 2
    n³ = n ^ 3
    result = b * F₀ * (
        (1.0 + n + 5.0 / 4.0 * n² + 5.0 / 4.0 * n³) * (ϕ - ϕ₀) -
        (3.0 * n + 3.0 * n² + 21.0 / 8.0 * n³) * sin(ϕ - ϕ₀) * cos(ϕ + ϕ₀) +
        (15.0 / 8.0 * n² + 15.0 / 8.0 * n³) * sin(2.0 * (ϕ - ϕ₀)) * cos(2.0 * (ϕ + ϕ₀)) -
        35.0 / 24.0 * n³ * sin(3.0 * (ϕ - ϕ₀)) * cos(3.0 * (ϕ + ϕ₀))
    )
    return result
end

"""
Perform Longitude, Latitude to ETRS89 conversion
"""
function convert_etrs89(longitude::T, latitude::T) where {T<:Real}
    # Input is restricted to the UK bounding box
    # Convert bounds-checked input to degrees, or return an Err
    λ = deg2rad(check(longitude, (MIN_LONGITUDE, MAX_LONGITUDE)))
    ϕ = deg2rad(check(latitude, (MIN_LATITUDE, MAX_LATITUDE)))
    # ellipsoid squared eccentricity constant
    e² = (GRS80_SEMI_MAJOR ^ 2 - GRS80_SEMI_MINOR ^ 2) / GRS80_SEMI_MAJOR ^ 2
    n = (GRS80_SEMI_MAJOR - GRS80_SEMI_MINOR) / (GRS80_SEMI_MAJOR + GRS80_SEMI_MINOR)
    
    sin²ϕ = sin(ϕ) ^ 2
    ν = GRS80_SEMI_MAJOR * F₀ * (1.0 - e² * sin²ϕ) ^ -0.5
    ρ = GRS80_SEMI_MAJOR * F₀ * (1.0 - e²) * (1.0 - e² * sin²ϕ) ^ -1.5
    η² = ν / ρ - 1.0
    
    m = compute_m(ϕ, GRS80_SEMI_MINOR, n)
    
    cosϕ = cos(ϕ)
    cos³ϕ = cosϕ ^ 3
    cos⁵ϕ = cosϕ ^ 5
    sinϕ = sin(ϕ)
    tanϕ = tan(ϕ)
    tan²ϕ = tanϕ ^ 2
    tan⁴ϕ = tanϕ ^ 4
    
    I = m + N₀
    II = ν / 2.0 * sinϕ * cosϕ
    III = ν / 24.0 * sinϕ * cos³ϕ * (5.0 - tan²ϕ + 9.0 * η²)
    @fastmath IIIA = ν / 720.0 * sinϕ * cos⁵ϕ * (61.0 - 58.0 * tan²ϕ + tan⁴ϕ)
    IV = ν * cosϕ
    V = ν / 6.0 * cos³ϕ * (ν / ρ - tan²ϕ)
    @fastmath VI = ν / 120.0 * cos⁵ϕ * (5.0 - 18.0 * tan²ϕ + tan⁴ϕ + 14.0 * η² - 58.0 * tan²ϕ * η²)
    l = λ - λ₀
    @fastmath north = I + II * l ^ 2 + III * l ^ 4 + IIIA * l ^ 6
    @fastmath east = E₀ + IV * l + V * l ^ 3 + VI * l ^ 5

    return [round_to_mm(east) round_to_mm(north)]
end

"""
Easting and Northing to Lon, Lat conversion using a Helmert transform
Note that either GRS80 or Airy 1830 ellipsoids can be passed
"""
@fastmath function convert_to_ll(eastings::T, northings::T, ℓ_a::T, ℓ_b::T) where {T<:Real}
    # ellipsoid squared eccentricity constant
    a = ℓ_a
    b = ℓ_b
    e² = (a ^ 2 - b ^ 2) / a ^ 2
    n = (a - b) / (a + b)

    dN = check(northings, (0.0, MAX_NORTHING)) - N₀
    ϕ = ϕ₀ + dN / (a * F₀)
    m = compute_m(ϕ, b, n)
    while (dN - m) >= 0.001
        ϕ += (dN - m) / (a * F₀)
        m = compute_m(ϕ, b, n)
    end
    
    sin²ϕ = sin(ϕ) ^ 2
    tanϕ = tan(ϕ)
    tan²ϕ = tanϕ ^ 2
    tan⁴ϕ = tanϕ ^ 4
    tan⁶ϕ = tan⁴ϕ * tan²ϕ    
    secϕ = 1.0 / cos(ϕ)
    
    ν = a * F₀ * (1.0 - e² * sin²ϕ) ^ -0.5
    ρ = a * F₀ * (1.0 - e²) * (1.0 - e² * sin²ϕ) ^ -1.5
    η² = ν / ρ - 1.0

    VII = tanϕ / (2.0 * ρ * ν)
    VIII = tanϕ / (24.0 * ρ * ν ^ 3) * (5.0 + 3.0 * tan²ϕ + η² - 9.0 * tan²ϕ * η²)
    IX = tanϕ / (720.0 * ρ * ν ^ 5) * (61.0 + 90.0 * tan²ϕ + 45.0 * tan⁴ϕ)
    X = secϕ / ν
    XI = secϕ / (6.0 * ν ^ 3) * (ν / ρ + 2.0 * tan²ϕ)
    XII = secϕ / (120.0 * ν ^ 5) * (5.0 + 28.0 * tan²ϕ + 24.0 * tan⁴ϕ)
    XIIA = secϕ / (5040.0 * ν ^ 7) * (61.0 + 662.0 * tan²ϕ + 1320.0 * tan⁴ϕ + 720.0 * tan⁶ϕ)

    e = check(eastings, (0.000, MAX_EASTING)) - E₀

    ϕ = ϕ - VII * e ^ 2 + VIII * e ^ 4 - IX * e ^ 6
    λ = λ₀ + X * e - XI * e ^ 3 + XII * e ^ 5 - XIIA * e ^ 7

    ϕ = rad2deg(ϕ)
    λ = rad2deg(λ)
    λ, ϕ = round_to_eight(λ, ϕ)
    return [λ ϕ]
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
    return [e_corr n_corr]
end

function convert_osgb36(lonlat::Array{T,2}) where{T<:Real}
    # Create empty output array
    en_arr = zeros(T, size(lonlat))
    # Transform
    #Threads.@threads
    for i in 1:size(en_arr)[1]
        en_arr[i,:] .= convert_osgb36(lonlat[i,1], lonlat[i,2])
    end
    return en_arr
end

"""
Perform ETRS89 to OSGB36 conversion, using [OSTN15] data
"""
function convert_etrs89_to_osgb36(eastings::T, northings::T) where {T<:Real}
    # Obtain OSTN15 corrections, and incorporate
    e_shift, n_shift, _ = ostn15_shifts(eastings, northings)
    e_corr = round_to_mm(eastings + e_shift)
    n_corr = round_to_mm(northings + n_shift)
    return [e_corr n_corr]
end

"""
Convert OSGB36 coordinates to Lon, Lat using OSTN15 data
"""
function convert_osgb36_to_ll(E::T, N::T) where {T<:Real}
    # Apply reverse OSTN15 adustments
    epsilon = 0.0001
    dx, dy, _ = ostn15_shifts(E, N)
    x, y = E - dx, N - dy
    last_dx, last_dy = dx, dy
    while true
        dx, dy = ostn15_shifts(x, y)
        x = E - dx
        y = N - dy
        # If the difference […] is more than 0.00010m (User Guide, p15)
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

function convert_osgb36_to_ll(en::Array{T,2}) where{T<:Real}
    # Create empty output array
    ll_arr = zeros(T, size(en))
    # Transform
    #Threads.@threads
    for i in 1:size(ll_arr)[1]
        ll_arr[i,:] .= convert_osgb36_to_ll(en[i,1], en[i,2])
    end
    return ll_arr
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
    return [x y]
end

"""
Try to get OSTN15 shift parameters, and calculate offsets
"""
function get_ostn_ref(x::T, y::T) where {T<:Integer}
    # Check inputs sensible
    if (0 <= x <= 700) && (0 <= y <= 1250)
        # Calculate unique index key
        index = x + (y * 701) + 1
        # Some or None, so convert to Result, which we can try!
        se = OSTN15_DATA[][index,1]
        sn = OSTN15_DATA[][index,2]
        sg = OSTN15_DATA[][index,3]
    else
        se = NaN
        sn = NaN
        sg = NaN
    end
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
    e_index = trunc(Int, x / 1000.0)
    n_index = trunc(Int, y / 1000.0)
    
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
