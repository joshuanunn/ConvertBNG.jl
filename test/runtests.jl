using ConvertBNG
using Test

const PREC_EN = 0.002       #  2 mm for easting, northing absolute accuracy
const PREC_LL = 0.0000001   #  8 dp for lat, lon absolute accuracy
const PREC_FL = 0.000000001 # 10 dp for floating point comparison

const AIRY_1830_SEMI_MAJOR = 6377563.396
const AIRY_1830_SEMI_MINOR = 6356256.909

approx(x, y, atol) = all(((xx, yy),) -> isapprox(xx, yy, atol=atol), zip(x, y))
approx_ll(x, y) = approx(x, y, PREC_LL)
approx_en(x, y) = approx(x, y, PREC_EN)
approx_fl(x, y) = approx(x, y, PREC_FL)

@testset "All tests" begin

    @testset "caister water tower" begin
        # From Transformations and OSGM15 User Guide - Annexe A
        # Caister Water Tower example
        # ETRS89 longitude, latitude: (1.716073973, 52.658007833)
        # ETRS89 eastings, northings: (651307.003, 313255.686)
        # Shifts eastings, northings: (102.801, -78.236)
        # OSGB36 eastings, northings: (651409.804, 313177.450)    
        @test approx_en(convert_etrs89_to_osgb36(651307.003, 313255.686), (651409.804, 313177.450))
        @test approx_en(convert_bng(1.716073973, 52.658007833), (651409.804, 313177.450))
        @test approx_en(convert_etrs89(1.716073973, 52.658007833), (651307.003, 313255.686))
        @test approx_ll(convert_etrs89_to_ll(651307.003, 313255.686), (1.716073973, 52.658007833))
        @test approx_ll(convert_lonlat(651409.804, 313177.450), (1.716073973, 52.658007833))
        # Check OSTN15 shift retrival
        @test approx_fl(get_ostn_ref(651, 313), (102.787, -78.242, 44.236))
        @test approx_fl(ostn15_shifts(651307.003, 313255.686), (102.800790098, -78.235984076, 44.228405219))
        
        # From Transformations and OSGM15 User Guide - Annexe C
        # Caister Water Tower example
        # E,N to lon,lat using Airy 1830 ellipsoid
        @test approx_ll(convert_to_ll(651409.903, 313177.270, AIRY_1830_SEMI_MAJOR, AIRY_1830_SEMI_MINOR), (1.717921583, 52.65757031))
    end

    @testset "min/max bounds tests" begin
        # failed ostn15 array retrieval
        @test get_ostn_ref(999, 999)[1] === NaN
        @test get_ostn_ref(999, 999)[2] === NaN
        @test get_ostn_ref(999, 999)[3] === NaN

        # Test min/max easting, northing, longitude, latitude
        # bad max easting
        @test convert_lonlat(700000.1, 1240000.0)[1] === NaN
        @test convert_lonlat(700000.1, 1240000.0)[2] === NaN
        # bad max northing
        @test convert_lonlat(600000.0, 1250000.1)[1] === NaN
        @test convert_lonlat(600000.0, 1250000.1)[2] === NaN
        # bad lon
        @test convert_bng(181., 51.44533267)[1] === NaN
        @test convert_bng(181., 51.44533267)[2] === NaN
        # bad lat
        @test convert_bng(-0.32824866, -90.01)[1] === NaN
        @test convert_bng(-0.32824866, -90.01)[2] === NaN

        # Test extent bounds of check function
        # below min lon
        max_lon = 1.768960
        min_lon = -6.379880
        @test check(-6.379881, (min_lon, max_lon)) === NaN
        # below min lat
        max_lat = 55.811741
        min_lat = 49.871159
        @test check(49.871158, (min_lat, max_lat)) === NaN
        # above max lon
        max_lon = 1.768960
        min_lon = -6.379880
        @test check(1.768961, (min_lon, max_lon)) === NaN
        # above max lat
        max_lat = 55.811741
        min_lat = 49.871159
        @test check(55.811742, (min_lat, max_lat)) === NaN
    end

    @testset "OS latlon to OSGB36 testset" begin
        # TP01 to TP40
        @test approx_en(convert_bng(-6.299777520, 49.92226394),  (  91492.146,   11318.804))
        @test approx_en(convert_bng(-5.203046100, 49.96006138),  ( 170370.718,   11572.405))
        @test approx_en(convert_bng(-4.108645636, 50.43885826),  ( 250359.811,   62016.569))
        @test approx_en(convert_bng(-1.297822772, 50.57563665),  ( 449816.371,   75335.861))
        @test approx_en(convert_bng(-1.450514337, 50.93127938),  ( 438710.920,  114792.250))
        @test approx_en(convert_bng(-3.551283492, 51.40078220),  ( 292184.870,  168003.465))
        @test approx_en(convert_bng( 1.444547304, 51.37447026),  ( 639821.835,  169565.858))
        @test approx_en(convert_bng(-2.544076183, 51.42754743),  ( 362269.991,  169978.690))
        @test approx_en(convert_bng(-0.119925572, 51.48936565),  ( 530624.974,  178388.464))
        @test approx_en(convert_bng(-4.308524770, 51.85890896),  ( 241124.584,  220332.641))
        @test approx_en(convert_bng( 0.897243270, 51.89436637),  ( 599445.590,  225722.826))
        @test approx_en(convert_bng(-2.154586144, 52.25529382),  ( 389544.190,  261912.153))
        @test approx_en(convert_bng(-0.912489570, 52.25160951),  ( 474335.969,  262047.755))
        @test approx_en(convert_bng( 0.401535471, 52.75136687),  ( 562180.547,  319784.995))
        @test approx_en(convert_bng(-1.197476559, 52.96219109),  ( 454002.834,  340834.943))
        @test approx_en(convert_bng(-2.640493208, 53.34480280),  ( 357455.843,  383290.436))
        @test approx_en(convert_bng(-4.289180698, 53.41628516),  ( 247958.971,  393492.909))
        @test approx_en(convert_bng(-4.289177929, 53.41630925),  ( 247959.241,  393495.583))
        @test approx_en(convert_bng(-3.040454907, 53.77911026),  ( 331534.564,  431920.794))
        @test approx_en(convert_bng(-1.663791682, 53.80021520),  ( 422242.186,  433818.701))
        @test approx_en(convert_bng(-4.634521682, 54.08666318),  ( 227778.330,  468847.388))
        @test approx_en(convert_bng(-0.077731332, 54.11685144),  ( 525745.670,  470703.214))
        @test approx_en(convert_bng(-4.388491181, 54.32919541),  ( 244780.636,  495254.887))
        @test approx_en(convert_bng(-2.938277411, 54.89542340),  ( 339921.145,  556034.761))
        @test approx_en(convert_bng(-1.616576852, 54.97912274),  ( 424639.355,  565012.703))
        @test approx_en(convert_bng(-4.296490163, 55.85399953),  ( 256340.925,  664697.269))
        @test approx_en(convert_bng(-3.294792193, 55.92478266),  ( 319188.434,  670947.534))
        @test approx_en(convert_bng(-5.828366919, 57.00606696),  ( 167634.202,  797067.144))
        @test approx_en(convert_bng(-2.048560307, 57.13902519),  ( 397160.491,  805349.736))
        @test approx_en(convert_bng(-4.219263986, 57.48625001),  ( 267056.768,  846176.972))
        @test approx_en(convert_bng(-8.578544561, 57.81351838),  (   9587.909,  899448.996))
        @test approx_en(convert_bng(-7.592555606, 58.21262247),  (  71713.132,  938516.404))
        @test approx_en(convert_bng(-6.260914555, 58.51560361),  ( 151968.652,  966483.780))
        @test approx_en(convert_bng(-3.726310221, 58.58120461),  ( 299721.891,  967202.992))
        @test approx_en(convert_bng(-3.214540011, 59.03743871),  ( 330398.323, 1017347.016))
        @test approx_en(convert_bng(-4.417576746, 59.09335035),  ( 261596.778, 1025447.602))
        @test approx_en(convert_bng(-5.827993398, 59.09671617),  ( 180862.461, 1029604.114))
        @test approx_en(convert_bng(-1.625169661, 59.53470794),  ( 421300.525, 1072147.239))
        @test approx_en(convert_bng(-1.274869104, 59.85409914),  ( 440725.073, 1107878.448))
        @test approx_en(convert_bng(-2.073828228, 60.13308092),  ( 395999.668, 1138728.951))
    end

    @testset "OS OSGB36 to latlon testset" begin
        # TP01 to TP40
        @test approx_ll(convert_lonlat(  91492.146,   11318.804), (-6.299777520, 49.92226394))
        @test approx_ll(convert_lonlat( 170370.718,   11572.405), (-5.203046100, 49.96006138))
        @test approx_ll(convert_lonlat( 250359.811,   62016.569), (-4.108645636, 50.43885826))
        @test approx_ll(convert_lonlat( 449816.371,   75335.861), (-1.297822772, 50.57563665))
        @test approx_ll(convert_lonlat( 438710.920,  114792.250), (-1.450514337, 50.93127938))
        @test approx_ll(convert_lonlat( 292184.870,  168003.465), (-3.551283492, 51.40078220))
        @test approx_ll(convert_lonlat( 639821.835,  169565.858), ( 1.444547304, 51.37447026))
        @test approx_ll(convert_lonlat( 362269.991,  169978.690), (-2.544076183, 51.42754743))
        @test approx_ll(convert_lonlat( 530624.974,  178388.464), (-0.119925572, 51.48936565))
        @test approx_ll(convert_lonlat( 241124.584,  220332.641), (-4.308524770, 51.85890896))
        @test approx_ll(convert_lonlat( 599445.590,  225722.826), ( 0.897243270, 51.89436637))
        @test approx_ll(convert_lonlat( 389544.190,  261912.153), (-2.154586144, 52.25529382))
        @test approx_ll(convert_lonlat( 474335.969,  262047.755), (-0.912489570, 52.25160951))
        @test approx_ll(convert_lonlat( 562180.547,  319784.995), ( 0.401535471, 52.75136687))
        @test approx_ll(convert_lonlat( 454002.834,  340834.943), (-1.197476559, 52.96219109))
        @test approx_ll(convert_lonlat( 357455.843,  383290.436), (-2.640493208, 53.34480280))
        @test approx_ll(convert_lonlat( 247958.971,  393492.909), (-4.289180698, 53.41628516))
        @test approx_ll(convert_lonlat( 247959.241,  393495.583), (-4.289177929, 53.41630925))
        @test approx_ll(convert_lonlat( 331534.564,  431920.794), (-3.040454907, 53.77911026))
        @test approx_ll(convert_lonlat( 422242.186,  433818.701), (-1.663791682, 53.80021520))
        @test approx_ll(convert_lonlat( 227778.330,  468847.388), (-4.634521682, 54.08666318))
        @test approx_ll(convert_lonlat( 525745.670,  470703.214), (-0.077731332, 54.11685144))
        @test approx_ll(convert_lonlat( 244780.636,  495254.887), (-4.388491181, 54.32919541))
        @test approx_ll(convert_lonlat( 339921.145,  556034.761), (-2.938277411, 54.89542340))
        @test approx_ll(convert_lonlat( 424639.355,  565012.703), (-1.616576852, 54.97912274))
        @test approx_ll(convert_lonlat( 256340.925,  664697.269), (-4.296490163, 55.85399953))
        @test approx_ll(convert_lonlat( 319188.434,  670947.534), (-3.294792193, 55.92478266))
        @test approx_ll(convert_lonlat( 167634.202,  797067.144), (-5.828366919, 57.00606696))
        @test approx_ll(convert_lonlat( 397160.491,  805349.736), (-2.048560307, 57.13902519))
        @test approx_ll(convert_lonlat( 267056.768,  846176.972), (-4.219263986, 57.48625001))
        @test approx_ll(convert_lonlat(   9587.909,  899448.996), (-8.578544561, 57.81351838))
        @test approx_ll(convert_lonlat(  71713.132,  938516.404), (-7.592555606, 58.21262247))
        @test approx_ll(convert_lonlat( 151968.652,  966483.780), (-6.260914555, 58.51560361))
        @test approx_ll(convert_lonlat( 299721.891,  967202.992), (-3.726310221, 58.58120461))
        @test approx_ll(convert_lonlat( 330398.323, 1017347.016), (-3.214540011, 59.03743871))
        @test approx_ll(convert_lonlat( 261596.778, 1025447.602), (-4.417576746, 59.09335035))
        @test approx_ll(convert_lonlat( 180862.461, 1029604.114), (-5.827993398, 59.09671617))
        @test approx_ll(convert_lonlat( 421300.525, 1072147.239), (-1.625169661, 59.53470794))
        @test approx_ll(convert_lonlat( 440725.073, 1107878.448), (-1.274869104, 59.85409914))
        @test approx_ll(convert_lonlat( 395999.668, 1138728.951), (-2.073828228, 60.1330809281))
    end

    @testset "Array input-output" begin
        # Check results from single- and multiple-point methods are the same
        let n = 100, lons = 5.0.*(rand(n) .- 0.5), lats = 55 .+ 5.0.*rand(n)
            en = vcat([hcat(convert_bng(lon, lat)...) for (lon, lat) in zip(lons, lats)]...)
            # Broadcast identity comparison to check NaNs are the same
            @test all(convert_bng(hcat(lons, lats)) .=== en)
            ll = vcat([hcat(convert_lonlat(east, north)...) for (east, north) in eachrow(en)]...)
            @test all(convert_lonlat(en) .=== ll)
        end
    end

end
