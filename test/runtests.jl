using ConvertBNG
using Test

const PREC_EN = 0.002       #  2 mm for easting, northing absolute accuracy
const PREC_LL = 0.0000001   #  8 dp for lat, lon absolute accuracy
const PREC_FL = 0.000000001 # 10 dp for floating point comparison

const AIRY_1830_SEMI_MAJOR = 6377563.396
const AIRY_1830_SEMI_MINOR = 6356256.909

@testset "All tests" begin

    @testset "caister water tower" begin
        # From Transformations and OSGM15 User Guide - Annexe A
        # Caister Water Tower example
        # ETRS89 longitude, latitude: (1.716073973, 52.658007833)
        # ETRS89 eastings, northings: (651307.003, 313255.686)
        # Shifts eastings, northings: (102.801, -78.236)
        # OSGB36 eastings, northings: (651409.804, 313177.450)    
        @test convert_etrs89_to_osgb36(651307.003, 313255.686) ≈ [651409.804 313177.450] atol=PREC_EN
        @test convert_osgb36(1.716073973, 52.658007833) ≈ [651409.804 313177.450] atol=PREC_EN
        @test convert_etrs89(1.716073973, 52.658007833) ≈ [651307.003 313255.686] atol=PREC_EN
        @test convert_etrs89_to_ll(651307.003, 313255.686) ≈ [1.716073973 52.658007833] atol=PREC_LL
        @test convert_osgb36_to_ll(651409.804, 313177.450) ≈ [1.716073973 52.658007833] atol=PREC_LL
        # Check OSTN15 shift retrival
        @test collect(get_ostn_ref(651, 313)) ≈ [102.787, -78.242, 44.236] atol=PREC_FL
        @test collect(ostn15_shifts(651307.003, 313255.686)) ≈ [102.801, -78.236, 44.228] atol=PREC_FL
        
        # From Transformations and OSGM15 User Guide - Annexe C
        # Caister Water Tower example
        # E,N to lon,lat using Airy 1830 ellipsoid
        @test convert_to_ll(651409.903, 313177.270, AIRY_1830_SEMI_MAJOR, AIRY_1830_SEMI_MINOR) ≈ [1.717921583 52.65757031] atol=PREC_LL
    end

    @testset "min/max bounds tests" begin
        # failed ostn15 array retrieval
        @test get_ostn_ref(999, 999)[1] === NaN
        @test get_ostn_ref(999, 999)[2] === NaN
        @test get_ostn_ref(999, 999)[3] === NaN

        # Test min/max easting, northing, longitude, latitude
        # bad max easting
        @test convert_etrs89_to_osgb36(700000.1, 1240000.0)[1] === NaN
        @test convert_etrs89_to_osgb36(700000.1, 1240000.0)[2] === NaN
        # bad max northing
        @test convert_etrs89_to_osgb36(600000.0, 1250000.1)[1] === NaN
        @test convert_etrs89_to_osgb36(600000.0, 1250000.1)[2] === NaN
        # bad lon
        @test convert_osgb36(181., 51.44533267)[1] === NaN
        @test convert_osgb36(181., 51.44533267)[2] === NaN
        # bad lat
        @test convert_osgb36(-0.32824866, -90.01)[1] === NaN
        @test convert_osgb36(-0.32824866, -90.01)[2] === NaN

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
        @test convert_osgb36(-6.299777520, 49.92226394) ≈ [  91492.146   11318.804] atol=PREC_EN
        @test convert_osgb36(-5.203046100, 49.96006138) ≈ [ 170370.718   11572.405] atol=PREC_EN
        @test convert_osgb36(-4.108645636, 50.43885826) ≈ [ 250359.811   62016.569] atol=PREC_EN
        @test convert_osgb36(-1.297822772, 50.57563665) ≈ [ 449816.371   75335.861] atol=PREC_EN
        @test convert_osgb36(-1.450514337, 50.93127938) ≈ [ 438710.920  114792.250] atol=PREC_EN
        @test convert_osgb36(-3.551283492, 51.40078220) ≈ [ 292184.870  168003.465] atol=PREC_EN
        @test convert_osgb36( 1.444547304, 51.37447026) ≈ [ 639821.835  169565.858] atol=PREC_EN
        @test convert_osgb36(-2.544076183, 51.42754743) ≈ [ 362269.991  169978.690] atol=PREC_EN
        @test convert_osgb36(-0.119925572, 51.48936565) ≈ [ 530624.974  178388.464] atol=PREC_EN
        @test convert_osgb36(-4.308524770, 51.85890896) ≈ [ 241124.584  220332.641] atol=PREC_EN
        @test convert_osgb36( 0.897243270, 51.89436637) ≈ [ 599445.590  225722.826] atol=PREC_EN
        @test convert_osgb36(-2.154586144, 52.25529382) ≈ [ 389544.190  261912.153] atol=PREC_EN
        @test convert_osgb36(-0.912489570, 52.25160951) ≈ [ 474335.969  262047.755] atol=PREC_EN
        @test convert_osgb36( 0.401535471, 52.75136687) ≈ [ 562180.547  319784.995] atol=PREC_EN
        @test convert_osgb36(-1.197476559, 52.96219109) ≈ [ 454002.834  340834.943] atol=PREC_EN
        @test convert_osgb36(-2.640493208, 53.34480280) ≈ [ 357455.843  383290.436] atol=PREC_EN
        @test convert_osgb36(-4.289180698, 53.41628516) ≈ [ 247958.971  393492.909] atol=PREC_EN
        @test convert_osgb36(-4.289177929, 53.41630925) ≈ [ 247959.241  393495.583] atol=PREC_EN
        @test convert_osgb36(-3.040454907, 53.77911026) ≈ [ 331534.564  431920.794] atol=PREC_EN
        @test convert_osgb36(-1.663791682, 53.80021520) ≈ [ 422242.186  433818.701] atol=PREC_EN
        @test convert_osgb36(-4.634521682, 54.08666318) ≈ [ 227778.330  468847.388] atol=PREC_EN
        @test convert_osgb36(-0.077731332, 54.11685144) ≈ [ 525745.670  470703.214] atol=PREC_EN
        @test convert_osgb36(-4.388491181, 54.32919541) ≈ [ 244780.636  495254.887] atol=PREC_EN
        @test convert_osgb36(-2.938277411, 54.89542340) ≈ [ 339921.145  556034.761] atol=PREC_EN
        @test convert_osgb36(-1.616576852, 54.97912274) ≈ [ 424639.355  565012.703] atol=PREC_EN
        @test convert_osgb36(-4.296490163, 55.85399953) ≈ [ 256340.925  664697.269] atol=PREC_EN
        @test convert_osgb36(-3.294792193, 55.92478266) ≈ [ 319188.434  670947.534] atol=PREC_EN
        @test convert_osgb36(-5.828366919, 57.00606696) ≈ [ 167634.202  797067.144] atol=PREC_EN
        @test convert_osgb36(-2.048560307, 57.13902519) ≈ [ 397160.491  805349.736] atol=PREC_EN
        @test convert_osgb36(-4.219263986, 57.48625001) ≈ [ 267056.768  846176.972] atol=PREC_EN
        @test convert_osgb36(-8.578544561, 57.81351838) ≈ [   9587.909  899448.996] atol=PREC_EN
        @test convert_osgb36(-7.592555606, 58.21262247) ≈ [  71713.132  938516.404] atol=PREC_EN
        @test convert_osgb36(-6.260914555, 58.51560361) ≈ [ 151968.652  966483.780] atol=PREC_EN
        @test convert_osgb36(-3.726310221, 58.58120461) ≈ [ 299721.891  967202.992] atol=PREC_EN
        @test convert_osgb36(-3.214540011, 59.03743871) ≈ [ 330398.323 1017347.016] atol=PREC_EN
        @test convert_osgb36(-4.417576746, 59.09335035) ≈ [ 261596.778 1025447.602] atol=PREC_EN
        @test convert_osgb36(-5.827993398, 59.09671617) ≈ [ 180862.461 1029604.114] atol=PREC_EN
        @test convert_osgb36(-1.625169661, 59.53470794) ≈ [ 421300.525 1072147.239] atol=PREC_EN
        @test convert_osgb36(-1.274869104, 59.85409914) ≈ [ 440725.073 1107878.448] atol=PREC_EN
        @test convert_osgb36(-2.073828228, 60.13308092) ≈ [ 395999.668 1138728.951] atol=PREC_EN
    end

    @testset "OS OSGB36 to latlon testset" begin
        # TP01 to TP40
        @test convert_osgb36_to_ll(  91492.146,   11318.804) ≈ [-6.299777520 49.92226394] atol=PREC_LL
        @test convert_osgb36_to_ll( 170370.718,   11572.405) ≈ [-5.203046100 49.96006138] atol=PREC_LL
        @test convert_osgb36_to_ll( 250359.811,   62016.569) ≈ [-4.108645636 50.43885826] atol=PREC_LL
        @test convert_osgb36_to_ll( 449816.371,   75335.861) ≈ [-1.297822772 50.57563665] atol=PREC_LL
        @test convert_osgb36_to_ll( 438710.920,  114792.250) ≈ [-1.450514337 50.93127938] atol=PREC_LL
        @test convert_osgb36_to_ll( 292184.870,  168003.465) ≈ [-3.551283492 51.40078220] atol=PREC_LL
        @test convert_osgb36_to_ll( 639821.835,  169565.858) ≈ [ 1.444547304 51.37447026] atol=PREC_LL
        @test convert_osgb36_to_ll( 362269.991,  169978.690) ≈ [-2.544076183 51.42754743] atol=PREC_LL
        @test convert_osgb36_to_ll( 530624.974,  178388.464) ≈ [-0.119925572 51.48936565] atol=PREC_LL
        @test convert_osgb36_to_ll( 241124.584,  220332.641) ≈ [-4.308524770 51.85890896] atol=PREC_LL
        @test convert_osgb36_to_ll( 599445.590,  225722.826) ≈ [ 0.897243270 51.89436637] atol=PREC_LL
        @test convert_osgb36_to_ll( 389544.190,  261912.153) ≈ [-2.154586144 52.25529382] atol=PREC_LL
        @test convert_osgb36_to_ll( 474335.969,  262047.755) ≈ [-0.912489570 52.25160951] atol=PREC_LL
        @test convert_osgb36_to_ll( 562180.547,  319784.995) ≈ [ 0.401535471 52.75136687] atol=PREC_LL
        @test convert_osgb36_to_ll( 454002.834,  340834.943) ≈ [-1.197476559 52.96219109] atol=PREC_LL
        @test convert_osgb36_to_ll( 357455.843,  383290.436) ≈ [-2.640493208 53.34480280] atol=PREC_LL
        @test convert_osgb36_to_ll( 247958.971,  393492.909) ≈ [-4.289180698 53.41628516] atol=PREC_LL
        @test convert_osgb36_to_ll( 247959.241,  393495.583) ≈ [-4.289177929 53.41630925] atol=PREC_LL
        @test convert_osgb36_to_ll( 331534.564,  431920.794) ≈ [-3.040454907 53.77911026] atol=PREC_LL
        @test convert_osgb36_to_ll( 422242.186,  433818.701) ≈ [-1.663791682 53.80021520] atol=PREC_LL
        @test convert_osgb36_to_ll( 227778.330,  468847.388) ≈ [-4.634521682 54.08666318] atol=PREC_LL
        @test convert_osgb36_to_ll( 525745.670,  470703.214) ≈ [-0.077731332 54.11685144] atol=PREC_LL
        @test convert_osgb36_to_ll( 244780.636,  495254.887) ≈ [-4.388491181 54.32919541] atol=PREC_LL
        @test convert_osgb36_to_ll( 339921.145,  556034.761) ≈ [-2.938277411 54.89542340] atol=PREC_LL
        @test convert_osgb36_to_ll( 424639.355,  565012.703) ≈ [-1.616576852 54.97912274] atol=PREC_LL
        @test convert_osgb36_to_ll( 256340.925,  664697.269) ≈ [-4.296490163 55.85399953] atol=PREC_LL
        @test convert_osgb36_to_ll( 319188.434,  670947.534) ≈ [-3.294792193 55.92478266] atol=PREC_LL
        @test convert_osgb36_to_ll( 167634.202,  797067.144) ≈ [-5.828366919 57.00606696] atol=PREC_LL
        @test convert_osgb36_to_ll( 397160.491,  805349.736) ≈ [-2.048560307 57.13902519] atol=PREC_LL
        @test convert_osgb36_to_ll( 267056.768,  846176.972) ≈ [-4.219263986 57.48625001] atol=PREC_LL
        @test convert_osgb36_to_ll(   9587.909,  899448.996) ≈ [-8.578544561 57.81351838] atol=PREC_LL
        @test convert_osgb36_to_ll(  71713.132,  938516.404) ≈ [-7.592555606 58.21262247] atol=PREC_LL
        @test convert_osgb36_to_ll( 151968.652,  966483.780) ≈ [-6.260914555 58.51560361] atol=PREC_LL
        @test convert_osgb36_to_ll( 299721.891,  967202.992) ≈ [-3.726310221 58.58120461] atol=PREC_LL
        @test convert_osgb36_to_ll( 330398.323, 1017347.016) ≈ [-3.214540011 59.03743871] atol=PREC_LL
        @test convert_osgb36_to_ll( 261596.778, 1025447.602) ≈ [-4.417576746 59.09335035] atol=PREC_LL
        @test convert_osgb36_to_ll( 180862.461, 1029604.114) ≈ [-5.827993398 59.09671617] atol=PREC_LL
        @test convert_osgb36_to_ll( 421300.525, 1072147.239) ≈ [-1.625169661 59.53470794] atol=PREC_LL
        @test convert_osgb36_to_ll( 440725.073, 1107878.448) ≈ [-1.274869104 59.85409914] atol=PREC_LL
        @test convert_osgb36_to_ll( 395999.668, 1138728.951) ≈ [-2.073828228 60.13308092] atol=PREC_LL
    end

end