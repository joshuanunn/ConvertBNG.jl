using ConvertBNG
using Test

@test convert_osgb36(1.716073973, 52.658007833)[1] ≈ 651409.804 atol=0.01
@test convert_osgb36(1.716073973, 52.658007833)[2] ≈ 313177.45 atol=0.01
