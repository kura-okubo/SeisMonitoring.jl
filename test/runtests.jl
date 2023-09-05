using Test
using SeisMonitoring

@testset "SeisMonitoring test" begin

include("dvvmeasurement_validation.jl")
include("run_seismonitoring_process.jl")
include("auxiliary_test.jl")

end
