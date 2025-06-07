using Aqua


@testset "Aqua.jl checks" begin
    Aqua.test_all(IPIcalculator)
end
