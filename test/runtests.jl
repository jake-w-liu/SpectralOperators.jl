using Test

@testset "SpectralOperators" begin
    include("test_operators.jl")
    include("test_spectral_extra.jl")
end
