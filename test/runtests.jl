using Test
using Aqua
using SpectralOperators

@testset "SpectralOperators" begin
    include("test_operators.jl")
    include("test_spectral_extra.jl")
    include("test_shared_fourier.jl")
    include("test_inverse_headroom.jl")
    @testset "Aqua package quality" begin
        Aqua.test_all(SpectralOperators)
    end
end
