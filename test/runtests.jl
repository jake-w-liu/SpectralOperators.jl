using Test
using Aqua
using SpectralOperators

@testset "SpectralOperators" begin
    include("test_operators.jl")
    include("test_spectral_extra.jl")
    include("test_shared_fourier.jl")
    include("test_inverse_headroom.jl")
    include("test_nyquist_channel.jl")
    include("test_scaled_sources.jl")
    include("test_scaled_floor.jl")
    include("test_scaled_zero.jl")
    @testset "Aqua package quality" begin
        Aqua.test_all(SpectralOperators)
    end
end
