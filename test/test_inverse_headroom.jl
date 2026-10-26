# The filter reference is an independent separable convolution. It uses neither
# an FFT nor SpectralOperators metadata/multiplier helpers.
function inverse_cutoff_response(n)
    kernel = setprecision(BigFloat, 256) do
        cutoff = (n - 1) ÷ 3
        [Float64(sum(cospi(BigFloat(2m) * j / n) for m = -cutoff:cutoff) / n) for j = 0:n-1]
    end
    signs = sign.(kernel)
    response = [sum(kernel[mod(x - j, n)+1] * signs[j+1] for j = 0:n-1) for x = 0:n-1]
    return signs, response
end

@testset "Inverse FFT range follows complete operation bound" begin
    for T in (Float32, Float64), n in ((32, 32, 32), (32, 24, 16))
        axes_reference = map(inverse_cutoff_response, n)
        signs = map(first, axes_reference)
        response = map(last, axes_reference)
        amplitude = floatmax(T) / T(100)
        retained = prod(2 * ((size - 1) ÷ 3) + 1 for size in n)
        # Exact Parseval/Cauchy bound certifies every mathematical output finite.
        @test Rational{BigInt}(amplitude)^2 * retained < Rational{BigInt}(floatmax(T))^2
        grid = FourierGrid(n, (T(2), T(3), T(5)))
        for shift in (0, 5)
            input = [
                amplitude * T(prod(signs[d][mod(I[d] - 1 - shift, n[d])+1] for d = 1:3)) for
                I in CartesianIndices(n)
            ]
            expected = [
                prod(response[d][mod(I[d] - 1 - shift, n[d])+1] for d = 1:3) for
                I in CartesianIndices(n)
            ]
            @test dealias_two_thirds!(input, grid) === input
            @test all(isfinite, input)
            error = maximum(abs, Float64.(input ./ amplitude) .- expected)
            @test error <= 128eps(T) * max(maximum(abs, expected), 1.0)
        end
    end
end

@testset "Shared scalar headroom includes the complete inverse" begin
    for T in (Float32, Float64),
        count in (1, 32^3, typemax(Int)),
        gain in (zero(T), one(T), T(1e20))

        safe = SpectralOperators._fft_safe_magnitude(T, count, gain)
        @test isfinite(safe) && safe >= zero(T)
        # Independent high-precision norm bound, with half of the native
        # factor-eight margin reserved for rounding of count/sqrt/products.
        bound =
            BigFloat(safe) *
            BigFloat(count) *
            sqrt(BigFloat(count)) *
            max(BigFloat(1), BigFloat(gain))
        @test bound <= BigFloat(floatmax(T)) / 4
    end
end
