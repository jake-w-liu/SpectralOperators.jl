using LinearAlgebra, Random

# Independent real collocation matrices, including the even-derivative Nyquist.
function channel_matrix(n, L, kind)
    h = 2pi / Float64(L)
    modes = [m <= n ÷ 2 ? m : m - n for m = 0:n-1]
    return [
        sum(modes) do m
            phase = 2pi * m * (i - j) / n
            kind === :lap ? -(h * m)^2 * cos(phase) / n :
            kind === :first ? (iseven(n) && m == n ÷ 2 ? 0.0 : -h * m * sin(phase) / n) :
            (iseven(n) && m == n ÷ 2 ? h * m * cos(phase) / n : 0.0)
        end for i = 0:n-1, j = 0:n-1
    ]
end
function channel_apply(M, a, d)
    result = zeros(Float64, size(a))
    for I in CartesianIndices(a), j in axes(a, d)
        K = CartesianIndex(ntuple(c -> c == d ? j : I[c], ndims(a)))
        result[I] += M[I[d], j] * a[K]
    end
    return result
end
channel_allocated(out, a, g, d) = @allocated nyquist_channel!(out, a, g, d)

@testset "Full-Laplacian Nyquist energy channel" begin
    rng = Xoshiro(0x53a28)
    for T in (Float32, Float64), D = 1:3, variant = 1:3
        n = ntuple(d -> variant == 3 && d == D ? 1 : (variant == 1 ? 6 : 5) + d % 2, D)
        L = ntuple(d -> T(3 + d / 3), D)
        g = FourierGrid(n, L)
        a = randn(rng, T, n)
        out = similar(a)
        first = similar(a)
        lap = similar(a)
        laplacian!(lap, a, g)
        heat = 0.0
        for d = 1:D
            reference = channel_apply(channel_matrix(n[d], L[d], :nyquist), a, d)
            @test nyquist_channel!(out, a, g, d) === out
            scale = max(1.0, maximum(abs, reference))
            @test maximum(abs, out .- reference) <= 512eps(T) * scale
            deriv!(first, a, g, d)
            @test maximum(
                abs,
                first .- channel_apply(channel_matrix(n[d], L[d], :first), a, d),
            ) <= 512eps(T) * max(1.0, maximum(abs, first))
            heat += sum(abs2, Float64.(first)) + sum(abs2, Float64.(out))
            inplace = copy(a)
            nyquist_channel!(inplace, inplace, g, d)
            @test isequal(inplace, out)
            nyquist_channel!(out, a, g, d)
            @test channel_allocated(out, a, g, d) == 0
        end
        loss = -dot(vec(Float64.(a)), vec(Float64.(lap)))
        @test abs(loss - heat) <= 512eps(T) * max(1.0, loss, heat)
        # Full Laplacian is checked independently, not reconstructed from channels.
        reference_lap =
            sum(channel_apply(channel_matrix(n[d], L[d], :lap), a, d) for d = 1:D)
        @test maximum(abs, lap .- reference_lap) <=
              512eps(T) * max(1.0, maximum(abs, reference_lap))
    end
end

@testset "Nyquist channel boundary and range contracts" begin
    for T in (Float32, Float64)
        g = FourierGrid((8, 5), (T(2pi), T(2pi)))
        a = ones(T, g.n)
        out = fill(T(7), g.n)
        for d in (0, 3)
            @test_throws ArgumentError nyquist_channel!(out, a, g, d)
            @test all(==(T(7)), out)
        end
        @test_throws DimensionMismatch nyquist_channel!(zeros(T, 7, 5), a, g, 1)
        @test_throws DimensionMismatch nyquist_channel!(out, zeros(T, 7, 5), g, 1)
        for d = 1:2
            a[1] = T(NaN)
            @test_throws ArgumentError nyquist_channel!(out, a, g, d)
            @test all(==(T(7)), out)
            a[1] = one(T)
        end
        # Every mixed-axis mode is retained; a pure Nyquist is not doubled.
        for I in CartesianIndices(a)
            a[I] = isodd(I[1]) ? T(0.125) : T(-0.125)
        end
        nyquist_channel!(out, a, g, 1)
        @test out == 4 .* a
        nyquist_channel!(out, a, g, 2)
        @test all(iszero, out)
        e = T === Float64 ? 600 : 70
        for (amplitude, k) in
            ((ldexp(one(T), e), T(4)), (ldexp(one(T), -e), ldexp(one(T), e ÷ 2)))
            gr = FourierGrid((8,), (T(8pi) / k,))
            ar = T[isodd(i) ? amplitude : -amplitude for i = 1:8]
            expected = (gr.kfull[1][5] * amplitude)
            result = similar(ar)
            nyquist_channel!(result, ar, gr, 1)
            @test all(isfinite, result)
            @test maximum(abs, result .- (ar ./ amplitude) .* expected) <=
                  8eps(T) * abs(expected)
        end
    end
end
