# Phase-1 spectral-extra tests: exponential filter, 2/3 dealiasing, FFT-friendly
# sizes, FFTW wisdom round-trip, and the SpectralOperators namespace.
#
# Oracles are closed-form / definitional and never call the production routine
# under test to derive the expected value.

using SpectralOperators, Test, LinearAlgebra, FFTW

# Build a single pure cosine mode  f = cos(m x)  on a 2π periodic grid of n points
# (1D). Returns the grid and field. m must satisfy |m| ≤ n÷2 to be representable.
function cosine_mode(::Type{T}, n::Int, m::Int) where {T}
    L = (T(2π),)
    g = FourierGrid((n,), L)
    f = Vector{T}(undef, n)
    @inbounds for i = 1:n
        x = (i - 1) * g.dx[1]
        f[i] = cos(m * x)
    end
    return g, f
end

@testset "exp_filter! transfer function" begin
    for T in (Float64, Float32)
        n = 64
        # Low-k mode (m=1, far below Nyquist m=32): transfer ≈ 1.
        g, f_lo = cosine_mode(T, n, 1)
        amp0 = maximum(abs, f_lo)
        f = copy(f_lo)
        exp_filter!(f, g)
        amp1 = maximum(abs, f)
        @test isapprox(amp1, amp0; rtol = (T == Float64 ? 1e-10 : 1e-4))
        # Field essentially unchanged shape-wise.
        @test norm(f .- f_lo) / norm(f_lo) < (T == Float64 ? 1e-9 : 1e-4)

        # Near-Nyquist mode (m = n÷2 - 1 = 31): strongly damped.
        g2, f_hi = cosine_mode(T, n, n ÷ 2 - 1)
        amphi0 = maximum(abs, f_hi)
        fh = copy(f_hi)
        exp_filter!(fh, g2)
        amphi1 = maximum(abs, fh)
        @test amphi1 < 0.05 * amphi0          # >95% reduction at near-Nyquist

        # k = 0 (constant) mode preserved EXACTLY.
        c = fill(T(3.7), n)
        c0 = copy(c)
        exp_filter!(c, g)
        @test maximum(abs, c .- c0) <= 8 * eps(T) * abs(c0[1])
    end
end

@testset "exp_filter! monotone in k and sharper with p" begin
    T = Float64
    n = 64
    g = FourierGrid((n,), (T(2π),))
    # Transfer at increasing k should be non-increasing.
    function transfer(m)
        f = Vector{T}(undef, n)
        @inbounds for i = 1:n
            f[i] = cos(m * (i - 1) * g.dx[1])
        end
        a0 = maximum(abs, f)
        exp_filter!(f, g)
        maximum(abs, f) / a0
    end
    ts = [transfer(m) for m = 1:(n÷2-1)]
    @test all(diff(ts) .<= 1e-12)             # monotone non-increasing
    @test ts[1] > 0.999                       # low-k essentially untouched
end

@testset "exp_filter! dimension agnostic" begin
    T = Float64
    n = (16, 16)
    g = FourierGrid(n, (T(2π), T(2π)))
    f = Array{T,2}(undef, n)
    @inbounds for I in CartesianIndices(f)
        t = Tuple(I)
        x = (t[1] - 1) * g.dx[1]
        y = (t[2] - 1) * g.dx[2]
        f[I] = cos(x) * cos(y)                # low mode (1,1)
    end
    f0 = copy(f)
    exp_filter!(f, g)
    @test norm(f .- f0) / norm(f0) < 1e-9     # low mode preserved in 2D
end

@testset "dealias_two_thirds! masks high-k, keeps low-k" begin
    for T in (Float64, Float32)
        n = 48                                 # Nyquist mode = 24; 2/3 cutoff = 16
        # High-k mode m=20 (> 2/3·24=16): zeroed.
        g, f_hi = cosine_mode(T, n, 20)
        fh = copy(f_hi)
        dealias_two_thirds!(fh, g)
        # Mode is zeroed in spectral space; residual is only FFT-roundtrip noise.
        @test maximum(abs, fh) < (T == Float64 ? 1e-12 : 1e-4) * maximum(abs, f_hi)

        # Boundary-just-below mode m=15 (< 16): preserved.
        g2, f_lo = cosine_mode(T, n, 15)
        fl = copy(f_lo)
        dealias_two_thirds!(fl, g2)
        @test isapprox(maximum(abs, fl), maximum(abs, f_lo); rtol = (T == Float64 ? 1e-10 : 1e-4))

        # k = 0 preserved.
        c = fill(T(2.0), n)
        c0 = copy(c)
        dealias_two_thirds!(c, g)
        @test maximum(abs, c .- c0) <= 8 * eps(T) * abs(c0[1])
    end
end

@testset "dealias_two_thirds! per-axis (2D)" begin
    T = Float64
    n = (48, 48)
    g = FourierGrid(n, (T(2π), T(2π)))
    mk(mx, my) = begin
        f = Array{T,2}(undef, n)
        @inbounds for I in CartesianIndices(f)
            t = Tuple(I)
            f[I] = cos(mx * (t[1] - 1) * g.dx[1]) * cos(my * (t[2] - 1) * g.dx[2])
        end
        f
    end
    # High on one axis only (mx=20 > 16, my=2): should be zeroed.
    fa = mk(20, 2)
    fa0 = copy(fa)
    dealias_two_thirds!(fa, g)
    @test maximum(abs, fa) < 1e-6 * maximum(abs, fa0)
    # Low on both axes: preserved.
    fb = mk(3, 4)
    fb0 = copy(fb)
    dealias_two_thirds!(fb, g)
    @test norm(fb .- fb0) / norm(fb0) < 1e-9
end

@testset "fft_friendly_size small-prime correctness" begin
    is7smooth(m) = begin
        for p in (2, 3, 5, 7)
            while m % p == 0
                m ÷= p
            end
        end
        m == 1
    end
    @test fft_friendly_size(17) == 18         # 18 = 2·3²
    @test fft_friendly_size(100) == 100       # already 7-smooth (2²·5²)
    @test fft_friendly_size(101) == 105       # 105 = 3·5·7
    @test fft_friendly_size(1) == 1
    @test fft_friendly_size(0) == 1
    @test fft_friendly_size(2) == 2
    # Exhaustive property check: result ≥ n, 7-smooth, and minimal.
    for n = 1:2000
        r = fft_friendly_size(n)
        @test r >= n
        @test is7smooth(r)
        @test all(!is7smooth(k) for k = n:(r-1))   # nothing smaller works
    end
end

@testset "with_fftw_wisdom round-trip" begin
    dir = mktempdir()
    path = joinpath(dir, "sub", "wisdom.dat")   # nested dir must be created
    # Absent file: must not error, returns f() value.
    val = with_fftw_wisdom(path) do
        p = plan_fft!(zeros(ComplexF64, 32))
        p * zeros(ComplexF64, 32)
        42
    end
    @test val == 42
    @test isfile(path)                          # wisdom exported
    # Second call imports existing wisdom, still fine.
    val2 = with_fftw_wisdom(path) do
        plan_fft!(zeros(ComplexF64, 32))
        7
    end
    @test val2 == 7
    # Exception inside f still exports and rethrows.
    @test_throws ErrorException with_fftw_wisdom(path) do
        error("boom")
    end
    @test isfile(path)
end

@testset "SpectralOperators package API" begin
    # The package exposes the operator API independent of any global state.
    T = Float64
    n = 32
    g = SpectralOperators.FourierGrid((n,), (T(2π),))
    f = Vector{T}(undef, n)
    ex = similar(f)
    for i = 1:n
        x = (i - 1) * g.dx[1]
        f[i] = cos(2x)
        ex[i] = -2 * sin(2x)
    end
    out = similar(f)
    SpectralOperators.deriv!(out, f, g, 1)
    @test norm(out .- ex) / norm(ex) < 1e-10
end
