# Operator benchmarks OP-001/002/003/005/006 (periodic Fourier).
# OP-004 (mixed FD/Fourier) is deferred to the mixed-operator phase.
# Oracles are closed-form and never call the production operator under test.

using SpectralOperators, Test, LinearAlgebra, Random, Statistics

relL2(a, b) = norm(a .- b) / norm(b)

struct OffsetSbpVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    data::V
    first_index::Int
end

Base.size(v::OffsetSbpVector) = size(v.data)
Base.axes(v::OffsetSbpVector) = (v.first_index:(v.first_index + length(v.data) - 1),)
Base.IndexStyle(::Type{<:OffsetSbpVector}) = IndexLinear()
Base.getindex(v::OffsetSbpVector, i::Int) = v.data[i - v.first_index + 1]
Base.setindex!(v::OffsetSbpVector, x, i::Int) = (v.data[i - v.first_index + 1] = x)

struct OffsetSbpMatrix{T,V<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::V
    first_i::Int
    first_j::Int
end

Base.size(A::OffsetSbpMatrix) = size(A.data)
Base.axes(A::OffsetSbpMatrix) = (
    A.first_i:(A.first_i + size(A.data, 1) - 1),
    A.first_j:(A.first_j + size(A.data, 2) - 1),
)
Base.IndexStyle(::Type{<:OffsetSbpMatrix}) = IndexCartesian()
Base.getindex(A::OffsetSbpMatrix, i::Int, j::Int) = A.data[i - A.first_i + 1, j - A.first_j + 1]
Base.setindex!(A::OffsetSbpMatrix, x, i::Int, j::Int) = (A.data[i - A.first_i + 1, j - A.first_j + 1] = x)

# Build f = sin(kx x) cos(ky y) cos(kz z) and exact ∂x f on a 2π^D grid (integer modes).
function sincos_field(::Type{T}, n::NTuple{D,Int}, modes::NTuple{D,Int}) where {T,D}
    L = ntuple(_ -> T(2π), D)
    g = FourierGrid(n, L)
    f = Array{T,D}(undef, n)
    dfx = similar(f)
    @inbounds for I in CartesianIndices(f)
        t = Tuple(I)
        x = (t[1] - 1) * g.dx[1]
        val = sin(modes[1] * x)
        dval = modes[1] * cos(modes[1] * x)
        if D >= 2
            y = (t[2] - 1) * g.dx[2]
            val *= cos(modes[2] * y)
            dval *= cos(modes[2] * y)
        end
        if D >= 3
            z = (t[3] - 1) * g.dx[3]
            val *= cos(modes[3] * z)
            dval *= cos(modes[3] * z)
        end
        f[I] = val
        dfx[I] = dval
    end
    return g, f, dfx
end

@testset "OP-001 periodic scalar derivative" begin
    for (T, tol) in ((Float64, 1e-11), (Float32, 1e-5))
        for (n, modes) in (((16,), (3,)), ((16, 16), (3, 2)), ((12, 16, 8), (2, 3, 2)))
            g, f, dfx = sincos_field(T, n, modes)
            out = similar(f)
            deriv!(out, f, g, 1)
            e = relL2(out, dfx)
            @test e < tol
        end
    end
end

@testset "OP-002 div(curl)=0" begin
    Random.seed!(1)
    for T in (Float64, Float32)
        tol = T == Float64 ? 1e-10 : 1e-3
        for n in ((16,), (16, 16), (8, 12, 16))
            D = length(n)
            L = ntuple(_ -> T(2π), D)
            g = FourierGrid(n, L)
            A = ntuple(_ -> randn(T, n...), 3)
            B = ntuple(_ -> similar(first(A)), 3)
            curl!(B, A, g)
            divB = similar(B[1])
            divergence!(divB, B, g)
            kmax = maximum(maximum(abs, g.kvec[d]) for d = 1:D)
            resid = norm(divB) / (kmax * norm(B[1]) + eps(T))
            @test resid < tol
        end
    end
end

@testset "OP-003 divergence-free projection" begin
    Random.seed!(2)
    for T in (Float64, Float32)
        tol = T == Float64 ? 1e-10 : 1e-3
        for n in ((16,), (16, 16), (8, 12, 16))
            D = length(n)
            L = ntuple(_ -> T(2π), D)
            g = FourierGrid(n, L)
            # transverse part (already divergence-free) + uniform background
            A = ntuple(_ -> randn(T, n...), 3)
            Btrans = ntuple(_ -> similar(first(A)), 3)
            curl!(Btrans, A, g)
            bg = (T(0.7), T(-0.3), T(1.1))
            target = ntuple(c -> Btrans[c] .+ bg[c], 3)
            # add a longitudinal (gradient) part that projection must remove
            φ = randn(T, n...)
            gradφ = ntuple(_ -> similar(φ), D)
            gradient!(gradφ, φ, g)
            B = ntuple(c -> copy(target[c]), 3)
            for c = 1:D
                B[c] .+= gradφ[c]
            end
            project_divfree!(B, g)
            # longitudinal removed → recover transverse + background
            err = norm(B[1] .- target[1]) + norm(B[2] .- target[2]) + norm(B[3] .- target[3])
            scale = norm(target[1]) + norm(target[2]) + norm(target[3])
            @test err / scale < tol
            # result is divergence-free
            divB = similar(B[1])
            divergence!(divB, B, g)
            kmax = maximum(maximum(abs, g.kvec[d]) for d = 1:D)
            @test norm(divB) / (kmax * scale + eps(T)) < tol
            # k=0 mean preserved on every component
            for c = 1:3
                @test isapprox(mean(B[c]), bg[c]; atol = (T == Float64 ? 1e-10 : 1e-3))
            end
        end
    end
end

@testset "OP-005 discrete integration by parts" begin
    Random.seed!(3)
    for T in (Float64, Float32)
        tol = T == Float64 ? 1e-10 : 1e-4
        for n in ((16,), (16, 16), (8, 12, 16))
            D = length(n)
            L = ntuple(_ -> T(2π), D)
            g = FourierGrid(n, L)
            f = randn(T, n...)
            h = randn(T, n...)
            for j = 1:D
                lhs = dot(f, deriv(h, g, j))
                rhs = -dot(deriv(f, g, j), h)
                @test abs(lhs - rhs) / (abs(lhs) + abs(rhs) + eps(T)) < tol
            end
        end
    end
end

@testset "OP-006 Nyquist handling" begin
    for T in (Float64, Float32)
        n = (16,)
        g = FourierGrid(n, (T(2π),))
        # pure Nyquist mode cos(N/2 · 2π x / L) = (-1)^i : unrepresentable derivative → 0
        f = T[(-1.0)^(i - 1) for i = 1:n[1]]
        out = deriv!(similar(f), f, g, 1)
        @test maximum(abs, out) < (T == Float64 ? 1e-10 : 1e-4)
        # real input → real output (deriv! returns a real array by construction)
        @test eltype(out) == T
        # a resolved mode is still differentiated correctly
        g2, fr, dfx = sincos_field(T, n, (3,))
        @test relL2(deriv(fr, g2, 1), dfx) < (T == Float64 ? 1e-11 : 1e-5)
    end
end

@testset "mixed SBP/Fourier operator" begin
    T = Float64
    n = 17
    Lx = T(2)
    s = SBP1D(n, Lx)
    x = range(zero(T), Lx; length = n)

    f = collect(sin.(π .* x ./ Lx))
    exact = collect((π / Lx) .* cos.(π .* x ./ Lx))
    out = similar(f)
    sbp_deriv!(out, f, s)
    @test norm(out[2:end-1] .- exact[2:end-1]) / norm(exact[2:end-1]) < 2e-2

    u = collect(cos.(2π .* x ./ Lx))
    v = collect(sin.(π .* x ./ Lx))
    Du = sbp_deriv(u, s)
    Dv = sbp_deriv(v, s)
    lhs = dot(u .* s.H, Dv) + dot(Du .* s.H, v)
    rhs = u[end] * v[end] - u[1] * v[1]
    @test isapprox(lhs, rhs; atol = 100eps(T), rtol = 100eps(T))

    f2 = repeat(f, 1, 3)
    out2 = similar(f2)
    sbp_deriv_x!(out2, f2, s)
    @test out2[:, 1] ≈ out
    @test out2[:, 2] ≈ out

    ny = 24
    Ly = T(2π)
    y = range(zero(T), Ly; length = ny + 1)[1:end-1]
    fy = repeat(reshape(collect(sin.(3 .* y)), 1, ny), 4, 1)
    dy = similar(fy)
    fourier_deriv_y!(dy, fy, Ly)
    @test norm(dy .- repeat(reshape(collect(3 .* cos.(3 .* y)), 1, ny), 4, 1)) /
          norm(dy) < 1e-12
end

@testset "defensive shape and alias checks" begin
    T = Float64
    g = FourierGrid((8, 8), (T(2π), T(2π)))
    f = randn(T, g.n...)
    out = similar(f)

    @test_throws ArgumentError FourierGrid((0,), (T(1),))
    @test_throws ArgumentError FourierGrid((8,), (zero(T),))
    @test_throws DimensionMismatch FourierGrid((8, 8), (T(1),))
    @test_throws DimensionMismatch FourierGrid((8,), (T(1), T(1)))
    @test_throws ArgumentError deriv!(out, f, g, 3)
    @test_throws DimensionMismatch deriv!(zeros(T, 7, 8), f, g, 1)
    @test_throws DimensionMismatch gradient!((zeros(T, 7, 8), similar(f)), f, g)
    @test_throws ArgumentError gradient!((out, out), f, g)
    @test_throws ArgumentError divergence!(similar(f), (f,), g)
    @test_throws DimensionMismatch laplacian!(zeros(T, 7, 8), f, g)

    A = ntuple(_ -> randn(T, g.n...), 3)
    B = ntuple(_ -> similar(f), 3)
    curl!(B, A, g)
    @test_throws ArgumentError curl!(A, A, g)
    @test_throws ArgumentError curl!((B[1], B[1], B[3]), A, g)
    @test_throws DimensionMismatch curl!((zeros(T, 7, 8), B[2], B[3]), A, g)

    P = ntuple(_ -> randn(T, g.n...), 3)
    project_divfree!(P, g)
    @test_throws ArgumentError project_divfree!((P[1], P[1], P[3]), g)
    @test_throws DimensionMismatch project_divfree!((zeros(T, 7, 8), P[2], P[3]), g)

    @test_throws ArgumentError SBP1D(2, T(1))
    @test_throws ArgumentError SBP1D(4, zero(T))
    @test_throws ArgumentError SBP1D(4, -one(T))
    s = SBP1D(4, T(1))
    @test_throws DimensionMismatch sbp_deriv!(zeros(T, 3), zeros(T, 4), s)
    @test_throws DimensionMismatch sbp_deriv!(zeros(T, 4), zeros(T, 3), s)
    v = zeros(T, 4)
    @test_throws ArgumentError sbp_deriv!(v, v, s)
    offset_v = OffsetSbpVector(zeros(T, 4), -2)
    @test_throws DimensionMismatch sbp_deriv!(zeros(T, 4), offset_v, s)
    @test_throws DimensionMismatch sbp_deriv!(offset_v, zeros(T, 4), s)
    @test_throws DimensionMismatch sbp_deriv_x!(zeros(T, 3, 2), zeros(T, 4, 2), s)
    @test_throws DimensionMismatch sbp_deriv_x!(zeros(T, 4, 2), zeros(T, 3, 2), s)
    m = zeros(T, 4, 2)
    @test_throws ArgumentError sbp_deriv_x!(m, m, s)
    offset_m = OffsetSbpMatrix(zeros(T, 4, 2), -1, 1)
    @test_throws DimensionMismatch sbp_deriv_x!(zeros(T, 4, 2), offset_m, s)
    @test_throws DimensionMismatch sbp_deriv_x!(offset_m, zeros(T, 4, 2), s)
    offset_cols = OffsetSbpMatrix(zeros(T, 4, 2), 1, 0)
    @test_throws DimensionMismatch sbp_deriv_x!(zeros(T, 4, 2), offset_cols, s)
    @test_throws ArgumentError fourier_deriv_y!(zeros(T, 4, 2), zeros(T, 4, 2), zero(T))
    @test_throws DimensionMismatch fourier_deriv_y!(zeros(T, 4, 1), zeros(T, 4, 2), T(1))
    @test_throws ArgumentError FourierDerivYWorkspace(0, 2, T(1))
    @test_throws ArgumentError FourierDerivYWorkspace(2, 0, T(1))
    @test_throws ArgumentError FourierDerivYWorkspace(2, 2, zero(T))
    wy = FourierDerivYWorkspace(4, 2, T(1))
    @test_throws DimensionMismatch fourier_deriv_y!(zeros(T, 4, 2), zeros(T, 3, 2), wy)
    @test_throws DimensionMismatch fourier_deriv_y!(zeros(T, 4, 1), zeros(T, 4, 2), wy)

    sw = BinomialSmoothWorkspace(g)
    @test_throws ArgumentError binomial_smooth!(copy(f), g, sw; passes = -1)
    @test_throws DimensionMismatch binomial_smooth!(zeros(T, 7, 8), g, sw; passes = 1)
    @test_throws DimensionMismatch binomial_smooth!(copy(f), g, BinomialSmoothWorkspace{T}(zeros(T, 1)); passes = 1)
end

@testset "operators allocate nothing in steady state" begin
    T = Float64
    n = (16, 16, 16)
    g = FourierGrid(n, ntuple(_ -> T(2π), 3))
    f = randn(T, n...)
    out = similar(f)
    grad = ntuple(_ -> similar(f), 3)
    A = ntuple(_ -> randn(T, n...), 3)
    B = ntuple(_ -> similar(f), 3)
    divB = similar(f)
    lap = similar(f)
    filt = copy(f)
    smooth = copy(f)
    smooth_work = BinomialSmoothWorkspace(g)
    my = randn(T, 16, 16)
    dy = similar(my)
    ywork = FourierDerivYWorkspace(my, T(2π))
    s = SBP1D(16, T(1))
    fx = randn(T, 16)
    dfx = similar(fx)
    mx = randn(T, 16, 8)
    dmx = similar(mx)
    deriv!(out, f, g, 1)
    gradient!(grad, f, g)
    curl!(B, A, g)
    divergence!(divB, B, g)
    laplacian!(lap, f, g)
    exp_filter!(filt, g)
    filt .= f
    dealias_two_thirds!(filt, g)
    binomial_smooth!(smooth, g, smooth_work; passes = 1)
    fourier_deriv_y!(dy, my, ywork)
    sbp_deriv!(dfx, fx, s)
    sbp_deriv_x!(dmx, mx, s)
    project_divfree!(B, g)   # warm up
    @test (@allocated deriv!(out, f, g, 1)) == 0
    @test (@allocated gradient!(grad, f, g)) == 0
    @test (@allocated curl!(B, A, g)) == 0
    @test (@allocated divergence!(divB, B, g)) == 0
    @test (@allocated laplacian!(lap, f, g)) == 0
    filt .= f
    @test (@allocated exp_filter!(filt, g)) == 0
    filt .= f
    @test (@allocated dealias_two_thirds!(filt, g)) == 0
    smooth .= f
    @test (@allocated binomial_smooth!(smooth, g, smooth_work; passes = 1)) == 0
    @test (@allocated binomial_smooth!(smooth, g, smooth_work; passes = 0)) == 0
    @test (@allocated fourier_deriv_y!(dy, my, ywork)) == 0
    @test (@allocated sbp_deriv!(dfx, fx, s)) == 0
    @test (@allocated sbp_deriv_x!(dmx, mx, s)) == 0
    @test (@allocated project_divfree!(B, g)) == 0
end

@testset "laplacian! includes the Nyquist mode (even N)" begin
    # Regression: ∇² is an even derivative, so the Nyquist mode must get −k_nyq²·f,
    # not 0. The first-derivative kvec zeroes Nyquist; laplacian! must not reuse it.
    g = FourierGrid((8,), (2π,))
    f = Float64[(-1)^(i - 1) for i = 1:8]          # pure Nyquist mode
    out = similar(f)
    laplacian!(out, f, g)
    knyq = 2π * (8 ÷ 2) / 2π                        # = 4
    @test out ≈ (-knyq^2) .* f rtol = 1e-10        # was identically 0 before the fix
end
