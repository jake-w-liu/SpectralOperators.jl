using Random, LinearAlgebra

const SOScaled = SpectralOperators
const scaled_error_receipts = Dict("local_bound_fraction" => 0.0, "DFT_source_error_eps" => 0.0)

function scaled_restore(data)
    return ldexp.(BigFloat.(real.(data.values)) .* BigFloat(data.input_scale), data.exponent)
end

# Independent real-space collocation matrices formed from actual stored metadata.
function scaled_matrix(g, axis, kind)
    n = g.n[axis]
    return [
        sum(1:n) do m
            phase = 2BigFloat(pi) * (m - 1) * (i - j) / n
            k = BigFloat(g.kfull[axis][m])
            kc = BigFloat(g.kvec[axis][m])
            kind === :laplacian ? -k * k * cos(phase) / n :
            kind === :derivative ? -kc * sin(phase) / n :
            (iseven(n) && m == n ÷ 2 + 1 ? k * cos(phase) / n : BigFloat(0))
        end for i = 1:n, j = 1:n
    ]
end

function scaled_apply_matrix(matrix, f, axis)
    result = zeros(BigFloat, size(f))
    for I in CartesianIndices(f), j in axes(f, axis)
        J = CartesianIndex(ntuple(d -> d == axis ? j : I[d], ndims(f)))
        result[I] += matrix[I[axis], j] * BigFloat(f[J])
    end
    return result
end

scaled_matrix_gain(matrix) = maximum(sum(abs, row) for row in eachrow(matrix))

function scaled_reference(g, f, kind, axis)
    if kind === :laplacian || kind === :derivative_laplacian
        matrices = [scaled_matrix(g, d, :laplacian) for d = 1:ndims(f)]
        result = sum(scaled_apply_matrix(matrices[d], f, d) for d = 1:ndims(f))
        gain = sum(scaled_matrix_gain, matrices)
        if kind === :derivative_laplacian
            first = scaled_matrix(g, axis, :derivative)
            result = scaled_apply_matrix(first, result, axis)
            gain *= scaled_matrix_gain(first)
        end
        return result, gain * maximum(abs, BigFloat.(f))
    end
    matrix = scaled_matrix(g, axis, kind)
    return scaled_apply_matrix(matrix, f, axis),
    scaled_matrix_gain(matrix) * maximum(abs, BigFloat.(f))
end

function scaled_call(g, f, kind, axis)
    kind === :laplacian && return SOScaled.laplacian_scaled_data!(g, f)
    kind === :derivative && return SOScaled.deriv_scaled_data!(g, f, axis)
    kind === :nyquist && return SOScaled.nyquist_scaled_data!(g, f, axis)
    return SOScaled.deriv_laplacian_scaled_data!(g, f, axis)
end

struct ScaledCountingPlan{P}
    plan::P
    calls::Base.RefValue{Int}
end

function Base.:*(plan::ScaledCountingPlan, data)
    plan.calls[] += 1
    return plan.plan * data
end

function scaled_counting_grid(g)
    calls = Ref(0)
    counted = SOScaled.FourierGrid(
        g.n,
        g.L,
        g.dx,
        g.midx,
        g.kfull,
        g.ik,
        g.kvec,
        ScaledCountingPlan(g.plan, calls),
        ScaledCountingPlan(g.iplan, calls),
        g.cbuf,
        g.tbuf,
        g.abuf,
        g.sbuf,
    )
    return counted, calls
end

# Exact products of the represented source, independent of the implementation's
# mantissa normalization. No FFT is involved in this local arithmetic oracle.
function scaled_exact_terms(value, g, I, kind, axis, power)
    x = BigFloat(value)
    factor = ldexp(BigFloat(1), -power)
    if kind === :laplacian || kind === :derivative_laplacian
        odd = kind === :derivative_laplacian ? BigFloat(g.kvec[axis][I[axis]]) : BigFloat(1)
        return [-x * BigFloat(g.kfull[d][I[d]])^2 * odd * factor for d = 1:length(g.n)]
    end
    k =
        kind === :nyquist ?
        (
            iseven(g.n[axis]) && I[axis] == g.n[axis] ÷ 2 + 1 ? BigFloat(g.kfull[axis][I[axis]]) :
            BigFloat(0)
        ) : BigFloat(g.kvec[axis][I[axis]])
    return [x * k * factor]
end

@testset "Scaled multiplier local dimensioned error bound" begin
    setprecision(BigFloat, 8192) do
        for T in (Float32, Float64), D = 1:3, span in (-1, 0, 1)
            e = T === Float32 ? 50 : 500
            g = FourierGrid(
                ntuple(_ -> 4, D),
                ntuple(d -> ldexp(T(2pi) * T(1 + d / 5), span * (isodd(d) ? e : -e)), D),
            )
            for (i, I) in enumerate(CartesianIndices(g.cbuf))
                ex =
                    i % 3 == 0 ? exponent(nextfloat(zero(T))) :
                    i % 3 == 1 ? (T === Float32 ? 100 : 1000) : -(T === Float32 ? 100 : 1000)
                g.cbuf[I] = complex(ldexp(T(1.3), ex), ldexp(T(-1.7), ex))
            end
            source = copy(g.cbuf)
            for kind in (:laplacian, :derivative, :nyquist, :derivative_laplacian), axis = 1:D

                copyto!(g.cbuf, source)
                power = SOScaled._apply_normalized_multiplier!(g, axis, Val(kind))
                u = BigFloat(eps(T)) / 2
                gamma = (D + 3) * u / (1 - (D + 3) * u)
                omission = D * BigFloat(nextfloat(zero(T))) / 2
                for I in CartesianIndices(source)
                    r = scaled_exact_terms(real(source[I]), g, I, kind, axis, power)
                    i = scaled_exact_terms(imag(source[I]), g, I, kind, axis, power)
                    expected_r, expected_i = sum(r), sum(i)
                    if kind === :derivative || kind === :derivative_laplacian
                        expected_r, expected_i = -expected_i, expected_r
                        r, i = i, r
                    end
                    for (actual, expected, terms) in
                        ((real(g.cbuf[I]), expected_r, r), (imag(g.cbuf[I]), expected_i, i))
                        bound = gamma * sum(abs, terms) + omission
                        error = abs(BigFloat(actual) - expected)
                        @test error <= bound
                        if !iszero(bound)
                            scaled_error_receipts["local_bound_fraction"] = max(
                                scaled_error_receipts["local_bound_fraction"],
                                Float64(error / bound),
                            )
                        end
                    end
                end
            end
        end
    end
end

@testset "Scaled Fourier sources against independent 8192-bit DFT" begin
    setprecision(BigFloat, 8192) do
        rng = Xoshiro(0x734c0)
        for T in (Float32, Float64), D = 1:3, variant = 1:3
            n = ntuple(
                d ->
                    variant == 1 ? (isodd(d) ? 4 : 5) :
                    variant == 2 ? (d == D ? 1 : 3) : (d == D ? 2 : 5),
                D,
            )
            g = FourierGrid(n, ntuple(d -> T(2pi) * T(1 + d / 4), D))
            f = randn(rng, T, n)
            for kind in (:laplacian, :derivative, :nyquist, :derivative_laplacian), axis = 1:D

                reference, scale = scaled_reference(g, f, kind, axis)
                data = scaled_call(g, f, kind, axis)
                actual = scaled_restore(data)
                @test data.values === g.cbuf
                @test isfinite(data.input_scale) && data.input_scale >= one(T)
                @test all(isfinite, data.values)
                error = maximum(abs, actual .- reference)
                @test error <= 2048eps(T) * scale
                if !iszero(scale)
                    scaled_error_receipts["DFT_source_error_eps"] = max(
                        scaled_error_receipts["DFT_source_error_eps"],
                        Float64(error / (BigFloat(eps(T)) * scale)),
                    )
                end
            end
            old = similar(f)
            laplacian!(old, f, g)
            data = SOScaled.laplacian_scaled_data!(g, f)
            @test data.ordinary
            @test data.exponent == 0
            restored = T.(real.(data.values)) .* data.input_scale
            @test isequal(restored, old)
            loss = -sum(BigFloat.(f) .* scaled_restore(data))
            heat = BigFloat(0)
            for d = 1:D
                heat += sum(abs2, scaled_restore(SOScaled.deriv_scaled_data!(g, f, d)))
                heat += sum(abs2, scaled_restore(SOScaled.nyquist_scaled_data!(g, f, d)))
            end
            @test abs(loss - heat) <= 512eps(T) * max(loss, heat)
        end
    end
end

@testset "Scaled ordinary Laplacian preserves complete legacy profiles" begin
    for T in (Float32, Float64),
        D = 1:3,
        n in (4, 5),
        input_power in (-10, 0, 10),
        grid_power in (-20, 0, 20)

        g = FourierGrid(ntuple(_ -> n, D), ntuple(_ -> ldexp(T(2pi), grid_power), D))
        rng = Xoshiro(0x772 + D + n)
        f = ldexp.(randn(rng, T, g.n), input_power)
        old = similar(f)
        laplacian!(old, f, g)
        data = SOScaled.laplacian_scaled_data!(g, f)
        @test data.ordinary
        @test data.exponent == 0
        restored = similar(f)
        SOScaled._store_fft_output!(restored, data.values, data.input_scale)
        @test isequal(restored, old)
    end
end

@testset "Scaled sources retain coefficient-rescued and anisotropic modes" begin
    setprecision(BigFloat, 8192) do
        for T in (Float32, Float64)
            e = T === Float32 ? 50 : 500
            coefficient = ldexp(BigFloat(1), 2e)
            pattern = T[0, 1, 0, -1]
            f = ldexp.(pattern, -e)
            g = FourierGrid((4,), (ldexp(T(2pi), e),))
            old = similar(f)
            laplacian!(old, f, g)
            @test all(iszero, old)
            data = SOScaled.laplacian_scaled_data!(g, f)
            @test !data.ordinary
            @test scaled_restore(data) .* coefficient == -BigFloat.(f)
            anisotropic = FourierGrid((4, 2), (ldexp(T(2pi), e), ldexp(T(2pi), -e)))
            broad = repeat(f, 1, 2)
            @test scaled_restore(SOScaled.laplacian_scaled_data!(anisotropic, broad)) .*
                  coefficient == -BigFloat.(broad)

            ek = T === Float32 ? 24 : 200
            ej = T === Float32 ? -101 : -674
            eh = T === Float32 ? -106 : -810
            stiff = FourierGrid((4,), (ldexp(T(2pi), -ek),))
            j = ldexp.(pattern, ej)
            laplacian!(old, j, stiff)
            @test all(iszero, ldexp.(old, eh))
            complete =
                scaled_restore(SOScaled.deriv_laplacian_scaled_data!(stiff, j, 1)) .*
                ldexp(BigFloat(1), eh)
            reference = ldexp.(BigFloat[-1, 0, 1, 0], ej + 3ek + eh)
            @test complete == reference
            @test any(!iszero, T.(complete))
        end
    end
end

@testset "Scaled source validation, range and borrowed lifetime" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(2pi),))
        f = T[0, 1, 0, -1]
        for kind in (:laplacian, :derivative, :nyquist, :derivative_laplacian)
            fill!(g.cbuf, complex(T(7), T(3)))
            saved = copy(g.cbuf)
            @test_throws DimensionMismatch scaled_call(g, zeros(T, 3), kind, 1)
            @test isequal(g.cbuf, saved)
            @test_throws ArgumentError scaled_call(g, fill(T(NaN), 4), kind, 1)
            @test isequal(g.cbuf, saved)
            @test_throws ArgumentError scaled_call(
                g,
                zeros(T === Float32 ? Float64 : Float32, 4),
                kind,
                1,
            )
            @test isequal(g.cbuf, saved)
            alias = @view reinterpret(T, g.cbuf)[1:4]
            @test_throws ArgumentError scaled_call(g, alias, kind, 1)
            @test isequal(g.cbuf, saved)
            if kind !== :laplacian
                @test_throws ArgumentError scaled_call(g, f, kind, 0)
                @test_throws ArgumentError scaled_call(g, f, kind, 2)
                @test isequal(g.cbuf, saved)
            end
        end
        metadata = map(copy, g.kfull)
        @test all(isfinite, SOScaled.laplacian_scaled_data!(g, g.kfull[1]).values)
        @test isequal(metadata, g.kfull)
        data = SOScaled.laplacian_scaled_data!(g, f)
        @test any(!iszero, data.values)
        SOScaled.nyquist_scaled_data!(g, f, 1)
        @test all(iszero, data.values) # the documented borrowed lifetime

        tiny = ldexp.(f, exponent(nextfloat(zero(T))))
        first = scaled_restore(SOScaled.deriv_scaled_data!(g, tiny, 1))
        @test first == ldexp.(BigFloat[1, 0, -1, 0], exponent(nextfloat(zero(T))))
        wide = T[floatmax(T)/2, nextfloat(zero(T)), 0, 0]
        @test_throws ArgumentError SOScaled.deriv_scaled_data!(g, wide, 1)
        @test_throws ArgumentError SOScaled.laplacian_scaled_data!(g, wide)

        e = T === Float32 ? 70 : 600
        invalid_square = FourierGrid((4,), (ldexp(T(2pi), -e),))
        saved = copy(invalid_square.cbuf)
        @test_throws ArgumentError SOScaled.laplacian_scaled_data!(invalid_square, f)
        @test isequal(saved, invalid_square.cbuf)
        @test_throws ArgumentError SOScaled.deriv_laplacian_scaled_data!(invalid_square, f, 1)
        @test isequal(saved, invalid_square.cbuf)
        @test all(isfinite, SOScaled.deriv_scaled_data!(invalid_square, f, 1).values)

        aliased_metadata = SOScaled.FourierGrid(
            g.n,
            g.L,
            g.dx,
            g.midx,
            g.kfull,
            g.ik,
            g.kvec,
            g.plan,
            g.iplan,
            g.ik[1],
            g.tbuf,
            g.abuf,
            g.sbuf,
        )
        saved = copy(g.ik[1])
        @test_throws ArgumentError SOScaled.laplacian_scaled_data!(aliased_metadata, f)
        @test isequal(saved, g.ik[1])
    end
end

@testset "Scaled source bounded transforms and exact null channels" begin
    for T in (Float32, Float64)
        ordinary, count = scaled_counting_grid(FourierGrid((4,), (T(2pi),)))
        pattern = T[0, 1, 0, -1]
        @test SOScaled.laplacian_scaled_data!(ordinary, pattern).ordinary
        @test count[] == 2
        count[] = 0
        tiny = ldexp.(pattern, exponent(nextfloat(zero(T))))
        @test !SOScaled.laplacian_scaled_data!(ordinary, tiny).ordinary
        @test count[] == 2 # pre-forward tiny-input routing, no discarded transform

        e = T === Float32 ? 50 : 500
        cold, count = scaled_counting_grid(FourierGrid((4,), (ldexp(T(2pi), e),)))
        @test !SOScaled.laplacian_scaled_data!(cold, ldexp.(pattern, -e)).ordinary
        @test count[] == 3 # one old forward attempt, then one normalized pair
        for kind in (:derivative, :nyquist, :derivative_laplacian)
            count[] = 0
            scaled_call(cold, pattern, kind, 1)
            @test count[] == 2
        end

        # A physical-gain scale would discard the tiny value. The normalized
        # path selects unit-gain scaling before that lossy copy or forward FFT.
        stiff, count = scaled_counting_grid(FourierGrid((4,), (ldexp(T(2pi), -e),)))
        broad = T[floatmax(T)/1024, nextfloat(zero(T)), 0, 0]
        @test_throws ArgumentError laplacian!(similar(broad), broad, stiff)
        count[] = 0
        data = SOScaled.laplacian_scaled_data!(stiff, broad)
        @test !data.ordinary
        @test all(isfinite, data.values)
        @test count[] == 2

        for n in (1, 2, 5)
            null, count = scaled_counting_grid(FourierGrid((n,), (T(2pi),)))
            finite = fill(floatmax(T), n)
            finite[end] = nextfloat(zero(T))
            if n <= 2
                for kind in (:derivative, :derivative_laplacian)
                    @test all(iszero, scaled_call(null, finite, kind, 1).values)
                    @test count[] == 0
                end
            end
            if isodd(n)
                @test all(iszero, SOScaled.nyquist_scaled_data!(null, finite, 1).values)
                @test count[] == 0
                finite[1] = T(NaN)
                @test_throws ArgumentError SOScaled.nyquist_scaled_data!(null, finite, 1)
                @test count[] == 0
            end
        end
    end
end

@info "Scaled-source error receipts" scaled_error_receipts
