const SOZero = SpectralOperators

# Compare the old canonical fill, including exact zero signs and private scalar
# compatibility. Public FFTW admission remains covered by test_scaled_floor.jl.
function scaled_zero_neighbors(g)
    return deepcopy((g.n, g.L, g.dx, g.midx, g.kfull, g.ik, g.kvec, g.tbuf, g.abuf, g.sbuf))
end

@testset "Scaled zero buffer preserves original storage and exact fill semantics" begin
    context = (precision(BigFloat), rounding(BigFloat))
    shapes = (
        (),
        (0,),
        (1,),
        (4,),
        (0, 2),
        (2, 0),
        (1, 1),
        (4, 4),
        (0, 1, 2),
        (1, 1, 1),
        (2, 2, 2),
        (1, 1, 1, 1),
        (2, 2, 2, 2),
    )
    for T in (Float16, Float32, Float64, BigFloat, Int), shape in shapes
        seeds =
            T <: AbstractFloat ?
            (complex(T(1), -T(2)), complex(-zero(T), zero(T)), complex(T(NaN), T(Inf))) :
            (complex(T(1), -T(2)), complex(zero(T), zero(T)))
        for seed in seeds
            buffer = fill(seed, shape)
            original = buffer
            alias = vec(buffer)
            expected = fill!(copy(buffer), zero(eltype(buffer)))
            result = SOZero._zero_scaled_buffer!(buffer)
            @test result === original
            @test typeof(result) === typeof(expected)
            @test size(result) == shape
            @test ndims(result) == length(shape)
            @test isequal(result, expected)
            @test all(x -> isequal(x, complex(zero(T), zero(T))), buffer)
            @test all(x -> isequal(x, complex(zero(T), zero(T))), alias)
            if !isempty(buffer)
                # Both directions establish that the original storage survived;
                # empty arrays have no element through which to witness sharing.
                alias[1] = complex(T(3), -T(4))
                @test isequal(buffer[1], complex(T(3), -T(4)))
                buffer[1] = complex(-T(5), T(6))
                @test isequal(alias[1], complex(-T(5), T(6)))
            end
        end
    end
    @test (precision(BigFloat), rounding(BigFloat)) == context
end

@testset "Both scaled zero owners preserve complete neighboring grid state" begin
    for T in (Float32, Float64), D = 1:4, n in (1, 4)
        g, calls = scaled_counting_grid(FourierGrid(ntuple(_ -> n, D), ntuple(_ -> T(2pi), D)))
        original = g.cbuf
        fill!(g.tbuf, complex(T(11), T(12)))
        fill!(g.abuf, complex(T(13), T(14)))
        fill!(g.sbuf, complex(T(15), T(16)))
        neighbors = scaled_zero_neighbors(g)
        fill!(g.cbuf, complex(T(NaN), T(Inf)))
        data = SOZero._zero_scaled_data!(g)
        @test data.values === original
        @test data.input_scale === one(T)
        @test data.exponent == 0 && !data.ordinary
        @test all(x -> isequal(x, complex(zero(T), zero(T))), data.values)
        @test isequal(scaled_zero_neighbors(g), neighbors)
        @test calls[] == 0

        for kind in (:laplacian, :derivative, :nyquist, :derivative_laplacian), axis = 1:D
            # n=1 has no active metadata even with a nonzero represented DC
            # source. At n=4 the zero source has no active products on any axis.
            fill!(g.cbuf, n == 1 ? complex(T(7), -T(3)) : zero(Complex{T}))
            @test SOZero._apply_normalized_multiplier!(g, axis, Val(kind)) == 0
            @test g.cbuf === original
            @test all(x -> isequal(x, complex(zero(T), zero(T))), g.cbuf)
            @test isequal(scaled_zero_neighbors(g), neighbors)
            @test calls[] == 0
        end

        fill!(g.cbuf, complex(T(7), -T(3)))
        g.cbuf[1] = complex(T(NaN), zero(T))
        before = copy(g.cbuf)
        @test_throws ArgumentError SOZero._apply_normalized_multiplier!(g, 1, Val(:laplacian))
        @test isequal(g.cbuf, before)
        @test isequal(scaled_zero_neighbors(g), neighbors)
        @test calls[] == 0
    end
end

@testset "Public null sources retain validation, plan counts and borrowed lifetime" begin
    for T in (Float32, Float64), D = 1:3, n in (1, 2, 4)
        g, calls = scaled_counting_grid(FourierGrid(ntuple(_ -> n, D), ntuple(_ -> T(2pi), D)))
        original = g.cbuf
        neighbors = scaled_zero_neighbors(g)
        for magnitude in (zero(T), nextfloat(zero(T))),
            kind in (:laplacian, :derivative, :nyquist, :derivative_laplacian),
            axis = 1:D

            f = fill(magnitude, g.n)
            input_before = copy(f)
            calls[] = 0
            data = scaled_call(g, f, kind, axis)
            @test data.values === original
            @test all(iszero, data.values)
            @test all(isfinite, data.values)
            @test data.input_scale === one(T)
            @test data.ordinary == (kind === :laplacian && iszero(magnitude))
            null_axis = (kind === :derivative || kind === :derivative_laplacian) && n <= 2
            null_nyquist = kind === :nyquist && isodd(n)
            @test calls[] == (null_axis || null_nyquist ? 0 : 2)
            @test isequal(f, input_before)
            @test isequal(scaled_zero_neighbors(g), neighbors)
        end

        calls[] = 0
        fill!(g.cbuf, complex(T(7), -T(3)))
        before = copy(g.cbuf)
        wrong_shape = zeros(T, ntuple(d -> d == 1 ? n + 1 : n, D))
        @test_throws DimensionMismatch SOZero.laplacian_scaled_data!(g, wrong_shape)
        @test_throws ArgumentError SOZero.deriv_scaled_data!(g, zeros(T, g.n), 0)
        @test_throws ArgumentError SOZero.nyquist_scaled_data!(g, fill(T(NaN), g.n), 1)
        @test isequal(g.cbuf, before)
        @test isequal(scaled_zero_neighbors(g), neighbors)
        @test calls[] == 0
    end

    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(2pi),))
        f = T[0, 1, 0, -1]
        data = SOZero.deriv_scaled_data!(g, f, 1)
        @test any(!iszero, data.values)
        @test SOZero._zero_scaled_buffer!(g.cbuf) === data.values
        @test all(iszero, data.values)
        @test f == T[0, 1, 0, -1]
    end
end
