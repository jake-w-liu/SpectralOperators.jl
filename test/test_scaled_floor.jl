const SOFloor = SpectralOperators

scaled_floor_bits(x::Float32) = reinterpret(UInt32, x)
scaled_floor_bits(x::Float64) = reinterpret(UInt64, x)

# Independent integer definition: the smallest power of two at least n.
function scaled_floor_reference_shift(n::Int)
    power, shift = big(1), 4
    while power < n
        power *= 2
        shift += 1
    end
    return shift
end

@testset "Scaled Laplacian input floor exact IEEE values" begin
    counts = Set{Int}((typemin(Int), -1, 0, 1, 2, typemax(Int) - 1, typemax(Int)))
    for power = 0:(8sizeof(Int)-2)
        n = one(Int) << power
        push!(counts, n - 1, n, n + 1)
    end
    for T in (Float32, Float64), n in sort!(collect(counts))
        shift = scaled_floor_reference_shift(n)
        exact = Rational{BigInt}(floatmin(T)) * big(2)^shift
        old = ldexp(floatmin(T), SOFloor._ceil_log2(n) + 4)
        floor = SOFloor._laplacian_input_floor(T, n)
        @test isfinite(floor)
        @test floor >= floatmin(T)
        @test Rational{BigInt}(floor) == exact
        @test scaled_floor_bits(floor) == scaled_floor_bits(old)
    end
end

@testset "Scaled Laplacian input floor routing neighbors" begin
    setprecision(BigFloat, 8192) do
        for T in (Float32, Float64), n in (1, 4, 5)
            g = FourierGrid((n,), (T(2pi),))
            floor = SOFloor._laplacian_input_floor(T, n)
            for magnitude in
                (zero(T), nextfloat(zero(T)), prevfloat(floor), floor, nextfloat(floor)),
                sign in (-one(T), one(T))

                f = zeros(T, n)
                f[1] = sign * magnitude
                reference, scale = scaled_reference(g, f, :laplacian, 1)
                counted, calls = scaled_counting_grid(g)
                data = SOFloor.laplacian_scaled_data!(counted, f)
                @test data.values === g.cbuf
                @test data.ordinary == (iszero(magnitude) || magnitude >= floor)
                @test calls[] == 2
                @test isfinite(data.input_scale) && data.input_scale >= one(T)
                @test all(isfinite, data.values)
                @test maximum(abs, scaled_restore(data) .- reference) <= 2048eps(T) * scale
            end
        end
    end
end

struct ScaledFloorForbiddenPlan
    calls::Base.RefValue{Int}
end

function Base.:*(plan::ScaledFloorForbiddenPlan, data)
    plan.calls[] += 1
    error("unsupported grid reached its forbidden FFT plan")
end

# The public struct constructor permits unsupported manually assembled types.
# These plans must never run: this fixture tests rejection, not FFT support.
function scaled_floor_unsupported_grid(::Type{T}) where {T}
    calls = Ref(0)
    plan = ScaledFloorForbiddenPlan(calls)
    g = FourierGrid(
        (2,),
        (T(2),),
        (one(T),),
        ([0, 1],),
        ([zero(T), one(T)],),
        (zeros(Complex{T}, 2),),
        (zeros(T, 2),),
        plan,
        plan,
        fill(complex(T(1), T(2)), 2),
        fill(complex(T(3), T(4)), 2),
        fill(complex(T(5), T(6)), 2),
        fill(complex(T(7), T(8)), 2),
    )
    return g, calls
end

function scaled_floor_grid_values(g)
    return deepcopy((g.n, g.L, g.dx, g.midx, g.kfull, g.ik, g.kvec, g.cbuf, g.tbuf, g.abuf, g.sbuf))
end

@testset "Unsupported manual grid rejection precedes scratch and FFT work" begin
    for T in (Float16, BigFloat, Int)
        g, calls = scaled_floor_unsupported_grid(T)
        f = ones(T, g.n)
        before, input_before = scaled_floor_grid_values(g), copy(f)
        if T <: AbstractFloat
            @test_throws ArgumentError SOFloor.laplacian_scaled_data!(g, f)
        else
            # The canonical FFTW validator is defined for AbstractFloat only.
            @test_throws MethodError SOFloor.laplacian_scaled_data!(g, f)
        end
        @test calls[] == 0
        @test isequal(scaled_floor_grid_values(g), before)
        @test isequal(f, input_before)
        @test_throws MethodError SOFloor._laplacian_input_floor(T, length(f))
    end
end
