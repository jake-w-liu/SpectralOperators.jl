module ScaledAxisConditioningTests
using Test, SpectralOperators
const SO = SpectralOperators
const Q = Rational{BigInt}

restore(data) = map(value -> Q(real(value)) * Q(data.input_scale) * Q(2)^data.exponent, data.values)

function derivative(k)
    return [mod(j-i, 4) == 1 ? k/2 : mod(j-i, 4) == 3 ? -k/2 : zero(k) for i = 1:4, j = 1:4]
end
function laplacian(k)
    row = Q[-3//2, 1, -1//2, 1]
    return [k^2 * row[mod(j-i, 4)+1] for i = 1:4, j = 1:4]
end
function directional_call(g, f, axis, kind)
    kind === :derivative && return SO.deriv_scaled_data!(g, f, axis)
    kind === :nyquist && return SO.nyquist_scaled_data!(g, f, axis)
    return SO.deriv_laplacian_scaled_data!(g, f, axis)
end

@testset "Directional source removes only a certified axis-null field" begin
    for T in (Float32, Float64), axis = 1:2, kind in (:derivative, :nyquist, :derivative_laplacian)
        a = T === Float32 ? 100 : 600
        tiny = ldexp(one(T), -a)
        lengths = ntuple(d -> d == axis ? ldexp(T(2pi), a) : T(2pi), 2)
        g = FourierGrid((4, 4), lengths)
        @test g.kvec[axis][2] == tiny
        f = ones(T, g.n)
        pattern = kind === :nyquist ? T[3, 1, 3, 1] : T[2, 3, 2, 1]
        axis == 1 ? (f[:, 1] .= tiny .* pattern) : (f[1, :] .= tiny .* pattern)
        saved = copy(f)
        original = Q.(f)
        k = Q(g.kvec[axis][2])
        if kind === :nyquist
            op = [isodd(i-j) ? -k/2 : k/2 for i = 1:4, j = 1:4]
        else
            op = derivative(k)
        end
        source = kind === :derivative_laplacian ?
            laplacian(Q(g.kfull[1][2])) * original + original * transpose(laplacian(Q(g.kfull[2][2]))) : original
        reference = axis == 1 ? op * source : source * transpose(op)
        @test maximum(abs, reference) > 0
        data = directional_call(g, f, axis, kind)
        actual = restore(data)
        @test data.values === g.cbuf && !data.ordinary
        @test isequal(f, saved)
        @test maximum(abs, actual - reference) <= 8Q(eps(T)) * maximum(abs, reference)

        # These large transverse values are constant along the consumed axis.
        # Their subtraction is exact zero before any operand downscaling.
        constants = T[floatmax(T), nextfloat(zero(T)), -floatmax(T), -nextfloat(zero(T))]
        f = axis == 1 ? repeat(reshape(constants, 1, 4), 4, 1) : repeat(reshape(constants, 4, 1), 1, 4)
        @test all(iszero, directional_call(g, f, axis, kind).values)
    end
end

@testset "Conditioning retains exact factors and unsafe cases use old preparation" begin
    for T in (Float32, Float64)
        g = FourierGrid((4,), (T(2pi),))
        for (a, b) in ((one(T), nextfloat(one(T))), (nextfloat(zero(T)), -nextfloat(zero(T))),
                       (floatmax(T)/T(1024), -floatmax(T)/T(1024)), (one(T), eps(T)/T(4)))
            value, exact = SO._axis_difference(a, b)
            @test exact == (Q(value) == Q(a) - Q(b))
        end
        extreme = T[floatmax(T), -floatmax(T), floatmax(T), -floatmax(T)]
        scale, power = SO._prepare_normalized_input!(g, extreme, floatmax(T), 1, Val(:nyquist))
        @test all(isfinite, g.cbuf)
        @test Q.(real.(g.cbuf)) .* Q(scale) .* Q(2)^power == Q.(extreme) .- Q(extreme[1])
        @test restore(SO.nyquist_scaled_data!(g, extreme, 1)) == 2 .* Q.(extreme)

        # The stored difference between 1 and eps/4 is not representable.
        # Certification must not round it away while removing the reference.
        inexact = T[1, eps(T)/4, -1, 0]
        scale, power = SO._prepare_normalized_input!(g, inexact, one(T), 1, Val(:derivative))
        @test scale == one(T) && power == 0
        @test real.(g.cbuf) == inexact
        wide = T[floatmax(T)/2, nextfloat(zero(T)), 0, 0]
        @test_throws ArgumentError SO.deriv_scaled_data!(g, wide, 1)

        # Full Laplacian cannot discard a background varying in another axis.
        two = FourierGrid((4, 4), (T(2pi), T(2pi)))
        transverse = repeat(reshape(T[1, 2, 1, 0], 1, 4), 4, 1)
        expected = similar(transverse)
        laplacian!(expected, transverse, two)
        data = SO.laplacian_scaled_data!(two, transverse)
        @test data.ordinary && restore(data) == Q.(expected)
        @test any(!iszero, expected)

        fill!(g.cbuf, Complex{T}(3, 4))
        aliased = view(reinterpret(T, g.cbuf), 1:2:2length(g.cbuf))
        saved = copy(g.cbuf)
        @test_throws ArgumentError SO.deriv_scaled_data!(g, aliased, 1)
        @test isequal(g.cbuf, saved)
    end
end
end
