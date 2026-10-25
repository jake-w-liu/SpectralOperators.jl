# A dependency-free logical global view over independently permuted storage.
# dataids matters: raw kernels must see overlapping parent storage.
struct FourierLogicalTile{T,D,A<:AbstractArray{T,D}} <: AbstractArray{T,D}
    data::A
    ranges::NTuple{D,UnitRange{Int}}
end
Base.size(a::FourierLogicalTile) = map(length, a.ranges)
Base.axes(a::FourierLogicalTile) = a.ranges
Base.IndexStyle(::Type{<:FourierLogicalTile}) = IndexCartesian()
Base.dataids(a::FourierLogicalTile) = Base.dataids(a.data)
# Preserve this adapter's concrete storage/permutation type when Base must unalias.
Base.unaliascopy(a::FourierLogicalTile)::typeof(a) = deepcopy(a)
function Base.similar(
    a::FourierLogicalTile,
    ::Type{T},
    ranges::NTuple{D,UnitRange{Int}},
) where {T,D}
    return FourierLogicalTile(Array{T}(undef, map(length, ranges)), ranges)
end
@inline Base.getindex(a::FourierLogicalTile{T,D}, I::Vararg{Int,D}) where {T,D} =
    a.data[ntuple(d->I[d]-first(a.ranges[d])+1, D)...]
@inline function Base.setindex!(a::FourierLogicalTile{T,D}, x, I::Vararg{Int,D}) where {T,D}
    a.data[ntuple(d->I[d]-first(a.ranges[d])+1, D)...] = x
    return x
end
function fourier_tile(a::Array{T,D}, ranges; permuted = false) where {T,D}
    shape = map(length, ranges)
    data = if permuted
        permutation = ntuple(d -> D + 1 - d, D)
        PermutedDimsArray(Array{T}(undef, reverse(shape)), permutation)
    else
        Array{T}(undef, shape)
    end
    tile = FourierLogicalTile(data, ranges)
    for I in CartesianIndices(tile)
        tile[I] = a[I]
    end
    return tile
end
function tile_equals(tile, full)
    all(I -> isequal(tile[I], full[I]), CartesianIndices(tile))
end

# Independent O(N^2) transform: no FFT or package mode arithmetic is called.
function dense_fourier(a, inverse = false)
    n = size(a)
    N = length(a)
    [
        sum(
            a[J] * cis(
                (inverse ? 2π : -2π) *
                sum((I[d] - 1) * (J[d] - 1) / n[d] for d = 1:ndims(a)),
            ) for J in CartesianIndices(a)
        ) / (inverse ? N : 1) for I in CartesianIndices(a)
    ]
end
function independent_modes(n, L, odd)
    ntuple(length(n)) do d
        [
            odd && iseven(n[d]) && i == n[d] ÷ 2 + 1 ? 0.0 :
            2π / L[d] * (i - 1 <= n[d] ÷ 2 ? i - 1 : i - 1 - n[d]) for i = 1:n[d]
        ]
    end
end

@testset "Shared Fourier metadata and logical kernels" begin
    SOShared = SpectralOperators
    for T in (Float32, Float64), n in ((7,), (8,), (1,), (6, 5), (4, 3, 2))
        D = length(n)
        L = ntuple(d -> T(2π + d / 3), D)
        g = FourierGrid(n, L)
        m = SOShared._fourier_metadata(n, L)
        for name in (:n, :L, :dx, :midx, :kfull, :ik, :kvec)
            @test isequal(getproperty(m, name), getproperty(g, name))
        end
        @test isapprox(
            SOShared._laplacian_gain(m.kfull),
            sum(k -> maximum(abs2, k), m.kfull);
            rtol = 4eps(T),
            atol = zero(T),
        )
        hats = ntuple(
            c -> [
                Complex{T}(T(sin(c * sum(Tuple(I)))), T(cos(c * prod(Tuple(I)))))
                for I in CartesianIndices(n)
            ],
            3,
        )
        ranges = ntuple(d -> max(1, n[d]÷2):n[d], D)
        for permuted in (false, true)
            tiles = map(a -> fourier_tile(a, ranges; permuted), hats)
            for tile in tiles
                copied = copy(tile)
                @test axes(copied) == axes(tile)
                @test !Base.mightalias(copied, tile)
                unaliased = Base.unaliascopy(tile)
                @test typeof(unaliased) === typeof(tile)
                @test !Base.mightalias(unaliased, tile)
                @test all(I -> isequal(unaliased[I], tile[I]), CartesianIndices(tile))
                @test all(I -> isequal(copied[I], tile[I]), CartesianIndices(tile))
            end
            for j = 1:D
                for (f, args) in (
                    (SOShared._apply_ik_store!, ()),
                    (SOShared._apply_ik_store_scaled!, (T(-2),)),
                    (SOShared._apply_ik_accum!, (T(-2),)),
                )
                    full = copy(hats[2])
                    tile = fourier_tile(full, ranges; permuted)
                    f(full, hats[1], m.ik[j], j, args...)
                    f(tile, tiles[1], m.ik[j], j, args...)
                    @test tile_equals(tile, full)
                    alias = copy(hats[1])
                    expected = copy(alias)
                    f(expected, hats[1], m.ik[j], j, args...)
                    @test f(alias, alias, m.ik[j], j, args...) === alias
                    @test isequal(alias, expected)
                end
            end
            full = copy(hats[1])
            tile = fourier_tile(full, ranges; permuted)
            SOShared._curl_component_store!(full, hats[1], m.ik[1], 1, hats[2], m.ik[D], D)
            SOShared._curl_component_store!(
                tile,
                tiles[1],
                m.ik[1],
                1,
                tiles[2],
                m.ik[D],
                D,
            )
            @test tile_equals(tile, full)
            for (f, args) in (
                (SOShared._apply_laplacian_multiplier!, (m.kfull,)),
                (SOShared._apply_exp_filter_multiplier!, (m.midx, n, T(5), 3)),
                (SOShared._apply_two_thirds_mask!, (m.midx, n)),
            )
                full = copy(hats[1])
                tile = fourier_tile(full, ranges; permuted)
                @test f(full, args...) === full
                @test f(tile, args...) === tile
                @test tile_equals(tile, full)
            end
            full = map(copy, hats)
            SOShared._project_fourier_modes!(full, m.kvec)
            @test SOShared._project_fourier_modes!(tiles, m.kvec) === tiles
            @test all(c -> tile_equals(tiles[c], full[c]), 1:3)
            for c = (D+1):3
                @test isequal(full[c], hats[c])
            end
        end
    end
end

@testset "Independent dense-DFT public operators" begin
    for T in (Float32, Float64), n in ((5,), (6,), (1,), (4, 3), (3, 2, 4))
        D = length(n)
        L = ntuple(d -> T(2π + d / 3), D)
        g = FourierGrid(n, L)
        f = [T(sin(sum(Tuple(I))) + cos(prod(Tuple(I)))) for I in CartesianIndices(n)]
        v = ntuple(c -> [T(sin(c * sum(Tuple(I)))) for I in CartesianIndices(n)], 3)
        fh = dense_fourier(f)
        vh = map(dense_fourier, v)
        k = independent_modes(n, L, true)
        kf = independent_modes(n, L, false)
        tol = 128length(f) * eps(T)
        out = similar(f)
        for j = 1:D
            ref = real.(
                dense_fourier([im * k[j][I[j]] * fh[I] for I in CartesianIndices(n)], true),
            )
            deriv!(out, f, g, j)
            @test isapprox(out, ref; atol = tol, rtol = tol)
        end
        ref = real.(
            dense_fourier(
                [-sum(kf[d][I[d]]^2 for d = 1:D) * fh[I] for I in CartesianIndices(n)],
                true,
            ),
        )
        laplacian!(out, f, g)
        @test isapprox(out, ref; atol = tol, rtol = tol)
        projected = map(copy, vh)
        for I in CartesianIndices(n)
            kv = ntuple(d -> k[d][I[d]], D)
            k2 = sum(abs2, kv)
            if k2 > 0
                dot = sum(kv[d] * vh[d][I] for d = 1:D)
                for d = 1:D
                    projected[d][I] -= kv[d] * dot / k2
                end
            end
        end
        actual = map(copy, v)
        project_divfree!(actual, g)
        for c = 1:3
            ref = real.(dense_fourier(projected[c], true))
            @test isapprox(actual[c], ref; atol = tol, rtol = tol)
        end
        modes = ntuple(d -> [i - 1 <= n[d] ÷ 2 ? i - 1 : i - 1 - n[d] for i = 1:n[d]], D)
        for filter in (:exp, :thirds)
            refhat = copy(fh)
            for I in CartesianIndices(n)
                if filter == :exp
                    r = maximum(
                        d -> n[d] == 1 ? 0.0 : abs(modes[d][I[d]]) / (n[d] ÷ 2),
                        1:D,
                    )
                    refhat[I] *= exp(-5 * r^6)
                elseif any(d -> 3abs(modes[d][I[d]]) >= n[d], 1:D)
                    refhat[I] = 0
                end
            end
            copyto!(out, f)
            filter == :exp ? exp_filter!(out, g; α = 5, p = 3) : dealias_two_thirds!(out, g)
            @test isapprox(out, real.(dense_fourier(refhat, true)); atol = tol, rtol = tol)
        end
    end
end

@testset "Shared kernels reject invalid ranges and aliases before writing" begin
    SO = SpectralOperators
    m = SO._fourier_metadata((8,), (8.0,))
    a = fill(1.0 + 2.0im, 8)
    b = copy(a)
    for j in (0, 2)
        @test_throws ArgumentError SO._apply_ik_store!(b, a, m.ik[1], j)
        @test isequal(b, a)
    end
    short = fill(3.0 + 4.0im, 7)
    saved = copy(short)
    @test_throws DimensionMismatch SO._apply_ik_store!(short, a, m.ik[1], 1)
    @test isequal(short, saved)
    bad = FourierLogicalTile(copy(a), (2:9,))
    saved = copy(bad.data)
    @test_throws DimensionMismatch SO._apply_ik_store!(bad, bad, m.ik[1], 1)
    @test isequal(bad.data, saved)
    # Equal logical axes but overlapping physical storage with a shifted mapping.
    storage = fill(1.0 + 2.0im, 9)
    x = FourierLogicalTile(view(storage, 1:8), (1:8,))
    y = FourierLogicalTile(view(storage, 2:9), (1:8,))
    saved = copy(storage)
    @test Base.mightalias(x, y)
    @test_throws ArgumentError SO._apply_ik_store!(x, y, m.ik[1], 1)
    @test isequal(storage, saved)
    @test_throws ArgumentError SO._apply_ik_store!(m.ik[1], a, m.ik[1], 1)
    @test_throws ArgumentError SO._project_fourier_modes!((a, a, b), m.kvec)
    @test isequal(a, b)
    for (f, args) in (
        (SO._apply_laplacian_multiplier!, (m.kfull,)),
        (SO._apply_exp_filter_multiplier!, (m.midx, m.n, 5.0, 3)),
        (SO._apply_two_thirds_mask!, (m.midx, m.n)),
    )
        @test_throws DimensionMismatch f(bad, args...)
        @test isequal(bad.data, fill(1.0 + 2.0im, 8))
    end
    for (alpha, p) in ((-1.0, 3), (NaN, 3), (1.0, 0))
        @test_throws ArgumentError SO._apply_exp_filter_multiplier!(
            a,
            m.midx,
            m.n,
            alpha,
            p,
        )
        @test isequal(a, b)
    end
    # Empty tiles do no reads/writes, even when their empty axes start outside metadata.
    empty = FourierLogicalTile(ComplexF64[], (10:9,))
    @test SO._apply_laplacian_multiplier!(empty, m.kfull) === empty
    @test SO._apply_ik_store!(empty, empty, m.ik[1], 1) === empty
end

function shared_kernel_allocations(f::F, args::A) where {F,A<:Tuple}
    for _ = 1:3
        f(args...)
    end
    return ntuple(_ -> @allocated(f(args...)), 3)
end
function shared_kernel_calls(T, n, permuted)
    D = length(n)
    m = SpectralOperators._fourier_metadata(n, ntuple(d -> T(2π + d), D))
    ranges = ntuple(d -> max(1, n[d]÷2):n[d], D)
    a = [
        Complex{T}(T(sin(sum(Tuple(I)))), T(cos(prod(Tuple(I))))) for
        I in CartesianIndices(n)
    ]
    x = fourier_tile(a, ranges; permuted)
    y = fourier_tile(a, ranges; permuted)
    z = fourier_tile(a, ranges; permuted)
    realfield = fourier_tile(real.(a), ranges; permuted)
    SO = SpectralOperators
    return (
        (SO._apply_ik_store!, (x, y, m.ik[D], D)),
        (SO._apply_ik_store_scaled!, (x, y, m.ik[D], D, T(-2))),
        (SO._apply_ik_accum!, (x, y, m.ik[D], D, T(-2))),
        (SO._curl_component_store!, (x, y, m.ik[1], 1, z, m.ik[D], D)),
        (SO._apply_laplacian_multiplier!, (x, m.kfull)),
        (SO._project_fourier_modes!, ((x, y, z), m.kvec)),
        (SO._apply_exp_filter_multiplier!, (x, m.midx, n, T(5), 3)),
        (SO._apply_two_thirds_mask!, (x, m.midx, n)),
        (SO._laplacian_gain, (m.kfull,)),
        (SO._copy_fft_input!, (x, realfield, T(2))),
        (SO._store_fft_output!, (realfield, x, T(2))),
        (SO._copy_fft_input!, (x, realfield, one(T))),
        (SO._store_fft_output!, (realfield, x, one(T))),
        (SO._fft_safe_magnitude, (T, prod(n), one(T))),
        (SO._fft_input_scale, (floatmax(T) / T(8), floatmax(T) / T(16))),
    )
end
@testset "Shared global-view kernels allocate nothing" begin
    for T in (Float32, Float64), n in ((8,), (6, 5), (4, 3, 2)), permuted in (false, true)
        for (f, args) in shared_kernel_calls(T, n, permuted)
            @test all(iszero, shared_kernel_allocations(f, args))
        end
    end
end

@testset "Shared kernel boundary controls" begin
    SO = SpectralOperators
    for T in (Float32, Float64)
        m = SO._fourier_metadata((6,), (T(6),))
        a = Complex{T}.(1:6, 6:-1:1)
        b = copy(a)
        expected = a .* m.ik[1]
        alias = copy(a)
        @test SO._curl_component_store!(alias, alias, m.ik[1], 1, b, m.ik[1], 1) === alias
        @test iszero(alias)
        # A local-index multiplier cannot equal this tile's global-index result.
        tile = fourier_tile(a, (4:6,))
        out = fourier_tile(b, (4:6,))
        SO._apply_ik_store!(out, tile, m.ik[1], 1)
        @test tile_equals(out, expected)
        wrong = [a[j+3] * m.ik[1][j] for j = 1:3]
        @test any(j -> !isequal(out[j+3], wrong[j]), 1:3)
        # Strict cutoff: a non-strict rule wrongly retains modes ±N/3.
        masked = copy(a)
        SO._apply_two_thirds_mask!(masked, m.midx, m.n)
        @test iszero(masked[3]) && iszero(masked[5])
        @test !iszero(a[3]) && !iszero(a[5])
        # Full Nyquist belongs to Laplacian, but not to the first derivative.
        nyquist = fill(zero(Complex{T}), 6)
        nyquist[4] = one(T)
        lap = copy(nyquist)
        odd = copy(nyquist)
        SO._apply_laplacian_multiplier!(lap, m.kfull)
        SO._apply_ik_store!(odd, odd, m.ik[1], 1)
        @test real(lap[4]) < zero(T) && iszero(imag(lap[4]))
        @test iszero(odd)
        # Every mixed DC/Nyquist corner is an exact projection null mode.
        m2 = SO._fourier_metadata((4, 6), (T(4), T(6)))
        hats = ntuple(c -> fill(Complex{T}(c, -c), 4, 6), 3)
        old = map(copy, hats)
        SO._project_fourier_modes!(hats, m2.kvec)
        for i in (1, 3), j in (1, 4), c = 1:3
            @test isequal(hats[c][i, j], old[c][i, j])
        end
        saved = map(copy, hats)
        @test_throws DimensionMismatch SO._project_fourier_modes!(hats, m.kvec)
        @test isequal(hats, saved)
        # Copy/store operate on global logical axes, including permuted storage.
        raw = T.(reshape(1:24, 4, 3, 2))
        ranges = (2:4, 2:3, 1:2)
        realtile = fourier_tile(raw, ranges; permuted = true)
        complex = fourier_tile(Complex{T}.(raw), ranges; permuted = true)
        SO._copy_fft_input!(complex, realtile, T(2))
        @test all(I -> complex[I] == raw[I] / T(2), CartesianIndices(complex))
        SO._store_fft_output!(realtile, complex, T(2))
        @test tile_equals(realtile, raw)
        SO._copy_fft_input!(complex, realtile, one(T))
        @test all(I -> complex[I] == raw[I], CartesianIndices(complex))
        fill!(realtile, zero(T))
        SO._store_fft_output!(realtile, complex, one(T))
        @test tile_equals(realtile, raw)
        extreme = T[floatmax(T)/T(2), nextfloat(zero(T))]
        scratch = zeros(Complex{T}, 2)
        @test_throws ArgumentError SO._copy_fft_input!(scratch, extreme, T(2))
        @test SO._fft_safe_magnitude(T, 24, one(T)) < SO._fft_safe_magnitude(T, 12, one(T))
    end
end
