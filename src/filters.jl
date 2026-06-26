# Filters.jl — spectral filters, dealiasing, FFT planning, binomial smoothing

# True (Nyquist-inclusive) signed mode wavenumber for axis of length N, domain L.
@inline function _mode_wavenumber(m::Int, N::Int, L::T) where {T<:AbstractFloat}
    mp = m <= N ÷ 2 ? m : m - N     # signed mode index in (−N/2, N/2]
    return T(2π) * mp / L
end

# Per-axis maximum representable |k| = 2π·(N÷2)/L.
@inline _kmax_axis(N::Int, L::T) where {T<:AbstractFloat} = T(2π) * (N ÷ 2) / L

"""
    exp_filter!(field::Array{T,D}, g::FourierGrid{D,T}; α=36, p=8)

Apply an exponential spectral low-pass filter to `field` in place. Each Fourier
mode is multiplied by the transfer function

    σ(k) = exp( -α · (|k| / k_max)^(2p) ),

where `|k| = sqrt(Σ_d k_d²)` is the mode's wavenumber magnitude and
`k_max = sqrt(Σ_d k_max_d²)` with `k_max_d = 2π·(n_d÷2)/L_d`. Low-|k| modes are
left essentially unchanged (σ → 1), modes near the Nyquist boundary are strongly
damped, and the `k = 0` mean is preserved exactly (σ(0)=1).

Larger `α` damps harder; larger `p` makes the filter sharper (closer to a brick
wall). Defaults `α=36, p=8` give σ ≈ e^{-36} ≈ 2.3e-16 at the corner of k-space,
i.e. machine-zero damping of the highest mode while barely touching low modes.
Returns `field`.
"""
function exp_filter!(
    field::AbstractArray{T,D},
    g::FourierGrid{D,T};
    α::Real = 36,
    p::Integer = 8,
) where {T,D}
    size(field) == g.n || throw(DimensionMismatch("field size $(size(field)) ≠ grid $(g.n)"))
    p >= 1 || throw(ArgumentError("p must be ≥ 1"))
    α >= 0 || throw(ArgumentError("α must be ≥ 0"))
    αT = T(α)
    # k_max² = Σ_d k_max_d²  (squared corner wavenumber). All axes have n>0.
    kmax2 = zero(T)
    @inbounds for d = 1:D
        km = _kmax_axis(g.n[d], g.L[d])
        kmax2 += km * km
    end
    # Precompute per-axis squared mode wavenumbers (Nyquist included).
    k2axis = ntuple(D) do d
        N = g.n[d]
        v = Vector{T}(undef, N)
        @inbounds for m = 0:N-1
            kk = _mode_wavenumber(m, N, g.L[d])
            v[m+1] = kk * kk
        end
        v
    end
    g.cbuf .= field
    g.plan * g.cbuf
    @inbounds for I in CartesianIndices(g.cbuf)
        t = Tuple(I)
        k2 = zero(T)
        for d = 1:D
            k2 += k2axis[d][t[d]]
        end
        if k2 == 0 || kmax2 == 0
            # σ(0) = 1 exactly; degenerate kmax2==0 (n==1 everywhere) ⇒ no damping.
            continue
        end
        ratio2 = k2 / kmax2                      # (|k|/k_max)²  ∈ (0,1]
        σ = exp(-αT * ratio2^Int(p))             # (ratio2)^p = (|k|/k_max)^{2p}
        g.cbuf[I] *= σ
    end
    g.iplan * g.cbuf
    field .= real.(g.cbuf)
    return field
end

"""
    dealias_two_thirds!(field::Array{T,D}, g::FourierGrid{D,T})

Apply Orszag's 2/3 de-aliasing rule in place: zero every Fourier mode whose
per-axis wavenumber exceeds `(2/3)·k_max_d` on ANY axis, where
`k_max_d = 2π·(n_d÷2)/L_d`. This removes the highest one-third of modes on each
axis so that a single quadratic nonlinearity cannot alias energy back into the
retained band. The `k = 0` mean is preserved (it is well below the cutoff).
Returns `field`.
"""
function dealias_two_thirds!(field::AbstractArray{T,D}, g::FourierGrid{D,T}) where {T,D}
    size(field) == g.n || throw(DimensionMismatch("field size $(size(field)) ≠ grid $(g.n)"))
    # Per-axis keep-mask: true if |k_d| ≤ (2/3) k_max_d.
    keep = ntuple(D) do d
        N = g.n[d]
        cut = T(2) / T(3) * _kmax_axis(N, g.L[d])
        m = Vector{Bool}(undef, N)
        @inbounds for i = 0:N-1
            m[i+1] = abs(_mode_wavenumber(i, N, g.L[d])) <= cut
        end
        m
    end
    g.cbuf .= field
    g.plan * g.cbuf
    @inbounds for I in CartesianIndices(g.cbuf)
        t = Tuple(I)
        ok = true
        for d = 1:D
            ok &= keep[d][t[d]]
        end
        ok || (g.cbuf[I] = zero(Complex{T}))
    end
    g.iplan * g.cbuf
    field .= real.(g.cbuf)
    return field
end

"""
    fft_friendly_size(n::Integer) -> Int

Return the smallest integer ≥ `n` that is 7-smooth, i.e. a product of only the
small primes 2, 3, 5 and 7 (the radices FFTW handles most efficiently). For
`n ≤ 1` returns 1. Examples: `fft_friendly_size(17) == 18` (=2·3²),
`fft_friendly_size(100) == 100`, `fft_friendly_size(101) == 105` (=3·5·7).
"""
function fft_friendly_size(n::Integer)
    n <= 1 && return 1
    is7smooth(m::Int) = begin
        for p in (2, 3, 5, 7)
            while m % p == 0
                m ÷= p
            end
        end
        m == 1
    end
    c = Int(n)
    while !is7smooth(c)
        c += 1
    end
    return c
end

"""
    with_fftw_wisdom(f, path::AbstractString)

Run `f()` with FFTW plan wisdom persisted at `path`. Before calling `f`, import
any wisdom already saved at `path` (silently skipped if the file is absent or
cannot be read). After `f` returns — even if it throws — export the current
accumulated wisdom to `path` so newly created plans are remembered for next time.
Returns whatever `f()` returns.

Wisdom is FFTW's record of empirically tuned plans; reusing it lets repeated
`plan_fft!`/`plan_ifft!` calls (e.g. across runs or many `FourierGrid`s) reuse
prior planning effort instead of re-measuring.
"""
function with_fftw_wisdom(f, path::AbstractString)
    if isfile(path)
        try
            FFTW.import_wisdom(path)
        catch
            # Corrupt/incompatible wisdom file: ignore and plan from scratch.
        end
    end
    try
        return f()
    finally
        try
            dir = dirname(path)
            isempty(dir) || isdir(dir) || mkpath(dir)
            FFTW.export_wisdom(path)
        catch
            # Non-fatal: failure to persist wisdom must not break the computation.
        end
    end
end

# --- binomial smoothing (from smoothing.jl) ---
"""
    smoothing_transfer(k, dx; passes::Int=1)

Analytic Fourier transfer function of the binomial `(1,2,1)/4` smoothing stencil
for a single axis, for a mode of angular wavenumber `k` on a grid of spacing
`dx`, applied `passes` times.

One pass of the stencil maps `f[i] ← (f[i-1] + 2 f[i] + f[i+1]) / 4`. Acting on
a grid sample of `cos(k x)` (or any pure mode `exp(i k x)`) this multiplies the
amplitude by `(1 + cos(k dx)) / 2 = cos²(k dx / 2)`. `passes` repeated
applications multiply by that factor raised to `passes`.

The factor is `1` at `k = 0` (uniform fields, hence total/mean content, are
preserved exactly) and `0` at the grid Nyquist `k dx = π` (the (-1)^i mode is
fully removed). It is always in `[0, 1]`, so the filter never amplifies.
"""
function smoothing_transfer(k, dx; passes::Int = 1)
    passes >= 0 || throw(ArgumentError("passes must be non-negative, got $passes"))
    c = cos(k * dx / 2)^2                       # per-pass per-axis response cos²(k dx/2)
    return c^passes
end

# Single in-place pass of the (1,2,1)/4 stencil along axis `j` with periodic
# wrap. `buf` is a length-n[j] scratch vector holding one line at a time, so the
# update reads original neighbour values (no in-place contamination along the
# line). All other axes are looped over via CartesianIndices of the complement.
function _binomial_pass_axis!(field::Array{T,D}, j::Int, buf::Vector{T}) where {T,D}
    N = size(field, j)
    if N < 3
        # With fewer than 3 points along an axis the periodic 3-point stencil
        # folds onto itself; the constant-preserving limit is the identity.
        return field
    end
    q = T(1) / T(4)
    half = T(2) * q                      # = 1/2, weight on the centre (2/4)
    # Iterate over every line parallel to axis j. The set of lines is indexed by
    # the Cartesian product of the other axes' index ranges.
    other = ntuple(d -> d == j ? Base.OneTo(1) : axes(field, d), D)
    @inbounds for base in CartesianIndices(other)
        bt = Tuple(base)
        # Copy the line into buf.
        for i = 1:N
            idx = ntuple(d -> d == j ? i : bt[d], D)
            buf[i] = field[CartesianIndex(idx)]
        end
        # Apply (1,2,1)/4 with periodic neighbours, writing back.
        for i = 1:N
            ip = i == N ? 1 : i + 1
            im = i == 1 ? N : i - 1
            val = q * buf[im] + half * buf[i] + q * buf[ip]
            idx = ntuple(d -> d == j ? i : bt[d], D)
            field[CartesianIndex(idx)] = val
        end
    end
    return field
end

"""
    binomial_smooth!(field::Array{T,D}, g::FourierGrid{D,T}; passes::Int=1)

Apply the binomial `(1,2,1)/4` smoothing stencil in place to `field`, `passes`
times along **each** axis, with periodic wrap at the boundaries (matching the
periodic Fourier grid `g`). Returns `field`.

This is an explicit, optional moment filter: it is never applied automatically.
Each axis-pass multiplies the amplitude of a mode `k` by
[`smoothing_transfer`](@ref)`(k, g.dx[axis]; passes=1)`; the full effect of
`passes` passes on all axes is the product over axes of
`smoothing_transfer(k_axis, g.dx[axis]; passes)`. The `k = 0` (uniform / mean /
total) content is preserved exactly, so smoothing a deposited density conserves
total charge.

`g` supplies the dimensionality and is the natural place the grid spacing lives;
`field` must have the grid's shape.
"""
function binomial_smooth!(field::Array{T,D}, g::FourierGrid{D,T}; passes::Int = 1) where {T,D}
    passes >= 0 || throw(ArgumentError("passes must be non-negative, got $passes"))
    size(field) == g.n ||
        throw(DimensionMismatch("field size $(size(field)) does not match grid size $(g.n)"))
    passes == 0 && return field
    # One reusable line buffer sized to the largest axis; each axis pass uses
    # only its first n[j] entries.
    buf = Vector{T}(undef, maximum(g.n))
    for j = 1:D
        for _ = 1:passes
            _binomial_pass_axis!(field, j, buf)
        end
    end
    return field
end
