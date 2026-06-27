# filters.jl — spectral filters, dealiasing, FFT planning, binomial smoothing

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

    σ(k) = exp( -α · max_d(|k_d| / k_max_d)^(2p) ),

where `k_max_d = 2π·(n_d÷2)/L_d`. Low-|k| modes are left essentially unchanged
(σ → 1), modes near the Nyquist boundary on any axis are strongly damped, and
the `k = 0` mean is preserved exactly (σ(0)=1).

Larger `α` damps harder; larger `p` makes the filter sharper (closer to a brick
wall). Defaults `α=36, p=8` give σ ≈ e^{-36} ≈ 2.3e-16 for modes that reach the
Nyquist boundary on any axis, i.e. machine-zero damping of the highest modes
while barely touching low modes.
Returns `field`.
"""
function exp_filter!(
    field::AbstractArray{T,D},
    g::FourierGrid{D,T};
    α::Real = 36,
    p::Integer = 8,
) where {T,D}
    _require_grid_array(:field, field, g)
    p >= 1 || throw(ArgumentError("p must be ≥ 1"))
    α >= 0 || throw(ArgumentError("α must be ≥ 0"))
    αT = T(α)
    kmax = ntuple(d -> _kmax_axis(g.n[d], g.L[d]), D)
    g.cbuf .= field
    g.plan * g.cbuf
    @inbounds for I in CartesianIndices(g.cbuf)
        ratio2 = zero(T)
        for d = 1:D
            km = kmax[d]
            if km > 0
                kk = g.kfull[d][I[d]]
                ratio2 = max(ratio2, (kk * kk) / (km * km))
            end
        end
        if ratio2 == 0
            # σ(0) = 1 exactly; degenerate all-singleton grids have no damping.
            continue
        end
        σ = exp(-αT * ratio2^Int(p))             # (ratio2)^p = max_d(|k_d|/kmax_d)^{2p}
        g.cbuf[I] *= σ
    end
    g.iplan * g.cbuf
    field .= real.(g.cbuf)
    return field
end

"""
    dealias_two_thirds!(field::Array{T,D}, g::FourierGrid{D,T})

Apply Orszag's 2/3 de-aliasing rule in place: retain only integer Fourier modes
with `3*abs(m_d) < n_d` on every axis, and zero the rest. The strict boundary is
important when `n_d` is divisible by 3: retaining `abs(m_d) == n_d/3` lets the
quadratic self-interaction alias back into the retained band. The `k = 0` mean is
preserved. Returns `field`.
"""
function dealias_two_thirds!(field::AbstractArray{T,D}, g::FourierGrid{D,T}) where {T,D}
    _require_grid_array(:field, field, g)
    g.cbuf .= field
    g.plan * g.cbuf
    @inbounds for I in CartesianIndices(g.cbuf)
        ok = true
        for d = 1:D
            ok &= 3 * abs(g.midx[d][I[d]]) < g.n[d]
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
    n <= typemax(Int) || throw(OverflowError("no Int fft-friendly size can be >= $n"))
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
        c == typemax(Int) && throw(OverflowError("no Int fft-friendly size can be >= $n"))
        c += 1
    end
    return c
end

"""
    with_fftw_wisdom(f, path::AbstractString)

Run `f()` while attempting to persist FFTW plan wisdom at `path`. Before calling
`f`, import any wisdom already saved at `path` (silently skipped if the file is
absent or cannot be read). After `f` returns — even if it throws — attempt to
export the current accumulated wisdom to `path` so newly created plans can be
remembered for next time. Import/export failures are non-fatal. Returns whatever
`f()` returns.

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

# --- binomial smoothing ---
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
    BinomialSmoothWorkspace(g::FourierGrid)

Reusable real line buffer for allocation-free [`binomial_smooth!`](@ref) calls
on fields matching `g`.
"""
struct BinomialSmoothWorkspace{T}
    buf::Vector{T}
end

function BinomialSmoothWorkspace(g::FourierGrid{D,T}) where {D,T}
    return BinomialSmoothWorkspace{T}(Vector{T}(undef, maximum(g.n)))
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
    work = BinomialSmoothWorkspace(g)
    return binomial_smooth!(field, g, work; passes)
end

"""
    binomial_smooth!(field, g, work::BinomialSmoothWorkspace; passes=1)

Workspace-backed variant of [`binomial_smooth!`](@ref). Reuse `work` across
calls on the same real element type to avoid allocating the line buffer.
"""
function binomial_smooth!(
    field::Array{T,D},
    g::FourierGrid{D,T},
    work::BinomialSmoothWorkspace{T};
    passes::Int = 1,
) where {T,D}
    passes >= 0 || throw(ArgumentError("passes must be non-negative, got $passes"))
    size(field) == g.n ||
        throw(DimensionMismatch("field size $(size(field)) does not match grid size $(g.n)"))
    passes == 0 && return field
    buf = work.buf
    length(buf) >= maximum(g.n) ||
        throw(DimensionMismatch("workspace length $(length(buf)) is smaller than required $(maximum(g.n))"))
    for j = 1:D
        for _ = 1:passes
            _binomial_pass_axis!(field, j, buf)
        end
    end
    return field
end
