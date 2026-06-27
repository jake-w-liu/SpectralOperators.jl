# periodic_fourier.jl — periodic Fourier spectral operators on a collocated grid.
#
# Single primitive: the first-derivative spectral multiplier  ∂_j ↔ i k_j , with
# the Nyquist mode zeroed (the correct choice for an odd-order derivative on an
# even grid: its sine partner is unrepresentable, so a real first derivative must
# drop it). gradient / divergence / curl / divergence-free projection are all
# built from that one multiplier, which makes the discrete identities
# div(curl)=0 and ⟨f,Dg⟩=−⟨Df,g⟩ hold to roundoff (the multipliers commute and
# are used consistently on both sides).
#
# Operators reuse preplanned in-place FFTs and complex workspaces, so they
# allocate nothing in steady state. The workspaces live in the grid object, so a
# FourierGrid is NOT thread-safe: give each thread / rank its own.
# Uses full-complex FFTs for a simple, backend-independent implementation. A
# future rFFT specialization can reduce memory for very large real-valued grids.

"""
    FourierGrid(n::NTuple{D,Int}, L::NTuple{D,T})

Collocated periodic Fourier grid in `D` dimensions with `n[d]` points and
physical length `L[d]` per axis. Precomputes wavenumbers, in-place FFT plans,
and complex scratch buffers.
"""
struct FourierGrid{D,T,P,PI}
    n::NTuple{D,Int}
    L::NTuple{D,T}
    dx::NTuple{D,T}
    midx::NTuple{D,Vector{Int}}        # signed integer Fourier mode per axis
    kfull::NTuple{D,Vector{T}}         # real k per axis, Nyquist included (even derivatives / filters)
    ik::NTuple{D,Vector{Complex{T}}}   # i*k per axis, Nyquist zeroed (1st-deriv multiplier)
    kvec::NTuple{D,Vector{T}}          # real k per axis, Nyquist zeroed (projection geometry)
    plan::P                            # in-place forward FFT over all dims
    iplan::PI                          # in-place inverse FFT over all dims
    cbuf::Array{Complex{T},D}          # scratch / forward transform of input
    tbuf::Array{Complex{T},D}          # scratch for per-axis multiply+inverse
    abuf::Array{Complex{T},D}          # complex accumulator
end

function FourierGrid(n::Tuple{Int,Vararg{Int}}, L::Tuple{T,Vararg{T}}) where {T<:AbstractFloat}
    D = length(n)
    length(L) == D || throw(DimensionMismatch("domain length tuple has length $(length(L)); expected $D"))
    all(>(0), n) || throw(ArgumentError("grid sizes must be positive"))
    all(l -> isfinite(l) && l > zero(T), L) ||
        throw(ArgumentError("domain lengths must be positive and finite"))
    dx = ntuple(d -> L[d] / n[d], D)
    midx = ntuple(D) do d
        N = n[d]
        modes = Vector{Int}(undef, N)
        for m = 0:N-1
            mp = m <= N ÷ 2 ? m : m - N          # signed mode index
            modes[m+1] = mp
        end
        modes
    end
    kfull = ntuple(D) do d
        N = n[d]
        k = Vector{T}(undef, N)
        scale = T(2π) / L[d]
        @inbounds for i = 1:N
            k[i] = scale * midx[d][i]
        end
        k
    end
    kvec = ntuple(D) do d
        N = n[d]
        k = copy(kfull[d])
        if iseven(N)
            k[N÷2+1] = zero(T)               # zero Nyquist for first derivative
        end
        k
    end
    ik = ntuple(d -> Complex{T}.(zero(T), kvec[d]), D)
    cbuf = zeros(Complex{T}, n)
    tbuf = zeros(Complex{T}, n)
    abuf = zeros(Complex{T}, n)
    plan = plan_fft!(cbuf)
    iplan = plan_ifft!(cbuf)
    FourierGrid{D,T,typeof(plan),typeof(iplan)}(n, L, dx, midx, kfull, ik, kvec, plan, iplan, cbuf, tbuf, abuf)
end

@inline function _require_grid_array(name::Symbol, a::AbstractArray, g::FourierGrid)
    size(a) == g.n ||
        throw(DimensionMismatch("$(name) size $(size(a)) does not match grid size $(g.n)"))
    axes(a) == axes(g.cbuf) ||
        throw(DimensionMismatch("$(name) axes $(axes(a)) must be one-based axes $(axes(g.cbuf))"))
    return nothing
end

@inline function _require_grid_array(name::Symbol, c::Int, a::AbstractArray, g::FourierGrid)
    size(a) == g.n ||
        throw(DimensionMismatch("$(name)[$c] size $(size(a)) does not match grid size $(g.n)"))
    axes(a) == axes(g.cbuf) ||
        throw(DimensionMismatch("$(name)[$c] axes $(axes(a)) must be one-based axes $(axes(g.cbuf))"))
    return nothing
end

# dst[I] = f̂[I] * (i k_j)   — apply the axis-j first-derivative multiplier.
# Explicit Cartesian loop (no reshape) so it is allocation-free for runtime j.
