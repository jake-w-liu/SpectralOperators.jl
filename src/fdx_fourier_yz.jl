# sbp.jl — mixed shock operator (§8.2): a non-periodic summation-by-parts (SBP)
# first derivative along the shock normal x, Fourier along the periodic
# transverse directions. Avoids imposing periodicity between upstream and
# downstream plasma states while keeping spectral accuracy on the shock surface.
#
# 2nd-order SBP-(2,1) operator (Strand 1994): central interior, one-sided
# boundaries, diagonal norm H = dx·diag(½,1,…,1,½). It satisfies
#   uᵀ H D v + (Du)ᵀ H v = u_N v_N − u_1 v_1        (summation by parts, OP-005),
# i.e. H D + (H D)ᵀ = diag(−1,0,…,0,1).
# The current implementation is the 2nd-order SBP-(2,1) operator used to verify
# the mixed scheme; higher-order SBP closures can be added as separate types.

"2nd-order SBP first-derivative operator on `n` nodes over [0,L] (endpoints included)."
struct SBP1D{T}
    n::Int
    dx::T
    H::Vector{T}        # diagonal quadrature weights (the SBP norm)
end

function SBP1D(n::Integer, L::T) where {T<:AbstractFloat}
    n >= 3 || throw(ArgumentError("SBP needs at least 3 nodes"))
    L > 0 || throw(ArgumentError("SBP domain length must be positive"))
    dx = L / (n - 1)
    H = fill(dx, n)
    H[1] = dx / 2
    H[n] = dx / 2
    return SBP1D{T}(n, dx, H)
end

"In-place SBP first derivative of a 1-D array."
function sbp_deriv!(out::AbstractVector{T}, f::AbstractVector{T}, s::SBP1D{T}) where {T}
    n = s.n
    dx = s.dx
    length(f) == n || throw(DimensionMismatch("length $(length(f)) ≠ $n"))
    length(out) == n || throw(DimensionMismatch("output length $(length(out)) ≠ $n"))
    Base.mightalias(out, f) && throw(ArgumentError("sbp_deriv! output must not alias input"))
    @inbounds begin
        out[1] = (f[2] - f[1]) / dx                      # one-sided
        for i = 2:n-1
            out[i] = (f[i+1] - f[i-1]) / (2dx)           # central
        end
        out[n] = (f[n] - f[n-1]) / dx                    # one-sided
    end
    return out
end

sbp_deriv(f::AbstractVector{T}, s::SBP1D{T}) where {T} = sbp_deriv!(similar(f), f, s)

"Apply the SBP derivative along x (dim 1) of a 2-D array, column by column."
function sbp_deriv_x!(out::AbstractMatrix{T}, f::AbstractMatrix{T}, s::SBP1D{T}) where {T}
    size(f, 1) == s.n || throw(DimensionMismatch("input first dimension $(size(f, 1)) ≠ $(s.n)"))
    size(out) == size(f) || throw(DimensionMismatch("output size $(size(out)) does not match input size $(size(f))"))
    Base.mightalias(out, f) && throw(ArgumentError("sbp_deriv_x! output must not alias input"))
    @inbounds for j in axes(f, 2)
        sbp_deriv!(view(out, :, j), view(f, :, j), s)
    end
    return out
end

"""
    fourier_deriv_y!(out, f, Ly)

Spectral derivative along the periodic transverse direction y (dim 2) of a 2-D
field, with the Nyquist mode zeroed (odd derivative).
"""
function fourier_deriv_y!(out::Matrix{T}, f::Matrix{T}, Ly::T) where {T}
    work = FourierDerivYWorkspace(size(f, 1), size(f, 2), Ly)
    return fourier_deriv_y!(out, f, work)
end

"""
    FourierDerivYWorkspace(nx, ny, Ly)
    FourierDerivYWorkspace(f, Ly)

Reusable FFT workspace for allocation-free [`fourier_deriv_y!`](@ref) calls on
`nx × ny` matrices with periodic length `Ly` along the second dimension.
"""
struct FourierDerivYWorkspace{T,P,PI}
    nx::Int
    ny::Int
    Ly::T
    ky::Vector{T}
    cbuf::Matrix{Complex{T}}
    plan::P
    iplan::PI
end

function FourierDerivYWorkspace(nx::Integer, ny::Integer, Ly::T) where {T<:AbstractFloat}
    nx >= 1 || throw(ArgumentError("nx must be positive"))
    ny >= 1 || throw(ArgumentError("ny must be positive"))
    Ly > 0 || throw(ArgumentError("Ly must be positive"))
    nxi = Int(nx)
    nyi = Int(ny)
    ky = Vector{T}(undef, nyi)
    for m = 0:nyi-1
        mp = m <= nyi ÷ 2 ? m : m - nyi
        ky[m+1] = T(2π) * mp / Ly
    end
    iseven(nyi) && (ky[nyi÷2+1] = zero(T))
    cbuf = zeros(Complex{T}, nxi, nyi)
    plan = plan_fft!(cbuf, 2)
    iplan = plan_ifft!(cbuf, 2)
    return FourierDerivYWorkspace{T,typeof(plan),typeof(iplan)}(nxi, nyi, Ly, ky, cbuf, plan, iplan)
end

function FourierDerivYWorkspace(f::AbstractMatrix{T}, Ly::T) where {T<:AbstractFloat}
    return FourierDerivYWorkspace(size(f, 1), size(f, 2), Ly)
end

function fourier_deriv_y!(out::Matrix{T}, f::Matrix{T}, work::FourierDerivYWorkspace{T}) where {T}
    size(f) == (work.nx, work.ny) ||
        throw(DimensionMismatch("input size $(size(f)) does not match workspace size $((work.nx, work.ny))"))
    size(out) == size(f) ||
        throw(DimensionMismatch("output size $(size(out)) does not match input size $(size(f))"))
    work.cbuf .= f
    work.plan * work.cbuf
    @inbounds for m = 1:work.ny
        ik = Complex{T}(zero(T), work.ky[m])
        for i = 1:work.nx
            work.cbuf[i, m] *= ik
        end
    end
    work.iplan * work.cbuf
    out .= real.(work.cbuf)
    return out
end
