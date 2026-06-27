# differential_operators.jl — periodic Fourier ∂/∇/∇·/∇×/∇².

@inline function _apply_ik_store!(
    dst::AbstractArray{Complex{T},D},
    fhat::AbstractArray{Complex{T},D},
    ikj::Vector{Complex{T}},
    j::Int,
) where {D,T}
    @inbounds for I in CartesianIndices(dst)
        dst[I] = fhat[I] * ikj[I[j]]
    end
    return dst
end

# acc[I] += s * f̂[I] * (i k_j)
@inline function _apply_ik_accum!(
    acc::AbstractArray{Complex{T},D},
    fhat::AbstractArray{Complex{T},D},
    ikj::Vector{Complex{T}},
    j::Int,
    s::T,
) where {D,T}
    @inbounds for I in CartesianIndices(acc)
        acc[I] += s * fhat[I] * ikj[I[j]]
    end
    return acc
end

# acc .+= s * ∂_j(field). cbuf is scratch.
@inline function _accum_deriv!(
    acc::AbstractArray{Complex{T},D},
    field::AbstractArray{T,D},
    g::FourierGrid{D,T},
    j::Int,
    s::T,
) where {D,T}
    g.cbuf .= field
    g.plan * g.cbuf
    _apply_ik_accum!(acc, g.cbuf, g.ik[j], j, s)
    return acc
end

"""
    deriv!(out, f, g, j)

In-place first derivative of scalar field `f` along axis `j` (1≤j≤D). `out` and
`f` may alias.
"""
function deriv!(
    out::AbstractArray{T,D},
    f::AbstractArray{T,D},
    g::FourierGrid{D,T},
    j::Int,
) where {D,T}
    1 <= j <= D || throw(ArgumentError("axis $j out of range 1:$D"))
    _require_grid_array(:input, f, g)
    _require_grid_array(:output, out, g)
    g.cbuf .= f
    g.plan * g.cbuf
    _apply_ik_store!(g.tbuf, g.cbuf, g.ik[j], j)
    g.iplan * g.tbuf
    out .= real.(g.tbuf)
    return out
end

"Allocating convenience wrapper for [`deriv!`](@ref)."
deriv(f::AbstractArray{T,D}, g::FourierGrid{D,T}, j::Int) where {D,T} = deriv!(similar(f), f, g, j)

"""
    gradient!(out::NTuple{D}, f, g)

Gradient of scalar field `f`; `out[j]` receives ∂_j f. One forward FFT, D inverses.
"""
function gradient!(
    out::Tuple{Vararg{AbstractArray{T,D},D}},
    f::AbstractArray{T,D},
    g::FourierGrid{D,T},
) where {D,T}
    _require_grid_array(:input, f, g)
    for j = 1:D
        _require_grid_array(:output, j, out[j], g)
    end
    for j = 1:D-1, k = j+1:D
        Base.mightalias(out[j], out[k]) &&
            throw(ArgumentError("gradient! output components must not alias each other"))
    end
    g.cbuf .= f
    g.plan * g.cbuf                      # cbuf = f̂, reused for every axis
    for j = 1:D
        _apply_ik_store!(g.tbuf, g.cbuf, g.ik[j], j)
        g.iplan * g.tbuf
        out[j] .= real.(g.tbuf)
    end
    return out
end

"""
    divergence!(out, v::NTuple{N}, g)

Divergence Σ_{j=1}^{D} ∂_j v[j] of a vector field. Only the first `D` (spatial)
components are used, so a 3-component velocity field works in any dimension.
"""
function divergence!(out::AbstractArray{T,D}, v, g::FourierGrid{D,T}) where {D,T}
    length(v) >= D || throw(ArgumentError("divergence! needs at least $D vector components"))
    _require_grid_array(:output, out, g)
    for j = 1:D
        _require_grid_array(:input, j, v[j], g)
    end
    fill!(g.abuf, zero(Complex{T}))
    for j = 1:D
        _accum_deriv!(g.abuf, v[j], g, j, one(T))
    end
    g.iplan * g.abuf
    out .= real.(g.abuf)
    return out
end

"""
    curl!(out::NTuple{3}, A::NTuple{3}, g)

Curl of a 3-component vector field `A` whose components vary over the `D`
spatial axes (derivatives along absent axes are zero — the standard reduced curl
for 1D/2D). `out` must be distinct from `A`.
"""
function curl!(
    out::Tuple{Vararg{AbstractArray{T,D},3}},
    A::Tuple{Vararg{AbstractArray{T,D},3}},
    g::FourierGrid{D,T},
) where {D,T}
    D <= 3 || throw(ArgumentError("curl! supports spatial dimension D ≤ 3; got $D"))
    for c = 1:3
        _require_grid_array(:output, c, out[c], g)
        _require_grid_array(:input, c, A[c], g)
    end
    for c = 1:2, d = c+1:3
        Base.mightalias(out[c], out[d]) &&
            throw(ArgumentError("curl! output components must not alias each other"))
    end
    for c = 1:3, d = 1:3
        Base.mightalias(out[c], A[d]) &&
            throw(ArgumentError("curl! output components must not alias input components"))
    end
    o = one(T)
    m = -one(T)
    # (curl A)_x = ∂_y A_z − ∂_z A_y
    fill!(g.abuf, zero(Complex{T}))
    D >= 2 && _accum_deriv!(g.abuf, A[3], g, 2, o)
    D >= 3 && _accum_deriv!(g.abuf, A[2], g, 3, m)
    g.iplan * g.abuf
    out[1] .= real.(g.abuf)
    # (curl A)_y = ∂_z A_x − ∂_x A_z
    fill!(g.abuf, zero(Complex{T}))
    D >= 3 && _accum_deriv!(g.abuf, A[1], g, 3, o)
    _accum_deriv!(g.abuf, A[3], g, 1, m)
    g.iplan * g.abuf
    out[2] .= real.(g.abuf)
    # (curl A)_z = ∂_x A_y − ∂_y A_x
    fill!(g.abuf, zero(Complex{T}))
    _accum_deriv!(g.abuf, A[2], g, 1, o)
    D >= 2 && _accum_deriv!(g.abuf, A[1], g, 2, m)
    g.iplan * g.abuf
    out[3] .= real.(g.abuf)
    return out
end

"""
    laplacian!(out, f, g)

Spectral Laplacian ∇²f = −|k|² f̂ (periodic). Used for hyperresistivity and
diffusion terms.
"""
function laplacian!(out::AbstractArray{T,D}, f::AbstractArray{T,D}, g::FourierGrid{D,T}) where {D,T}
    _require_grid_array(:input, f, g)
    _require_grid_array(:output, out, g)
    g.cbuf .= f
    g.plan * g.cbuf
    @inbounds for I in CartesianIndices(g.cbuf)
        k2 = zero(T)
        for d = 1:D
            kk = g.kfull[d][I[d]]
            k2 += kk * kk
        end
        g.cbuf[I] *= -k2
    end
    g.iplan * g.cbuf
    out .= real.(g.cbuf)
    return out
end
