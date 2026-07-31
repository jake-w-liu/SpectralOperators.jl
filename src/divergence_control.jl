# divergence_control.jl — divergence-free projection.

"""
    project_divfree!(B::NTuple{3}, g)

Project a 3-component field onto its divergence-free part in place:
B̂ ← (I − k kᵀ/k²) B̂ for derivative-resolved k≠0, leaving modes whose
Nyquist-zeroed derivative wavenumber is zero untouched. After projection the
discrete divergence (same wavenumbers as [`divergence!`](@ref)) is zero to
roundoff.
"""
function project_divfree!(B::Tuple{Vararg{AbstractArray{T,D},3}}, g::FourierGrid{D,T}) where {D,T}
    D <= 3 || throw(ArgumentError("project_divfree! supports spatial dimension D ≤ 3; got $D"))
    for c = 1:3
        _require_grid_array(:component, c, B[c], g)
    end
    for c = 1:2, d = c+1:3
        Base.mightalias(B[c], B[d]) &&
            throw(ArgumentError("project_divfree! components must not alias each other"))
    end
    Bx = g.cbuf
    By = g.tbuf
    Bz = g.abuf
    Bx .= B[1]
    g.plan * Bx
    if D >= 2
        By .= B[2]
        g.plan * By
    end
    if D >= 3
        Bz .= B[3]
        g.plan * Bz
    end
    kx = g.kvec[1]
    ky = D >= 2 ? g.kvec[2] : kx
    kz = D >= 3 ? g.kvec[3] : kx
    @inbounds for I in CartesianIndices(Bx)
        wx = kx[I[1]]
        wy = D >= 2 ? ky[I[2]] : zero(T)
        wz = D >= 3 ? kz[I[3]] : zero(T)
        scale = max(abs(wx), abs(wy), abs(wz))
        if scale > zero(T)
            # Normalize before forming |k|². Squaring a finite, representable
            # Float32 wavenumber can still underflow or overflow on very large or
            # small physical domains; the projection itself is dimensionless.
            sx = wx / scale
            sy = wy / scale
            sz = wz / scale
            s2 = sx * sx + sy * sy + sz * sz
            f = sx * Bx[I]
            D >= 2 && (f += sy * By[I])
            D >= 3 && (f += sz * Bz[I])
            f /= s2
            Bx[I] -= sx * f
            D >= 2 && (By[I] -= sy * f)
            D >= 3 && (Bz[I] -= sz * f)
        end
    end
    g.iplan * Bx
    B[1] .= real.(Bx)
    if D >= 2
        g.iplan * By
        B[2] .= real.(By)
    end
    if D >= 3
        g.iplan * Bz
        B[3] .= real.(Bz)
    end
    return B
end
