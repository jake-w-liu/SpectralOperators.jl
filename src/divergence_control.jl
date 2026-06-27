# divergence_control.jl — divergence-free projection.

"""
    project_divfree!(B::NTuple{3}, g)

Project a 3-component field onto its divergence-free part in place:
B̂ ← (I − k kᵀ/k²) B̂ for k≠0, leaving the k=0 mean untouched. After projection
the discrete divergence (same wavenumbers) is zero to roundoff.
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
    Bx .= B[1]
    g.plan * Bx
    if D >= 2
        By = g.tbuf
        By .= B[2]
        g.plan * By
    end
    if D >= 3
        Bz = g.abuf
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
        k2 = wx * wx + wy * wy + wz * wz
        if k2 > 0
            f = wx * Bx[I]
            D >= 2 && (f += wy * By[I])
            D >= 3 && (f += wz * Bz[I])
            f /= k2
            Bx[I] -= wx * f
            D >= 2 && (By[I] -= wy * f)
            D >= 3 && (Bz[I] -= wz * f)
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
