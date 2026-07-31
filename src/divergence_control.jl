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
    input_scale = _fft_input_scale(B, D, g)
    Bx = g.cbuf
    By = g.tbuf
    Bz = g.abuf
    _copy_fft_input!(Bx, B[1], input_scale)
    g.plan * Bx
    if D >= 2
        _copy_fft_input!(By, B[2], input_scale)
        g.plan * By
    end
    if D >= 3
        _copy_fft_input!(Bz, B[3], input_scale)
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
            # Normalize the field too: the dot product can overflow for finite
            # high-amplitude components even when the projection is finite.
            # Real/imaginary component scaling remains finite in cases where
            # the complex magnitude itself would overflow.
            bx = Bx[I]
            by = D >= 2 ? By[I] : zero(Complex{T})
            bz = D >= 3 ? Bz[I] : zero(Complex{T})
            bscale = max(abs(real(bx)), abs(imag(bx)))
            D >= 2 && (bscale = max(bscale, abs(real(by)), abs(imag(by))))
            D >= 3 && (bscale = max(bscale, abs(real(bz)), abs(imag(bz))))
            iszero(bscale) && continue
            tx = bx / bscale
            ty = by / bscale
            tz = bz / bscale
            f = sx * tx
            D >= 2 && (f += sy * ty)
            D >= 3 && (f += sz * tz)
            f /= s2
            Bx[I] = bscale * (tx - sx * f)
            D >= 2 && (By[I] = bscale * (ty - sy * f))
            D >= 3 && (Bz[I] = bscale * (tz - sz * f))
        end
    end
    g.iplan * Bx
    _store_fft_output!(B[1], Bx, input_scale)
    if D >= 2
        g.iplan * By
        _store_fft_output!(B[2], By, input_scale)
    end
    if D >= 3
        g.iplan * Bz
        _store_fft_output!(B[3], Bz, input_scale)
    end
    return B
end
