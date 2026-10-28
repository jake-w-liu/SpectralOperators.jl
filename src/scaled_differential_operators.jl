# A borrowed Fourier result carries its binary exponent through downstream
# constitutive arithmetic. Public physical output buffers are never repurposed.

function _require_scaled_input(f, g::FourierGrid)
    _require_grid_array(:input, f, g)
    Base.mightalias(f, g.cbuf) &&
        throw(ArgumentError("scaled Fourier input must not alias its mutated scratch"))
    _require_scaled_metadata(g.cbuf, g.kfull)
    _require_scaled_metadata(g.cbuf, g.kvec)
    _require_scaled_metadata(g.cbuf, g.ik)
    _require_scaled_metadata(g.cbuf, g.midx)
    return nothing
end

function _require_scaled_metadata(scratch, vectors)
    for metadata in vectors
        Base.mightalias(scratch, metadata) &&
            throw(ArgumentError("scaled Fourier scratch must not alias metadata"))
    end
    return nothing
end

function _scaled_input_magnitude(f::AbstractArray{T}) where {T}
    magnitude = zero(T)
    for value in f
        isfinite(value) || throw(ArgumentError("spectral operator input values must be finite"))
        magnitude = max(magnitude, abs(value))
    end
    return magnitude
end

function _scaled_input_loses_nonzero(f, scale)
    isone(scale) && return false
    for value in f
        _fft_input_loses_nonzero(value, value / scale) && return true
    end
    return false
end

function _require_scaled_values(values)
    for value in values
        isfinite(value) || throw(ArgumentError("normalized Fourier data must be finite"))
    end
    return nothing
end

# Unit-gain headroom applies to the normalized multiplier. A tiny entire input
# is copied upward exactly before the forward FFT; its inverse shift is retained.
function _prepare_scaled_input!(g::FourierGrid{D,T}, f, magnitude) where {D,T}
    safe = _fft_safe_magnitude(T, length(f), one(T))
    input_scale = _fft_input_scale(magnitude, safe)
    _copy_fft_input!(g.cbuf, f, input_scale)
    prepared_magnitude = zero(T)
    for value in g.cbuf
        prepared_magnitude = max(prepared_magnitude, abs(real(value)))
    end
    input_exponent = 0
    if !iszero(prepared_magnitude) && prepared_magnitude < one(T)
        input_exponent = exponent(prepared_magnitude)
        for I in eachindex(g.cbuf)
            g.cbuf[I] = ldexp(real(g.cbuf[I]), -input_exponent)
        end
    end
    return input_scale, input_exponent
end

@inline function _scaled_axis(g::FourierGrid{D}, axis::Integer) where {D}
    1 <= axis <= D || throw(ArgumentError("axis $axis out of range 1:$D"))
    return Int(axis)
end

@inline _ceil_log2(n::Int) = n <= 1 ? 0 : ndigits(n - 1; base = 2)

# Keep the concrete IEEE type through ldexp during abstract package inference.
# These are the unchanged ordinary-path floor expressions for the FFTW types.
@inline _laplacian_input_floor(::Type{Float32}, n::Int) =
    ldexp(floatmin(Float32), _ceil_log2(n) + 4)
@inline _laplacian_input_floor(::Type{Float64}, n::Int) =
    ldexp(floatmin(Float64), _ceil_log2(n) + 4)

function _ordinary_laplacian_products(g::FourierGrid{D,T}) where {D,T}
    emin = exponent(floatmin(T))
    emax = exponent(floatmax(T))
    volume_power = _ceil_log2(length(g.cbuf))
    active_axes = count(>(1), g.n)
    upper = emax - _ceil_log2(max(1, active_axes)) - 4
    lower = emin + volume_power + 4
    for I in CartesianIndices(g.cbuf)
        h = g.cbuf[I]
        for value in (real(h), imag(h))
            iszero(value) && continue
            ex = exponent(abs(value))
            for d = 1:D
                k = g.kfull[d][I[d]]
                iszero(k) && continue
                ek = exponent(abs(k))
                ex + ek >= emin && lower <= ex + 2ek <= upper || return false
            end
        end
    end
    return true
end

@inline function _scaled_odd_wavenumber(g, I, axis, ::Val{K}) where {K}
    if K === :nyquist
        return iseven(g.n[axis]) && I[axis] == g.n[axis] ÷ 2 + 1 ? g.kfull[axis][I[axis]] :
               zero(eltype(g.L))
    end
    return g.kvec[axis][I[axis]]
end

# One exponent is selected over every active real/imaginary component and axis.
# Local metadata are normalized separately, so an inactive stiff axis cannot
# underflow an active small-axis squared ratio.
function _normalized_mode_power(g::FourierGrid{D}, axis, kind::Val{K}) where {D,K}
    power = typemin(Int)
    for I in CartesianIndices(g.cbuf)
        odd = K === :laplacian ? one(eltype(g.L)) : _scaled_odd_wavenumber(g, I, axis, kind)
        iszero(odd) && continue
        odd_power = K === :derivative_laplacian ? exponent(abs(odd)) : 0
        h = g.cbuf[I]
        for value in (real(h), imag(h))
            iszero(value) && continue
            ex = exponent(abs(value))
            if K === :laplacian || K === :derivative_laplacian
                for d = 1:D
                    k = g.kfull[d][I[d]]
                    iszero(k) || (power = max(power, ex + 2exponent(abs(k)) + odd_power))
                end
            else
                power = max(power, ex + exponent(abs(odd)))
            end
        end
    end
    return power
end

@inline function _normalized_component(
    value,
    g::FourierGrid{D,T},
    I,
    axis,
    power,
    kind::Val{K},
) where {D,T,K}
    iszero(value) && return value
    ex = exponent(abs(value))
    unit = ldexp(value, -ex)
    if K === :laplacian || K === :derivative_laplacian
        odd = K === :derivative_laplacian ? g.kvec[axis][I[axis]] : one(T)
        iszero(odd) && return zero(T)
        eo = K === :derivative_laplacian ? exponent(abs(odd)) : 0
        odd_unit = K === :derivative_laplacian ? ldexp(odd, -eo) : one(T)
        result = zero(T)
        for d = 1:D
            k = g.kfull[d][I[d]]
            iszero(k) && continue
            ek = exponent(abs(k))
            ku = ldexp(k, -ek)
            term = -(unit * ku) * ku
            K === :derivative_laplacian && (term *= odd_unit)
            result += ldexp(term, ex + 2ek + eo - power)
        end
        return result
    end
    k = _scaled_odd_wavenumber(g, I, axis, kind)
    iszero(k) && return zero(T)
    ek = exponent(abs(k))
    return ldexp(unit * ldexp(k, -ek), ex + ek - power)
end

# Keep scalar and linear index types correlated at the two cold zero owners.
# Base.vec shares storage; return the original buffer with its original shape.
function _zero_scaled_buffer!(buffer::Array{Complex{T},D}) where {T,D}
    fill!(vec(buffer), zero(Complex{T}))
    return buffer
end

# Keep the stored array and its CartesianIndex dimensions correlated.
# The value conversion and bounds checks remain owned by Base.setindex!.
@inline function _store_normalized_cartesian!(
    buffer::Array{Complex{T},D},
    value,
    I::CartesianIndex{D},
) where {T,D}
    buffer[I] = value
    return nothing
end

function _apply_normalized_multiplier!(g, axis, kind::Val{K}) where {K}
    _require_scaled_values(g.cbuf)
    power = _normalized_mode_power(g, axis, kind)
    if power == typemin(Int)
        _zero_scaled_buffer!(g.cbuf)
        return 0
    end
    for I in CartesianIndices(g.cbuf)
        h = g.cbuf[I]
        r = _normalized_component(real(h), g, I, axis, power, kind)
        i = _normalized_component(imag(h), g, I, axis, power, kind)
        _store_normalized_cartesian!(
            g.cbuf,
            K === :derivative || K === :derivative_laplacian ? complex(-i, r) : complex(r, i),
            I,
        )
    end
    return power
end

function _normalized_scaled_data!(g, f, magnitude, axis, kind)
    input_scale, input_exponent = _prepare_scaled_input!(g, f, magnitude)
    g.plan * g.cbuf
    power = _apply_normalized_multiplier!(g, axis, kind)
    g.iplan * g.cbuf
    _require_scaled_values(g.cbuf)
    return (; values = g.cbuf, input_scale, exponent = power + input_exponent, ordinary = false)
end

function _zero_scaled_data!(g::FourierGrid{D,T}) where {D,T}
    _zero_scaled_buffer!(g.cbuf)
    return (; values = g.cbuf, input_scale = one(T), exponent = 0, ordinary = false)
end

"""
    laplacian_scaled_data!(g, f)

Return borrowed normalized full-Laplacian data. Its physical real value is
`real(values) * input_scale * 2^exponent`. `ordinary=true` preserves the existing
serial Laplacian arithmetic. Other cases normalize individual mode/axis products
before they can underflow. Input rescaling rejects lost nonzero values.

The result aliases grid scratch and is valid only until the next grid operator.
Input must be disjoint from mutated scratch. The grid remains one-caller-only.
See `docs/FOURIER_BACKEND_CONTRACT.md` for the complete qualified API contract.
"""
function laplacian_scaled_data!(g::FourierGrid{D,T}, f::AbstractArray) where {D,T}
    _require_fftw_float(T, :laplacian_scaled_data!)
    _require_scaled_input(f, g)
    gain = _laplacian_gain(g.kfull)
    magnitude = _scaled_input_magnitude(f)
    scale = _fft_input_scale(magnitude, _fft_safe_magnitude(T, length(f), gain))
    floor = _laplacian_input_floor(T, length(f))
    if (iszero(magnitude) || magnitude >= floor) && !_scaled_input_loses_nonzero(f, scale)
        _copy_fft_input!(g.cbuf, f, scale)
        g.plan * g.cbuf
        _require_scaled_values(g.cbuf)
        if _ordinary_laplacian_products(g)
            _apply_laplacian_multiplier!(g.cbuf, g.kfull)
            g.iplan * g.cbuf
            _require_scaled_values(g.cbuf)
            return (; values = g.cbuf, input_scale = scale, exponent = 0, ordinary = true)
        end
    end
    return _normalized_scaled_data!(g, f, magnitude, 1, Val(:laplacian))
end

"""
    deriv_scaled_data!(g, f, axis)

Return borrowed normalized data for the canonical centered first derivative,
zeroing even-grid Nyquist. Representation/lifetime match `laplacian_scaled_data!`;
this cold operation always returns `ordinary=false`.
"""
function deriv_scaled_data!(g::FourierGrid, f::AbstractArray, axis::Integer)
    _require_scaled_input(f, g)
    j = _scaled_axis(g, axis)
    magnitude = _scaled_input_magnitude(f)
    g.n[j] <= 2 && return _zero_scaled_data!(g)
    return _normalized_scaled_data!(g, f, magnitude, j, Val(:derivative))
end

"""
    nyquist_scaled_data!(g, f, axis)

Return borrowed normalized data for the real positive Nyquist projector channel.
Odd axes produce zero, after input validation. Representation/lifetime match
`laplacian_scaled_data!`; this cold operation returns `ordinary=false`.
"""
function nyquist_scaled_data!(g::FourierGrid, f::AbstractArray, axis::Integer)
    _require_scaled_input(f, g)
    j = _scaled_axis(g, axis)
    magnitude = _scaled_input_magnitude(f)
    isodd(g.n[j]) && return _zero_scaled_data!(g)
    return _normalized_scaled_data!(g, f, magnitude, j, Val(:nyquist))
end

"""
    deriv_laplacian_scaled_data!(g, f, axis)

Return borrowed normalized data for the centered derivative of the full-k
Laplacian. Products are composed before either physical intermediate is rounded.
The existing full squared-wavenumber admission applies. Representation/lifetime
match `laplacian_scaled_data!`; this cold operation returns `ordinary=false`.
"""
function deriv_laplacian_scaled_data!(g::FourierGrid, f::AbstractArray, axis::Integer)
    _require_scaled_input(f, g)
    j = _scaled_axis(g, axis)
    _laplacian_gain(g.kfull)
    magnitude = _scaled_input_magnitude(f)
    g.n[j] <= 2 && return _zero_scaled_data!(g)
    return _normalized_scaled_data!(g, f, magnitude, j, Val(:derivative_laplacian))
end
