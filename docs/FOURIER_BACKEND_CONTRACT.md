# Qualified Fourier backend API

These names are supported through `SpectralOperators.<name>` and are marked
`public`; they are intentionally not exported. They are the same numerical
owners used by the serial operators. They do not plan, execute or normalize an
FFT, communicate, allocate field workspaces, or validate distributed ownership.

`_fourier_metadata(n, L)` accepts the same nonempty tuple of positive `Int` sizes
and tuple of positive finite Float32 or Float64 lengths as `FourierGrid`. It
performs the same spacing/wavenumber representability checks and returns a
NamedTuple with `n, L, dx, midx, kfull, ik, kvec`. Its vectors are full one-based
per-axis metadata. `midx` uses positive even-grid Nyquist; `kfull` retains it;
`ik` and `kvec` zero it for odd derivatives and projection. It allocates only the
metadata, not plans or field storage. Treat these vectors as read-only.

## Mode kernels

For `T = Float32` or `Float64`, fields have type `AbstractArray{Complex{T},D}`.
Their axes are **global logical Fourier indices**, even when storage is permuted.
Do not pass local array indices for a tile that starts at a nonzero mode offset.
All arrays in a call have identical logical axes. Nonempty axes must be subsets
of the relevant full metadata vectors. Empty tiles are valid no-ops. Metadata
vectors must come from a compatible `_fourier_metadata` result; kernels do not
rescan their numerical values. Their shapes/ranges and aliases are checked once
before the multiplier writes any output.

| Function | Operation and additional arguments |
| --- | --- |
| `_apply_ik_store!(dst, hat, ikj, j)` | `dst[I] = hat[I] * ikj[I[j]]`; `1 ≤ j ≤ D` |
| `_apply_ik_store_scaled!(dst, hat, ikj, j, s)` | `s * hat[I] * ikj[I[j]]` in that order; `s::T` |
| `_apply_ik_accum!(dst, hat, ikj, j, s)` | Adds the preceding ordered product |
| `_curl_component_store!(dst, plus, plus_ik, plus_axis, minus, minus_ik, minus_axis)` | Ordered difference of the two derivative products |
| `_apply_laplacian_multiplier!(hat, kfull)` | In-place sum of `-(hat[I]*k[d])*k[d]`, retaining full Nyquist and representable products whose isolated `k²` underflows |
| `_project_fourier_modes!(hats, kvec)` | In-place projection of a tuple of three complex arrays for `1 ≤ D ≤ 3`; absent-axis components and the derivative-zero nullspace remain unchanged |
| `_apply_exp_filter_multiplier!(hat, midx, n, alphaT, p)` | In-place exponential mode-ratio filter; `alphaT::T ≥ 0` (Inf allowed), `p::Integer ≥ 1`; DC retained exactly |
| `_apply_two_thirds_mask!(hat, midx, n)` | Retains precisely `3abs(midx[d][I[d]]) < n[d]` on every axis |

Each mutating kernel returns its destination (the original tuple for projection).
Exact object identity between output and input is safe for the derivative and
curl component kernels and is allowed. Other potentially overlapping storage is
rejected. Projection's three components must be pairwise disjoint. Outputs must
not alias metadata. Array adapters must implement `Base.dataids`/`Base.mightalias`
correctly; hidden overlapping storage violates the contract. These rules concern
complex mode arrays; the existing real-space serial API retains its own alias
rules and checks.

These primitives do not establish finite transformed values. Their caller owns
input validation, range scaling and FFT normalization. For Laplacian, call
`_laplacian_gain(kfull)` before transforming: it returns the existing global
squared-wavenumber bound or throws `ArgumentError` if that bound would overflow.
Projection normalizes both the derivative wavenumber and complex field components
before forming its dot product, preserving the serial extreme-range behavior.

## Range and conversion primitives

The following scalar/array helpers retain their existing arithmetic and failure
behavior. These are caller preconditions, not additional per-element scans:

1. `_fft_safe_magnitude(T, global_count::Int, gain::T)` requires a positive global
   transform volume and a finite nonnegative gain. It uses the global volume,
   never a rank's owned-cell count. Invalid gain throws `ArgumentError`.
2. `_fft_input_scale(magnitude::T, safe_magnitude::T)` consumes a finite
   nonnegative global maximum and finite nonnegative safe magnitude. Zero safe
   magnitude can result from underflow. The helper returns one when no scaling
   is needed, otherwise the existing upward-biased ratio, capped at `floatmax(T)`.
   Only this scalar overload is the transport-independent backend contract.
3. `_copy_fft_input!(dst, field, scale)` copies real `T` values to complex `T`
   storage with identical logical axes and disjoint storage, dividing by finite
   `scale::T ≥ 1` when needed. It rejects loss of any nonzero value to zero.
   A failure may leave partial scratch writes; do not publish that scratch.
4. `_store_fft_output!(out, transformed, scale)` extracts real values after a
   **normalized inverse** and restores the same scale once. Matching logical
   axes, real/complex `T` types and disjoint storage are caller preconditions.
   Unrepresentable results follow the existing floating-point behavior.

A distributed transport must validate ranks/components, reduce finite input
magnitudes globally, agree a common scale, and agree any lost-nonzero error before
entering another transpose. Those collectives and transactional output publication
belong to the transport implementation. This package has no MPI dependency.

The serial `FourierGrid` layout and operation order remain unchanged. Its scratch
buffers still belong to one caller at a time; this API does not make shared
workspaces thread-safe. The extraction does not promise GPU kernels or distributed
solver/checkpoint support.
