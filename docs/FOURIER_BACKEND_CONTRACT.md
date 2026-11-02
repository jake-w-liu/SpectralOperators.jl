# Qualified Fourier backend API

These names are supported through `SpectralOperators.<name>` and are marked
`public`; they are intentionally not exported. The mode and range primitives
below are the numerical owners used by serial operators; they do not plan,
execute or normalize an FFT, communicate, allocate field workspaces, or validate
distributed ownership. The borrowed scaled-data operations described at the
end execute existing grid FFT plans.

`_fourier_metadata(n, L)` accepts the same nonempty tuple of positive `Int` sizes
and tuple of positive finite Float32 or Float64 lengths as `FourierGrid`. It
performs the same spacing/wavenumber representability checks and returns a
NamedTuple with `n, L, dx, midx, kfull, ik, kvec`. Its vectors are full one-based
per-axis metadata. `midx` uses positive even-grid Nyquist; `kfull` retains it;
`ik` and `kvec` zero it for odd derivatives and projection. It allocates only the
metadata, not plans or field storage. Treat these vectors as read-only.

## Serial real-field storage

The serial `FourierGrid` operators accept real fields with the grid's type,
shape and one-based axes. Inputs and outputs must not overlap the internal
complex buffers that the complete operation mutates, including through ordinary
real reinterpretations or views. Writable outputs must not overlap any axis of
read-only `midx`, `kfull`, `ik` or `kvec` metadata. A shared entry preflight rejects
these conflicts with `ArgumentError` before scratch, FFT or output writes. Array
adapters must report aliases truthfully through `Base.dataids`/`Base.mightalias`.

| Operator | Mutated scratch |
| --- | --- |
| `deriv!`, `gradient!` | `cbuf`, `tbuf` |
| `laplacian!`, `exp_filter!`, `dealias_two_thirds!` | `cbuf` |
| `nyquist_channel!` | `cbuf` on an even axis; none on an odd axis |
| `divergence!` | `cbuf`, `abuf`; only the first `D` input components participate |
| `curl!` | D1: `cbuf`, `abuf`; D2: also `tbuf`; D3: also `sbuf` |
| `project_divfree!` | First `D` of `cbuf`, `tbuf`, `abuf` |

Projection protects all three field components from scratch writes, including
absent-axis components that must remain unchanged. Only its first `D` components
are writable outputs for the metadata rule. Odd-axis Nyquist uses no scratch but
still rejects a metadata output. Inputs may alias read-only metadata; real fields
may alias scratch not mutated by the operation. These checks assume the grid
constructor's internally disjoint storage; they do not validate manually assembled
overlapping grids. Each concurrent caller needs its own grid.

Existing user-field sharing rules remain: scalar input/output overlap is allowed;
`gradient!` outputs are pairwise disjoint but one may overlap its input;
`divergence!` permits input/output and inter-input sharing; filters are in place.
`curl!` outputs remain pairwise disjoint and disjoint from all inputs;
`project_divfree!` components remain pairwise disjoint. The qualified complex
kernels and borrowed operations below retain their separate contracts. This
storage preflight is a source candidate; remote numerical, allocation and package
gates remain pending.

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
| `_apply_exp_filter_multiplier!(hat, midx, n, alphaT, p)` | In-place exponential mode-ratio filter; `alphaT::T ≥ 0` (Inf allowed), `p::Integer ≥ 1`; DC retained exactly; positive infinite strength removes every non-DC mode, including underflowed mode powers |
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

The following scalar/array helpers define the shared transform-range policy.
These are caller preconditions, not additional per-element scans:

1. `_fft_safe_magnitude(T, global_count::Int, gain::T)` requires a positive global
   transform volume and a finite nonnegative gain. It uses the global volume,
   never a rank's owned-cell count. Headroom includes the unnormalized inverse
   before its final normalization: the spectral 1-norm is bounded by
   `gain * N * sqrt(N)` times the maximum input, with the existing factor-eight
   margin and a minimum unit gain. For the three-component orthogonal projector,
   the combined input bound adds at most sqrt(3), covered by that margin.
   Invalid gain throws `ArgumentError`.
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

The serial `FourierGrid` layout and transform/multiplier order remain unchanged.
D1 serial `curl!` chooses its common input scale from active components `A[2]`
and `A[3]`; unused `A[1]` retains finite/type/shape/alias validation. D2/3 retain
all-component scaling. The scalar/backend primitives above are unchanged, and
active copy loss still rejects. A changed D1 scale may change rounded output;
this does not establish exact FFT results or cross-version trajectory identity.
The selector's source implementation still needs exact-revision remote gates.
Scratch buffers still belong to one caller at a time; this API does not make shared
workspaces thread-safe. The extraction does not promise GPU kernels or distributed
solver/checkpoint support.

## Borrowed scaled Fourier data

`laplacian_scaled_data!(g, f)`, `deriv_scaled_data!(g, f, axis)`,
`nyquist_scaled_data!(g, f, axis)` and
`deriv_laplacian_scaled_data!(g, f, axis)` accept a serial `FourierGrid` and
matching real Float32/Float64 input. Axis is an integer in `1:ndims(f)`.
`laplacian_scaled_data!` also checks this floating-point domain at entry, before
scratch writes or FFT execution, for manually assembled grids that bypass the
public constructor. The canonical FFTW check rejects Float16/BigFloat with
`ArgumentError`; non-`AbstractFloat` manual grids are outside its method domain.
They return the same concrete NamedTuple:

| Field | Meaning |
| --- | --- |
| `values` | Borrowed complex grid scratch, after a normalized inverse FFT |
| `input_scale::T` | Finite common input scale, at least one |
| `exponent::Int` | Additional binary power retained outside floating-point storage |
| `ordinary::Bool` | True only when the Laplacian used unchanged serial arithmetic |

The physical real result is `real(values)*input_scale*2^exponent`. Keep these
factors through the downstream constitutive expression when an intermediate
physical result could underflow or overflow. Converting the factors separately
can lose the range these operations preserve. The caller owns physical output
conversion/admission; the operations require finite normalized returned data,
not representability of a separately restored physical intermediate.

Consume each result before another operator on the same grid. `values` aliases
`g.cbuf` and has no retained lifetime, reentrancy or thread-sharing guarantee.
Input must be disjoint from that mutated scratch. Scratch must be disjoint from
readonly Fourier metadata. Matching type, shape, one-based axes and these aliases
are checked before writes. A readonly input which aliases metadata is safe if
it is disjoint from mutated scratch. Failed input rescaling can leave partial
scratch writes; no physical output has been published there.

The Laplacian retains full even-grid Nyquist. The first derivative uses centered
metadata, zeroing even Nyquist. The Nyquist operation is the existing positive
real projector channel, zero on odd axes after finite-input validation. The
derivative-Laplacian composes full-k squares with the centered odd derivative
before either physical intermediate is rounded. Laplacian and
derivative-Laplacian require the existing `_laplacian_gain` full-grid squared
wavenumber admission even when an output channel is zero. No metadata/domain
rescaling or change of Nyquist convention occurs.

Individual mode/axis products use exact power-of-two normalization. This avoids
losing a tiny active axis through a squared ratio against an inactive stiff one.
Normalized transforms use the existing unit-gain input headroom and reject any
nonzero input lost during necessary common downscaling. A tiny entire input is
copied upward by an exact binary power before its forward FFT. Zero derivative
axes and odd Nyquist axes need no transforms after validation. Ordinary
Laplacian inputs retain the exact serial transform/multiplier operation order;
conservative exponent headroom checks select normalization otherwise.

Before a directional cold FFT, a certified preparation can subtract each
axis-line's first value. The removed field is constant in the differentiated
axis, so centered derivative, Nyquist projector and derivative-Laplacian
annihilate it. The full Laplacian does not use this preparation. Exact binary
operand scaling and a zero error-free subtraction residual are required in
nearest/gradual Float32/64 arithmetic; otherwise the original preparation runs.
Huge axis-constant lines become exact zero before scaling. The physical result
still uses the same returned scale/exponent and borrowed scratch lifetime.
See [axis conditioning and verification](SCALED_AXIS_CONDITIONING_PLAN.md).
This removes an unrelated axis-null background, not arbitrary active-mode FFT
conditioning or cancellation error.

Normalized scalar multiplier error is bounded at its positive operation scale
by `gamma_(D+3)*sum(abs(exact terms))`, plus at most half a least-subnormal unit
per shifted axis term. This is distinct from FFTW error and cancellation in the
complete real-space result. High-precision DFT tests check the FFT/source-scale
behavior; these operations do not promise exact Fourier coefficients or a
relative error bound on every arbitrarily canceled output.

## Borrowed scaled curl terms

`foreach_curl_scaled_term!(consume, A, g, component)` is a qualified, nonexported
serial API for Float32/64 Fourier grids in dimensions 1–3. It requires nearest
rounding and gradual underflow throughout the call. `A` is a three-tuple of real
arrays matching the grid type, shape and one-based axes. The integer component
index must be in `1:3`; it is checked before conversion to `Int`.

The canonical curl geometry supplies these terms, in positive-then-negative order:

| Component | Positive term | Negative term |
| --- | --- | --- |
| 1 | ∂₂ A₃ | ∂₃ A₂ |
| 2 | ∂₃ A₁ | ∂₁ A₃ |
| 3 | ∂₁ A₂ | ∂₂ A₁ |

Absent axes are omitted. For each present axis the existing `deriv_scaled_data!`
owner produces a receipt and immediately calls `consume(term, data)`. `term`
contains integer fields `component,input,axis,sign`, with sign +1 or -1. `data`
has the borrowed fields described above and `ordinary=false`. Present axes of
length one/two emit their exact zero receipt without a transform; even Nyquist
uses the centered first-derivative convention. Return value is the number of
callbacks, 0–2, independent of callback return values. D1 counts are 0/1/1, D2
counts are 1/1/2 and D3 counts are 2/2/2. A count is geometric, never a nonzero
certificate. No physical floating-point curl is formed or summed by this API.

All three inputs, including inactive A₁ in D1, pass finite/type/shape/axis and
scaled-input storage checks before any scratch write, transform or callback.
Only `g.cbuf` is mutated; it must be disjoint from all inputs and Fourier
metadata. Readonly input sharing and aliases to genuinely unused scratch remain
valid. Storage checks use `Base.mightalias`; callers must use arrays whose
overlap those interfaces report and constructor-compatible internally disjoint
grid storage. Initial admission errors leave scratch unchanged and invoke no
callback. Later derivative rescaling/FFT failures can follow earlier callbacks
and can leave scratch modified. Such errors and callback exceptions propagate
without retries or replacement results.

Consume each receipt before returning: the next term overwrites its values.
Callbacks must not mutate inputs, metadata or arithmetic context, reenter an
operator on this grid, retain borrowed values, or share it concurrently. SO does
not roll back arbitrary callback side effects. Transactional callers must keep
effects private and defer physical publication until all required components
and final finite conversions succeed. A component argument permits reuse of one
private accumulator across three calls; this API supplies no accumulator or
solver publication policy.

Keep each term's input scale and binary exponent through the complete downstream
expression. Restoring terms to the field type before multiplying a coefficient
can lose both overflow and underflow recovery. The result is a sum of represented
normalized derivative channels, with the same input-loss, FFT conditioning and
source-scale error qualifications as the borrowed directional owner. It does
not promise exact DFT coefficients or relative accuracy under arbitrary
cancellation. Separate directional transforms may round differently from
ordinary `curl!`; a nonfinite-only ordinary fallback cannot detect an already
underflowed channel. Downstream map selection therefore requires its own review.

Ordinary `curl!` shares only the component/axis descriptor and retains its
existing transform, multiplier and store order. For generic active axes, three
new component calls require 4/8/12 total forward/inverse transforms in D1/D2/D3,
compared with ordinary curl's 4/6/6, plus repeated all-input preflight. This is
additional source-recovery work, not a speed claim. The standalone implementation
and essential controls in `test/test_scaled_sources.jl` and
`test/test_operators.jl` require exact-revision remote checks; no downstream
consumer or new release is implied.
