# Canonical scaled-Laplacian floor repair: pre-edit work order V1

Sealed design and verification record (2026-10-29): the concrete floor, canonical public type admission,
zero-buffer and Cartesian store repairs are integrated and pass the complete
both-runtime normal campaign and independent root postflight. Actual results and
release provenance and the separate consumer boundary are owned by the
[verification record](SCALED_FOURIER_VERIFICATION.md). Earlier failed candidates,
including the first docstring AST fixture and public zero-fill reports, remain
preserved. Their failures do not describe the final candidate.

The reviewed pre-edit work order follows as the design and fixed gate contract;
its future-tense source proposal records the authorized change, not an extra scope.

## Verified finding and limits

`probe-v3-actual-classification-v1.json` binds both actual diagnostic calls, all
34 source/output inputs and the complete captures. Julia1.12.7 and1.13.3 each
passed13 target checks plus3051 inference/exact-value assertions. Each runtime
retains eleven captures: exact abstract ldexp and unbounded generic floor have
two reports each; supported-union generic/caller and seven concrete entries have
zero. Eight complete Base frame/signature comparisons reproduce the old public
package reports. The760 exact floor records additionally passed independent
Python integer-bit-length, Fraction and IEEE-bit references, shifts4..67.

This establishes the abstract inference path and a report-free supported scalar
prototype. That diagnostic experiment alone did not establish final public/package acceptance.
The later complete local disposition is linked above; native peak memory and a
compiler-clean-package claim remain outside it.

## Canonical owners and proposed source change

At595, `src/scaled_differential_operators.jl:208–213` owns public Laplacian
admission and the floor. `:72` owns the unchanged ceil-log2. The floor controls
whether the unchanged ordinary FFT path is attempted; no operator formula,
multiplier, FFT, scaling, borrowed storage or physical composition changes.

`FourierGrid` has unconstrained T (`periodic_fourier.jl:25`). Its canonical
`_require_fftw_float` (`:41–45`) admits only Float32/Float64 AbstractFloat types,
and the public constructor already uses it. `_require_scaled_input` currently
validates matching arrays/aliases/metadata, but does not call that type owner.
The borrowed API expressly supports matching real F32/F64 only
(`docs/FOURIER_BACKEND_CONTRACT.md:95–98`).

In a NEW worktree from595, propose exactly:

1. Two private `_laplacian_input_floor` methods, one for `Type{Float32}` and one
   for `Type{Float64}`, each containing the existing
   `ldexp(floatmin(Float32/Float64), _ceil_log2(n)+4)` expression, with n::Int.
   No generic numeric fallback or new scalar implementation.
2. At public `laplacian_scaled_data!` entry, call the existing
   `_require_fftw_float(T, :laplacian_scaled_data!)` before array/scratch work.
   Preserve the public signature with unconstrained T. Replace only the floor
   expression by `_laplacian_input_floor(T, length(f))` at its existing position.
   Existing supported-array validation and all subsequent operation order remain.
3. Update the candidate's live Markdown tracker and document explicit admission;
   do not modify old frozen design/evidence files or global trackers root owns.

Unsupported manually built Float16/BigFloat grids will receive the canonical
ArgumentError before scratch writes or plan calls. Manually assembled non-
AbstractFloat grids remain outside the backend contract: the canonical validator
has no such method, so MethodError is the declared rejection for that category.
Do not add generic arithmetic fallback, claim support, mutate buffers, or require
unsupported fake plans to execute. Supported F32/F64 behavior must preserve bits.
A new universal unsupported-type error policy is outside this minimal repair.

## Tests and fixed acceptance gates before measurements

- Actual full unbounded public signature: capture default JET.report_call on
  `Tuple{typeof(laplacian_scaled_data!),FourierGrid{D,T},AbstractArray} where {D,T}`
  after unique-target preflight. Retain all results. The supported scalar union
  alone cannot satisfy this gate. Capture both concrete private methods too.
- Full report_package using unchanged normal harness configuration: compare
  complete595 parent package to repair on both runtimes. Expected only the two
  floor reports disappear; the existing constructor report must match under an
  explicit root and checked source-line map. Preserve every unlimited render,
  frame field, specialization/static parameter, definition and full source.
  No filtering or changing analyzer settings. All60 previous concrete captures
  and the new public/helper captures must be empty.
- Permanent exact Rational/IEEE-bit helper tests use all original n/type controls,
  powers-of-two neighbors, extremes and n<=0 scalar compatibility. Run normal
  supported calls at adjacent floor thresholds, including subnormal, odd/even,
  constant/zero, anisotropic and ordinary/cold routing cases. Compare floor bits
  and complete borrowed result fields/values to595. Existing DFT/local gamma+
  omission, range/alias/lifetime/error and all package tests stay unchanged.
- Unsupported manual-grid tests use typed finite metadata/scratch and a counting
  sentinel plan. Float16 and BigFloat must throw canonical ArgumentError; one
  non-AbstractFloat fixture must throw MethodError. All scratch, input and
  metadata bytes/values and plan-call counts remain unchanged. These fixtures
  test rejection only and do not assert unsupported FFT support. Supported bad
  shape/type/alias/error precedence remains tested against595.
- Mutants must fail if either concrete type is swapped, floor shift changes,
  admission is removed, unsupported plan runs, ordinary floor is bypassed or the
  unbounded public/report-package gate is replaced with supported-only capture.
  Mutants run in separate disposable evidence, never alter frozen parent.
- Normal warmed Julia allocation remains exactly0B for every existing API/resource
  profile, with grid-owned storage equal595. Private helper calls at all scalar
  n/type boundaries must also allocate0B. No gate loosening, resource cap growth,
  timing acceptance under contention or native-peak claim. Existing native peak
  accounting stays explicitly open.
- Run full package/Aqua, concrete/full quality, formatter1.0.62, AST/caller audit
  and numerical/resource gates on min1.12.7/latest1.13.3 in the serial slot only
  after a frozen exact command/input packet is reviewed. Lock dependency/provider
  graphs: only the declared SO worktree path/compiled-image identity may change;
  all raw graphs/DLLs and complete image/source/CRC proofs remain retained.
  Pkg.develop/resolve must not silently update dependencies. No retry overwrites.

## Execution and handoff sequence

- [x] Both V3 diagnostic calls and independent complete-frame/exact-bit checks.
- [x] This source-supported pre-edit brief and explicit unsupported-T semantics.
- [x] Root review/release; create new worktree from595 and write minimal patch.
- [x] Permanent tests/docs written; formatter12 and corrected whole-source AST20
checks pass. Both exact provider-graph/hash preparations pass.
- [x] Final floor/zero/store normal campaign and complete root postflight.
- [x] Classify every full compiler and resource delta under unchanged gates.
- [x] Dependency publication/main push, four exact-SHA CI jobs and registration/release.
- [ ] Separate HPIC36 consumer integration and its complete acceptance gates.

Any added report, changed supported bit pattern, weakened zero-allocation result
or undeclared provider/source drift holds promotion; retain the failed attempt.
Local SO acceptance does not imply HPIC36 acceptance.
