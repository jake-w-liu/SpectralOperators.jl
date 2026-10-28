# Shared Fourier kernels: CAP-003A

Parent: `cc515cfc3396e616b8e42b9539e1effcab2a069c` (0.1.2).
Scope: shared numerical primitives only. MPI transport, reductions, ownership,
Pencil storage, distributed solvers and HPIC precision remain separate work.

## Pre-edit brief and compatibility

- Mode: bounded extraction and supported backend API.
- Requirement: reuse the serial Fourier arithmetic on globally indexed complex
  arrays without constructing FFT plans or global field storage.
- Source evidence: `periodic_fourier.jl` owns metadata, FFT range and rescaling;
  `differential_operators.jl`, `divergence_control.jl` and `filters.jl` own the
  multiplier loops. `fdx_fourier_yz.jl` also consumes range helpers. Existing
  tests exercise odd/even modes, extreme ranges, aliases and zero allocations.
- Gap: metadata construction and Laplacian/projection/filter multiplication are
  currently coupled to serial FFT orchestration. This is missing backend reuse,
  not a demonstrated scientific defect in the serial implementation.
- Invariants: unchanged FourierGrid fields, exported serial signatures, FFT
  ordering, floating-point expression ordering, Nyquist rules, dynamic-range
  rejection, real-array alias behavior and output publication on existing errors.
- Strategy: extract the existing loops once. Keep transport orchestration in its
  current owner; no grid hierarchy, duplicate physics, dependency or wrapper API.
- Verification: immutable parent source, both Julia 1.12.7 and 1.13.2, existing
  suite/Aqua, independent dense DFT and analytic modes, bitwise serial output and
  workspace parity, global-offset/permuted tiles, failure/alias controls,
  ordinary allocation parity and concrete/full compiler reports.

## Qualified backend contract

The following names are supported through `SpectralOperators.name`; exports stay
unchanged. Underscores preserve the existing canonical names, not private status.

| Name | Inputs and result |
| --- | --- |
| `_fourier_metadata(n,L)` | Same nonempty Int/Float32-or-Float64 tuples and checks as FourierGrid; returns `(n,L,dx,midx,kfull,ik,kvec)` NamedTuple without plans or field arrays |
| `_apply_ik_store!(dst,hat,ik,j)` | Pointwise `hat*ik[j]` |
| `_apply_ik_store_scaled!(dst,hat,ik,j,s)` | Pointwise `s*hat*ik[j]`, preserving order |
| `_apply_ik_accum!(dst,hat,ik,j,s)` | Adds that product |
| `_curl_component_store!(dst,plus,plus_ik,plus_axis,minus,minus_ik,minus_axis)` | Existing two-product difference |
| `_apply_laplacian_multiplier!(hat,kfull)` | Existing ordered `(hat*k)*k` sum, in place |
| `_laplacian_gain(kfull)` | Existing global finite squared-wavenumber bound/rejection |
| `_project_fourier_modes!(hats,kvec)` | Three complex component arrays, spatial dimension 1–3; existing normalized projection, absent components untouched |
| `_apply_exp_filter_multiplier!(hat,midx,n,alphaT,p)` | Existing mode-ratio filter, including exact DC retention |
| `_apply_two_thirds_mask!(hat,midx,n)` | Existing strict `3abs(m)<n` mask |
| `_fft_safe_magnitude(T,global_count,gain)` | Existing transform/multiplier headroom bound |
| `_fft_input_scale(magnitude,safe_magnitude)` | Existing scalar scale selection |
| `_copy_fft_input!(dst,field,scale)` | Existing real-to-complex copy and lost-nonzero rejection |
| `_store_fft_output!(out,transformed,scale)` | Existing real extraction and rescaling |

Mode arrays use `AbstractArray{Complex{T},D}` with Cartesian global **logical**
indices, not local memory-order indices. All participating arrays have identical
axes. Axes may be non-one-based subsets of metadata indices, including permuted
storage; nonempty axes must lie within the relevant full one-based metadata
vectors. Empty raw tiles are legal no-ops (distributed construction may reject
them). Metadata is the immutable-by-contract output of `_fourier_metadata`, kept
full and replicated; callers must not alter it while operators run. Once-per-call
axis/range/alias checks precede mutation; no per-cell metadata validation.

Exact object identity between destination and a multiplier input is safe and
allowed; other potentially overlapping storage rejects. Projection components
must be pairwise disjoint, including absent components. Output must not alias
metadata. Matching axes and element types are required for copy/store; their real
and complex storage must be disjoint. Backend array implementations must implement
`Base.mightalias`/`Base.dataids` truthfully. The existing real serial public alias
contract is unchanged.

Use Float32 or Float64 consistently, metadata vectors with matching scalar type,
and scalar `s`/`alphaT`/`scale` of that type. Filter `p` is an Integer at least one;
`alphaT >= 0` (including Inf), preserving the released filter contract. Numerical
values are not scanned in pointwise kernels. Transport owners establish finite
input magnitudes and a finite nonnegative gain before transforming; bound scales
and lossless input conversion follow the existing range helpers. The count in
`_fft_safe_magnitude` is positive **global** transform volume. Safe magnitude is
finite and nonnegative; zero is permitted when the bound underflows. Input
magnitude is finite/nonnegative. Copy/store scale is finite and at least one.

No collective communication is performed here. A backend must agree validation
failures and global magnitudes/scales across ranks before entering a transpose,
and use a normalized inverse before storing output. This extraction alone does
not establish distributed failure atomicity or multi-node scaling.

## Frozen validation and resource gates

Evidence directory: `HybridPlasmaPIC-audit-evidence/2026-10-26/shared-fourier-kernels-v1`.
Baseline metadata records every tracked source hash. Retain all failures.

1. Existing full test suite and Aqua on each runtime, before and after.
2. Deterministic Float32/64 serial families for dimensions 1/2/3, odd/even and
   singleton shapes: compare output, unchanged input, all four scratch arrays,
   metadata and repeated-operation histories bitwise against the frozen parent.
3. Independent dense DFT and analytic modes for derivatives, Laplacian, projection
   and filters; global-offset and permuted-storage tiles must equal full-kernel
   slices. Include DC, Nyquist, strict thirds, extreme k/field ranges, safe exact
   alias and rejected shifted overlap/invalid-axis before-write checks.
4. Ordinary warmed serial and raw-kernel allocations must be zero wherever parent
   is zero, with no additional persistent storage except returned metadata that
   the serial constructor already owned. No latency acceptance while shared heavy
   jobs run. A future paired cost gate must retain identical inputs and parents.
5. Concrete inference/JET plus complete raw package reports; review every changed
   report, do not suppress findings. No new dependency or broad dispatch.
6. Wrong local-index, odd/even Nyquist, reordered Laplacian, non-strict thirds and
   alias overwrite patches are rejected by the independent/mutation controls.

## Release and reproducible consumption

Current registry check (2026-10-29): General contains 0.1.3 through 0.1.6.
The original 0.1.3 release preparation and failed preflight below are historical.
The current scaled-source release gates live in the
[verification record](SCALED_FOURIER_VERIFICATION.md).

Prepare SpectralOperators 0.1.3 for the supported additions, with minimum/latest
Julia Linux/macOS CI added to the release gates. HPIC must raise its minimum
compat to `0.1.3` in a separate parent-owned integration commit. The accepted
source must be publicly available and General must contain its exact git tree
before clean registry resolution can pass. A local path override or immutable
git revision can support isolated integration tests; neither satisfies this
registry-release gate. Record `Pkg.status`, the manifest version/tree and the
qualified API checks from a fresh environment without a sources/develop override
on both runtimes before HPIC promotion.

The handoff prepares the exact candidate tree, an inventory of its files, and
the minimal General `Versions.toml` patch against pinned General commit
`a514710d144c744795f8eab8e4bfe7c8bec647ed`. Existing General dependency/compat
ranges cover 0.1.3. Root coordinates commit/push, exact-commit CI, actual registry
publication and tag/tree consistency. Preparation sends no registration-bot or
maintainer message. These were the original preparation limits. General now
contains 0.1.3; current release status must be checked independently of local tests.

## Progress

- [x] Inspect canonical owners, callers, contracts, tests and parent CAP-003 plan.
- [x] Agree minimal qualified API and pre-edit compatibility contract with root.
- [x] Freeze clean parent tracked sources and hashes.
- [x] Collect two-runtime serial/resource/compiler baseline: 6,284 existing assertions,
  238 contract checks, 448 state records and 224 zero-byte profiles per runtime.
  One full-package JET report per runtime retained for complete comparison.
- [x] Extract shared numerical owners and document supported API.
- [x] Verify raw logical-index/error/alias contract and both-runtime parity.
- [x] Review complete production/test diff, compiler reports and resource evidence.
- [x] Complete final evidence audit: 1,309 source/resource/artifact checks.
- [x] Prepare the 0.1.3 publication path and retained registry-failure evidence.
- [x] Shared primitives are registered in General starting at 0.1.3.
- [ ] Current release/consumer CI gates are tracked in the linked verification record.

Final V6 verification on **each** Julia 1.12.7 and 1.13.2: 7,620 full-package
assertions, 238 serial quality/resource checks and 2,245 complete-state replay
checks. All 448 state records match the parent bitwise. All 224 serial profiles
and 180 raw logical-view profiles have three zero-byte samples; the latter include
ordinary scale=1 and scale=2 copy/store on global-offset/permuted storage. All
372 concrete raw checks pass. Six constructor profiles with five allocation
samples each match parent bytes, retained storage and field layout exactly.
The one full-package compiler report per runtime is unchanged after a checked
constructor source-line mapping; none is suppressed. Fifteen canonical-owner
AST comparisons and fourteen actual source-mutation controls per runtime pass.
Production remains byte-identical to V4. Latency acceptance has not run.

## Retained findings

- Candidate V1: both runtimes passed 7,023 package assertions and failed one Aqua
  unbound-type-parameter assertion. Three new signatures admitted an empty
  metadata tuple that could leave `T` unbound. V2 narrows those metadata arguments
  to the already documented nonempty tuple contract and checks tuple dimension
  before mutation. No numerical loop or test limit changed. Raw V1 logs/source
  remain in the evidence directory; the serial comparison stage did not run after
  the full suite failed and is not claimed passed for V1.

- Candidate V2: 7,231 pass / one Aqua failure on each runtime. All 156 raw global
  view allocation profiles passed. The remaining projection parameter `D` was
  flagged unconstrained (`unbound-probe-v1.jl` records Julia's exact predicate).
  V3 derives `D` from the nonempty full metadata tuple and retains explicit
  dimension checks on every component before writing. A mismatched component
  dimension still rejects; no Aqua suppression or numerical change is used.
- Source review removed the obsolete private grid-only Laplacian guard adapter;
  all owned callers now call the sole `_laplacian_gain` owner directly. No
  external callers were found in sibling HPIC/FourierPSTD source. This avoids
  introducing an unnecessary compatibility wrapper for a private implementation.
- Historical registry preflight: the then-current General commit
  `a514710d144c744795f8eab8e4bfe7c8bec647ed` ends at 0.1.2. A clean project with
  minimum 0.1.3 fails resolution visibly on Julia 1.12.7 and 1.13.2. Local verification is
  not a release. This failure is retained as historical evidence; it is not the
  current registry state or a current consumer-resolution result.

- Candidate V3 failed to load on both runtimes: JuliaFormatter 1.0.62 changed
  `public _fourier_metadata` into the identifier `public_fourier_metadata` (and
  three other declarations). V3 source/logs are retained. A minimal syntax probe
  verifies JuliaFormatter 2.14.0 preserves the `public` AST; V4 restores the
  declarations and uses that formatter with an AST-preservation check before
  final verification. This was a tooling-induced candidate failure, not a parent
  package defect; no numerical test result is claimed for V3.

- Replay checker V1 completed 453 state/resource assertions, then errored because
  its mutation control indexed the gradient tuple instead of scalar output.
  V2 corrected the index but its `complex_scratch += 1` mutation rounded away for
  large fields (454 pass / one failed control). V3 uses `nextfloat` on the real
  component and verifies the mutation itself for every one of the 448 records.
  Original comparisons and limits remain unchanged; both failed checker versions
  and logs are retained. Production V4 source is unchanged by these corrections.
- Mutation harness V1 detected all five wrong-source variants (14 assertions per
  runtime) but emitted Julia world-age warnings while retrieving newly included
  module bindings. V2 uses `Base.invokelatest` at that dynamic loader boundary and
  passes the same 14 assertions with `--depwarn=error` on both runtimes.
- Full compiler auditor V1 had a Python string/compiled-regex API error after its
  source-line check. V2 corrects that auditor call; complete parent/candidate
  report text matches after one checked `plan = plan_fft!(cbuf)` location mapping
  per runtime. The sole parent report is retained, not suppressed or reclassified
  as a newly proven upstream defect.

- V4 raw concrete JET gate: 300 pass / 24 fail on each runtime. Every report
  traces Base broadcast's possible `unaliascopy` branch to missing `similar` for
  the private `FourierLogicalTile` test adapter (copy/store helpers only).
  `adapter-copy-baseline-v1.jl` reproduces six actual `copy(tile)` MethodErrors.
  V5 completes that adapter's `similar` interface and tests copied logical axes,
  values and independent storage. Production source remains byte-identical to V4;
  no correct broadcasting/range arithmetic is rewritten to suppress diagnostics.

- V5 passed 7,592 full-package, 238 serial quality, 2,245 complete-state replay,
  and 324 concrete raw quality assertions per runtime. All 156 raw profiles had
  zero reports and three zero-byte samples. Constructor receipt V1 then failed
  because TOML does not serialize a tuple of field names; V2 writes an array and
  retains the failed receipt. Six constructor profiles with five samples each
  match parent allocations/storage/layout exactly on both runtimes.
- Final audit expands V6 raw copy/store coverage from scale=2 to both scale=1
  (the ordinary broadcast branch) and scale=2. This adds 24 resource/compiler
  profiles under the same zero-byte/no-report gates and four value assertions.
  Production source remains byte-identical to V4. No numerical tolerance or
  resource gate is relaxed.
