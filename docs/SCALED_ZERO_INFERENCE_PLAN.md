# Canonical scaled-buffer zero inference repair: implementation tracker

Current status (2026-10-29): the helper and both zero owners, together with the
Cartesian store repair, pass permanent tests and the full both-runtime normal
campaign. Independent root postflight accepts bounded local gates; the
[verification record](SCALED_FOURIER_VERIFICATION.md) owns actual results and
remaining CI/registration/HPIC consumer work. All failed prior candidates and
their full captures remain preserved.

The source-supported pre-edit brief below is preserved as the reviewed design;
its future tense is design history. Root owns global trackers and publication.

# Canonical scaled-buffer zero inference repair: pre-edit brief V1

Historical pre-edit status (2026-10-29): design recorded before production edits. Root authorized a
new isolated worktree from bb9d01df2789f8729e680cf36d15a3cefd15351d. Old
595/bb9 source trees, all 937 historical pins, failed normal campaign and raw
compiler captures remain frozen. No normal campaign, publication or HPIC36 release.

## Verified finding and bounded evidence

The failed normal candidate retained two reports for the actual unbounded public
Laplacian. The complete cold diagnosis reproduced both reports at inherited
zero-fill owners: `src/scaled_differential_operators.jl:176,199` in bb9. Generic
`g.cbuf` is `Array{Complex{T},D} where {D,T}`; generic `eltype` loses scalar
correlation in the captured analysis. Concrete supported calls did not reproduce
a runtime failure. This is an inference repair, not a claim of invalid runtime
indexing, a cache cause, or a new physics defect. Direct typed zero plus unchanged
rank-generic fill removed one report but retained one; that failed proposal is
preserved in `zero-buffer-prototype-v1`.

The separately released normal-compilation Julia1.12.7 alias experiment completed
2504 assertions: both canonical Base.vec and Base.view alternatives had empty
full unbounded buffer/grid captures. All 840 individual warmed allocation samples
across 42 profiles were exactly0B, with exact owned storage and 364 private semantic
rows plus 32 grid-state rows. The actual unchanged public entry still had two
reports. At that prototype stage, latest Julia and integrated public/package quality were untested; the final campaign is recorded above.
`zero-buffer-alias-v1/actual-results-v2.json` SHA
8953251519a671ab4b5a98ca201d80d3175d754c7faf5ba3d60ac0a08a5640fc
retains all full reports and five literal historical/current public render deltas;
no normalization or literal-equality claim is made. Root's bounded alias review is
`so017-zero-buffer-alias-actual-root-review-v1.json` SHA
173c91618f68ef42491f2d4f02785e3f65d7d7ac0ac6b1dcb330d756ea01a685.

## Canonical owner, source proof and smallest change

`FourierGrid.cbuf` is an Array{Complex{T},D}
(`src/periodic_fourier.jl:25–39`). Julia1.12.7 and1.13.3 Base vec is documented to
share storage; AbstractVector returns itself, other arrays call reshape with
length (`base/abstractarraymath.jl:12–17,41–42`). Array reshape creates a rank-one
header using the original data reference after checking equal length
(`base/reshapedarray.jl:41–50`). It does not change the original header/shape.
Thus `fill!(vec(buffer), zero(Complex{T}))` writes exactly the same zero value to
all elements, including rank0 and empty arrays; returning buffer preserves the
original object, type, rank and dimensions. The alias is ephemeral and never
retained or resized. Private Float16/BigFloat/Int semantic comparisons passed;
this does not expand public FFTW admission.

Add exactly one private `_zero_scaled_buffer!(buffer::Array{Complex{T},D}) where
{T,D}` using canonical `fill!(vec(buffer), zero(Complex{T}))`, then `return buffer`.
Reuse it only at the two existing zero-fill sites: inactive normalized multiplier
and `_zero_scaled_data!`. Select vec because its checked Base source is simpler
than view; no timing superiority claim. No public signatures/admission, fields,
FFT count, operator arithmetic, thresholds, ordinary path or resource caps change.
No handwritten fill loops, unsafe memory, cached alias, supported-only analysis,
new fallback or analyzer configuration changes. Removing the helper and restoring
the two calls must reproduce bb9 executable source exactly.

## Permanent proof and fixed gates before measurements

- Permanent tests compare unchanged Base.fill! reference values, positive complex
  zero bits, original return identity/type/rank/shape and same storage for real
  scalar Float32/64 rank0–4, including empty/singleton/nontrivial shapes. Retain
  private non-FFTW scalar compatibility tests. Exercise both canonical callers,
  inactive modes, null axes and public zero/constant/tiny inputs with exact plan
  counts, metadata/neighbor/input preservation and borrowed lifetime checks.
  Existing unsupported manual-grid rejection fixtures remain unchanged.
- Keep actual default JET full unbounded public Laplacian hard0, actual complete
  report_package, all60 concrete captures and both floor entries. Add full
  unbounded zero-owner and both-caller captures; do not replace the public gate
  with supported-only analysis. Every full frame/render remains retained.
- Preserve full-package 595→candidate report comparison. The existing constructor
  tilde-root spelling issue is an independent evidence/comparator defect: a new
  exact root correction requires source-specific proof and negative controls.
  Five literal unrelated Bool-render deltas remain raw and are not retrofitted.
- Existing342 API resource rows, 380 floor profiles, 432 bit profiles and seven
  floor mutants keep their exact gates. Add actual canonical zero-owner profiles
  using the already declared 42-shape/grid matrix, all individual warmed samples
  exactly0B, exact owned storage. No timing or full-native-peak acceptance.
- Source/AST/inverse checks, formatter1.0.62, locked dependency/provider preparation,
  full min1.12.7/latest1.13.3 package/Aqua, resources, compatibility and complete
  compiler captures require a fresh reviewed paired campaign. Preserve old failed
  runs and unexecuted plans; no retry/overwrite. Source/runtime/dependency/native
  image identity remains checked, with only explicitly proved package path/cache
  changes. No inferred cache/render cause.
- A declared source/stdlib preflight will retain actual ===, == and both subtype
  results plus raw renders for the two observed UnionAll expressions. It cannot
  replace the actual integrated public hard0 gate or erase historical differences.

## Tracker

- [x] Source-supported cause classification and bounded normal alias experiment.
- [x] Seal full evidence, raw literal differences and this pre-edit brief.
- [x] Isolated helper/two-call patch, store repair and permanent tests.
- [x] Static inverse/AST/caller checks, formatter and locked preparation.
- [x] Reviewed paired normal campaign and complete independent root postflight.
- [ ] Publication/CI/registration and separate HPIC36 integration.
