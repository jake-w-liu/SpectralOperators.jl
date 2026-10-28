# S2 source-aware range repair — pre-edit brief and implementation tracker

Sealed design and verification record (2026-10-29): the final scaled-source implementation passes both
runtime normal campaigns and independent root postflight. See the
[verification record](SCALED_FOURIER_VERIFICATION.md) for actual gates, release
provenance and the separate consumer boundary. The original clean595 package had two extra abstract
floor reports; subsequent floor, zero-buffer and Cartesian store repairs are
integrated. All failed intermediate source and evidence remain preserved.

The approved pre-edit design follows as a historical contract. Its external
V3 brief (e241a456...) and233-pin evidence receipt remain sealed; the frozen
HPIC eab73e8, original S2 4dd7243 and SO e92f326 are unchanged. V1/V2 designs and failed evidence
harnesses remain retained. This document supersedes V2's stored-E Faraday replay,
63-FFT estimate and 4 MiB/cell allocation proposal. Root owns main integration,
dependency version allocation, release/registration and exact-SHA GitHub CI.

## Pre-Edit Brief

- **Mode:** repair plus extension of the canonical Fourier numerical representation.
- **Requirement:** retain representable coefficient-scaled Laplacian fields and
  complete pressure/Faraday rates even when a separately rounded intermediate
  would disappear. Preserve the full-k Laplacian, centered odd derivatives and
  positive Nyquist-completed heating. Reject an exactly nonzero complete required
  physical RHS norm that has no nonzero storage-type representative. Preserve
  ordinary arithmetic and etaH=0 behavior, public force/J/gradp/lapJ meaning,
  transactions and restart configuration.
- **Baseline:** primary public steps lose all magnetic damping while adding
  positive heat (approximately 2% common-time energy increase); downstream
  dt-rescued witnesses lose complete physical source norms; a separate constant-E
  witness loses varying hyperresistive induction with exactly canceled total
  pressure. Source and raw-buffer packets are listed below.
- **Cause:** the actual canonical ordered Fourier multiplier loses nonzero
  products before etaH is applied. Additional rounding can occur at E publication
  or physical RHS publication. Differentiating rounded E cannot recover a source
  that was already lost there. This was reproduced, not inferred from small norms.
- **Invariants:** normalized data never enter public physical buffers; complete
  signed scalar sums precede final range classification; exact cancellation is
  accepted; a small intermediate is a routing condition, not a rejection proof.
  Exact scalar composition refers to represented canonical derivative data,
  not a claim of symbolic/exact FFT coefficients.
- **Owner:** SO owns Fourier input admission, centered/full metadata, normalized
  multiplier/FFT representation. HPIC Ohm/Faraday owns complete transport field
  and induction composition; the pressure equation owns all heat/advection terms.
  The integrator coordinates those owners and retains existing RK storage/update
  and outer transaction. No second integrator physics formula.
- **Compatibility:** etaH=0 takes the unchanged S1 path. Ordinary positive-etaH
  inputs retain old operation order/bits and allocations. Cold evaluation may
  differ from curl of separately rounded E: it completes the constitutive/operator
  expression before E publication. This is an explicit range contract, not an
  added correctly-rounded-RHS guarantee for every ordinary cancellation.
- **Strategy:** borrowed normalized Fourier channels, exact scalar composition
  through existing rational/conversion owners, one streamed rational accumulator
  array and three transient ideal-base arrays. No new persistent storage/schema.
- **Verification:** fixed numerical/resource gates below; both 1.12.7/1.13.3,
  independent DFT/SBP/energy/refinement, rollback/retry/restart, ordinary exact
  comparison, mutants, Aqua/concrete JET/full compiler frames, formatter/AST.

## Verified evidence and limits

All paths below are under `underflow-public-v1` unless specified otherwise.

| Evidence | Observed result |
| --- | --- |
| `candidate-112`, `candidate-1133`, `parent-1133` | Stable public primary cases: B unchanged, positive heat, about 2% common-time energy rise. Full represented state/term buffers retained. |
| `canonical-scale-v1` | Actual nonzero hats become all zero inside canonical full-k multiplication. Eight assertions/runtime. |
| `block-prototype-commands-v2.json` | Per-mode/per-axis normalized Laplacian retains 108 ordinary profiles/324 checks/runtime and four extreme cases. |
| `block-extremes-commands-v1.json` | Anisotropic inactive stiff axis and broad input span: global normalization error=1; per-mode result agrees with independent 8192-bit DFT within 4.34e-8 F32 /8.48e-17 F64. |
| `block-gradp-commands-v2.json` | Positive pressure gradient cancels separately overflowing hyper product to finite E, preserving J/gradp; 12 assertions/runtime. Complete Faraday overflows in this raw-Ohm case and is not claimed admissible for stepping. |
| `block-downstream-underflow-commands-v1.json` | Rescuing E alone still leaves pressure/Faraday rates all zero and misses .00995017 relative decay; eight pass/two expected failures/runtime. |
| `block-rhs-commands-v3.json` | Power-up before the cold forward FFT improves subnormal-source accuracy; complete source/cancellation classification passes 32 checks/runtime. |
| `source-aware-boundary-commands-v4.json` | Prepared total pressure is exactly zero, E3 is exactly constant one, stored-E curl target is zero, but public By misses .01309707 decay. Both below-minimum and normal hidden-source variants; 24 pass/four expected failures/runtime. Independent midpoint control avoids assuming deposition round trips are exact. |
| `source-aware-certificate-commands-v2.json` | Source-aware linear completion gives represented target rate; exact dyadic modal reference error zero; ten checks/runtime. |
| `source-aware-public-after-commands-v1.json` | Evidence-only stage replay repairs admitted normal hidden-source public steps; 14 checks/runtime, relative By error 2.31e-8 F32 /zero F64. |
| `centered-kernel-bound-commands-v2.json` | 2,742 checks/runtime: actual odd/even stored centered/full metadata, n1–65, largest/unit/short domains, 8192-bit circulant row sums and metadata rounding. |
| `warm-routing-controls-commands-v2.json` | 852 checks/runtime: ordinary controls, inclusive threshold neighbors, independent canceled-pressure/constant-E routing and multiplier normal-headroom neighbors. |
| `complete-field-scalar-controls-commands-v1.json` | Complete ideal/Ohmic/hyper overflow-cancellation, exact zero and minsubnormal neighbors against 8192-bit references. Scalar-owner evidence; public Fourier composition remains a permanent implementation gate. |
| `source-budget-derivation-v3.json` | Finite-representation bit envelope and explicit static operation/FFT model; this is not a measured implementation or GMP allocation theorem. |

Failed harness revisions remain available. In particular centered-kernel V1
completed 2,652 checks and reported 90 undefined `isnormal` harness errors; V2
uses the explicit finite/abs>=floatmin condition. Old boundary V1/V3 assumptions
about manual midpoint moments/deposition and vanishing *other* curl components
were corrected and retained. None is presented as a production failure.

The dt-rescued **unrepresentable required RHS** cases remain unsupported in this
slice. Their finite update would need normalized dt-scaled RK/pressure storage,
a separately tracked implementation. Do not claim full exponent-range physics.

## Canonical interfaces and state boundaries

SO e92f326 sources checked: `periodic_fourier.jl` metadata, finite input scaling,
loss-of-nonzero admission and output restoration; `differential_operators.jl`
first derivative/gradient/curl, full-k Laplacian and real Nyquist channel;
`docs/FOURIER_BACKEND_CONTRACT.md` qualified backend and one-caller scratch.

Add public, unexported qualified operations:

- `laplacian_scaled_data!(g, f)`
- `deriv_scaled_data!(g, f, axis)`
- `nyquist_scaled_data!(g, f, axis)`
- `deriv_laplacian_scaled_data!(g, f, axis)`

Each returns the same concrete NamedTuple: `values` (borrowed complex grid
buffer), `input_scale::T`, `exponent::Int`, `ordinary::Bool`. Its physical real
value is `real(values)*input_scale*2^exponent`. The caller must consume the data
before another grid operator. Validate matching type/shape/one-based axes and
input disjointness from every mutated scratch/metadata before writes. No new
thread safety, retained borrowed result, GPU or MPI claim. Plain serial methods
keep their existing interfaces and arithmetic. Metadata remains readonly and
must satisfy existing full squared-wavenumber admission.

HPIC eab73e8 sources checked: shared `_add_hyperresistive_field!`; canonical
`_ohm_Efield!` and `faraday_rhs!`; exact/conversion owners; pressure rate, heating
coefficient `(gamma-one(T))` convention and positive channels; stage ordering;
outer init/step rollback and checkpoint configuration. The raw Ohm function may
leave partial work on later error. Coupled owners retain their current rollback
payload; derived scratch is not newly checkpoint-owned. Final carried ideal
force/J/pressure-gradient state must be tested explicitly.

The canonical ideal-base evaluation remains the existing Ohm owner with eta=etaH=0,
including its finite represented input/output requirements. This slice does not
replace nonlinear ideal-force arithmetic, ion velocity, current, public pressure
or the existing inverse-density owner with arbitrary precision. These domains
must remain explicit rather than silently clipping or reconstructing force by
subtracting resistive fields.

## Normalized operators and ordinary admission

For nonzero real/imaginary hat component x and metadata k, use exact binary
normalization x=mx*2^ex, k=mk*2^ek, with magnitudes of mx,mk in [1,2).
Choose one common integer P as the maximum active product exponent. Compute each
axis contribution with normalized products followed by `ldexp` and the original
axis order. Laplacian uses `-((mx*mk_full)*mk_full)`; derivative uses mx*mk_center
and canonical i rotation; derivative-Laplacian uses the ordered full-k square
followed by the centered derivative. No squared tiny global ratios. Nyquist uses
its real full-k factor only on its plane, and returns exact zero on odd axes.

The ideal normalized spectral scalar bounds are 4, 8D and 16D respectively.
For Int64 volume at most 63 axes are active: 4,504,1008. Root independently
checked these against unit-gain FFT headroom for Int64/Int32 volumes, including
the omission bound. This does not prove FFTW rounding. Normalized multiplier
roundoff is bounded by gamma_(D+3) at its positive operation scale; with D<=63,
the local gate is gamma_(D+3) times the sum of exact positive normalized term magnitudes, plus the explicit omission envelope. In physical D1–3 this is a few rounding units, not permission to spend a blanket global epsilon multiple. Per-axis ldexp loss is at most half-minsubnormal; the
volume-scaled omitted-mode envelope is below eps/512. The runtime must validate
all consumed normalized values are finite.

For cold input copying, use canonical **unit-gain** headroom and loss-of-nonzero
rejection. When its maximum is below one, an exact power-of-two copy raises it
to [1,2); carry the inverse power in the returned exponent. Do not normalize a
broad large input down to one. Physical-gain input scaling must not discard a
nonzero input before the normalized branch is selected. Reuse the canonical
loss predicate/scale policy rather than catching an error and retrying blindly.
A maximum of one discarded forward attempt/component is allowed and counted.

The ordinary Laplacian branch uses identical old input scale, FFT, ordered
multiplier, inverse and restoration. It is selected only after normal headroom
checks. With emin=exponent(floatmin(T)), emax=exponent(floatmax(T)),
L=ceil(log2(volume)), A=ceil(log2(max(1,active axes))), every nonzero component/
axis term must satisfy ex+ek>=emin and

    emin+L+4 <= ex+2ek <= emax-A-4.

The lower margin is sixteen times the inverse-volume normal floor; the upper
margin is sixteen times the axis-sum headroom. A tiny nonzero input maximum
below `ldexp(floatmin(T),L+4)` also selects cold before the forward FFT. These
are conservative sufficient **routing** conditions, not admission failures and
not a relative bound on canceled real-space cells. Whole-result FFT cancellation
continues to have the standard declared source/operator-scale accuracy contract.

## Source-aware stage completion and routing

The stage retains the existing canonical ideal-base computation and current.
One transport field mutation owner receives ideal E, eta, J, etaH and normalized Laplacian data together, and reports scalar allocation-free range information. It evaluates the exact old finite order: first ideal+eta*J (or the eta=0 copy), then subtraction of etaH*lapJ. If either intermediate loses range, its cold expression is the complete ideal+eta*J-etaH*LJ sum through the existing exact owners. In particular it must not accept an already rounded/overflowing Ohmic initial field as its cold source. The old hyper-only helper delegates with eta=0; there is one field formula owner. Raw E may alias its ideal input cellwise; the carried force retains its separate owner.
A normalized Laplacian is restored independently to physical lapJ; nonfinite
physical lapJ rejects. Cold scalar E uses existing `_exact_rational`/`_from_exact`
with the complete ideal/Ohmic/normalized-hyper source, never physical rounded lapJ. Public E may be zero
although a derivative of the complete source remains representable.

Routing checks happen before an unsafe ordinary heating/spectral operation can
raise from an unrepresentable intermediate, and after ordinary complete rates
when no earlier risk exists. They may request one bounded cold replay; no catch
or masked DomainError. A conservative preflight uses canonical input headroom,
full/centered metadata and positive source bounds to route physical derivative
range risk before heat's exponent extraction. For each nonconstant input, compute
the upward physical channel bound Gd*max|input| (kNyq*max|input| for a Nyquist
channel). Route if a nonzero bound lies below 32floatmin(T), exceeds
floatmax(T)/32, or is nonfinite. A cold request also precedes any ordinary
velocity-divergence copy that needs nonunit FFT input scaling, and any tiny
input maximum below the inverse-volume normal floor. These checks are
conservative sufficient routing conditions; exact constants/null axes are
handled without treating a zero upper bound as a failure. Nonfinite complete ordinary rates
route before their usual finite rejection. The cold path then either publishes
the represented complete rate or rejects it explicitly.

When etaH>0 and current is not exactly componentwise constant, route cold on:

1. normalized-Laplacian/replay use or nonnormal/nonfinite etaH-Laplacian products;
2. ordinary derivative input/multiplier or heating intermediate range risk;
3. pressure-rate norm below 32floatmin(T), including zero;
4. magnetic-rate norm at most 2^-10 G*S.

S is an upward bound on max|Eideal|+eta*max|J|+etaH*max|LJ|, formed before lossy E
publication (use the normalized source when needed). Positive Float64 upward
arithmetic may produce Inf, which selects cold. It must never clamp or turn a
nonzero product into a zero bound. G=sum of actual centered-derivative row-sum
bounds. For axis n<=2, Gd=0. For n>2 and stored fundamental k0,

    Gd <= k0*(n/2)*(1+ceil(log2(floor(n/2))))
          + gamma2*k0*n^2/4.

The first term bounds both the even cotangent and odd cosecant differentiation
matrix; the second bounds generated integer conversion/multiplication metadata
rounding. All positive implementation arithmetic is rounded upward. The actual
kernel bound and odd/even controls are independently checked in the receipt
above. Primary sources: Oxford's Trefethen p2.m and Chebfun trigcolloc/diffmat.m,
whose downloaded artifacts/URL/hash receipts are retained.

The factors 32, sixteen and 2^-10 are conservative engineering routing choices,
not newly asserted mathematical admission boundaries or error allowances.
Neighbor tests freeze inclusive comparisons; they do not establish a uniform
relative bound on arbitrarily canceled RHS. All routing depends on complete
source information, independently of pressure failure. A large constant field
cannot certify a hidden varying hyper source. No allzero/E-norm rejection rule.

Cold Faraday evaluates, from the same stage data,

    -curl(Eideal) - eta*curl(J) + etaH*curl(LJ),

using canonical scaled derivative/derivative-Laplacian channels and a complete
signed exact scalar sum before conversion. Recompute ideal E through the same
Ohm owner into three transient T arrays; do not copy warm fields every stage or
reconstruct force by subtraction. Confirm current recomputation is bitwise equal.
Pressure exact composition includes advection, compression, Ohmic heat and every
positive gradient/Nyquist channel, in its canonical equation owner.

Stream each borrowed channel immediately into one Rational{BigInt} N-array.
Process pressure and magnetic components sequentially. Do not retain rational
channel arrays. Keep complete pressure and full three-component magnetic norms
for classification: exact zero passes; a nonzero complete norm with all rounded
values zero rejects; nonfinite final values reject. Individual zero cells/channels
are valid and never sufficient for rejection. Public init/step must be atomic.

## Fixed numerical and resource gates

Accuracy is normwise at the explicitly declared canonical operation/source
scale, not relative to an arbitrarily canceled small RHS. Before measurement,
use 4096eps solely as the complete stage/physics-oracle engineering gate.
Enforce and report the tighter local normalized-multiplier gate
gamma_(D+3)*sum|exact terms| plus D*half-minsubnormal omission, including its
normalization/inverse-volume conversion when comparing physical units.
For accounting only, 128eps bounds this multiplier portion even at D=63;
it does not replace the local dimensioned gate. The empirical FFT comparison
portion is at most2048eps; restoration/conversion and ordinary local constitutive
arithmetic each have their own at most8eps operation-scale gates. The maximum
operator/RHS accounting sum is2192eps; the remaining1904eps is reserved for
stage/time-composition oracle error. Report observed worst local, source, FFT
and time/refinement errors separately. Do not spend the FFT portion again on local arithmetic. The FFT
portion is an empirical acceptance gate against independent high-precision DFT,
not a proof of FFTW error for every transform size. Local ordinary expressions
retain their existing few-ulp source-scale behavior; cold scalar sums/conversion
use canonical exact arithmetic. Refinement/energy tests separately check the
actual RK order and paired heat/field conservation.

Warm ordinary allocations, owned storage/state and outputs must match eab73e8
exactly; etaH=0 replay is bitwise. Retain 256 KiB/cell field-only cold total gate.
For complete cold RHS, static preimplementation accounting is in
`source-budget-derivation-v3.json`: numerator <19,958 bits and denominator <32,768
bits under finite consumed data and admitted canonical metadata; planned112
rational multiplications,51 additions,60 finite conversions,10 final conversions
per cell, including the complete three-term cold field evaluation. Proposed total allocation gate is 8 MiB*N+2 MiB fixed. It is an
engineering budget, **not a GMP allocation theorem**. The model predicts
6,057,984 B/cell before measurement. Record actual operation counts as implemented.

Report peak-live bytes separately: one rational N-array bounded by a proposed
16 KiB/cell retained-object gate, three T arrays (at most24 B/cell), plus a
1 MiB scalar transient allowance; existing FFT grid storage is reported separately.
Use separate total-allocation and peak-live measurements. No RSS/timing ceiling
claim under host contention, and no claim of full optimization from these caps.

At most 97 additional FFT executions per cold stage: 15 first-derivative and9
Nyquist pressure calls, six ideal/six Ohmic/six derivative-Laplacian Faraday calls
(each two transforms), up to three discarded Laplacian forward attempts, plus
4 pressure-gradient and6 current recomputation transforms. eta=0 omits12. Tests
must count actual paths/retries, not just a symbolic happy path. If implementation
requires more calls/retention/bit operations, revise the declared plan before
measurement rather than widening a failed gate afterward.

## Implementation and verification tracker

This cross-package breakdown retains the original design boundary. The HPIC
checkboxes describe acceptance scope; current implementation progress belongs to
[the consumer capability plan](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/blob/main/docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-011--electron-transport-and-complete-energy-accounting).
Dependency completion and release evidence are recorded above.

- [x] Preserve clean overflow candidate and failing/public typed evidence.
- [x] Verify primary field loss, anisotropic/full-Nyquist and broad-span failures.
- [x] Verify downstream required-rate failure and exact-cancellation distinction.
- [x] Verify source-aware normal hidden-source public-step feasibility.
- [x] Derive and independently check actual odd/even metadata kernel bound.
- [x] Replace 63-FFT/4 MiB assumptions with explicit streaming/97-FFT/8 MiB model.
- [x] Finish routing-neighbor and mixed-source scalar controls; preserve failed harnesses.
- [x] Freeze final source/evidence pin inventory before production edits (`underflow-design-evidence-v3.json`).
- [x] Root design review/version allocation and isolated HPIC/SO worktrees.
- [x] Implement and locally verify SO scaled-data owners and permanent DFT/symbol/alias/range tests.
- [ ] Implement canonical field/pressure/Faraday completion and stage coordination.
- [ ] Permanent primary energy, hidden normal source/canceled pressure, exact zero,
      subnormal neighbors, unsupported dt-rescued rollback/retry/restart regressions.
- [ ] Mixed-axis SBP, Nyquist, cancellation, constants, input/output range and
      final carried force/J/gradp checks; keep public scratch physically unscaled.
- [ ] Mutants: coefficient after rounded lapJ, global k ratio, centered Laplacian,
      omitted Nyquist heat, rounded-E-only replay, incomplete pressure sum,
      blanket zero rejection and normalized public scratch.
- [x] SO source/harness freeze, both-runtime normal numerical/resource/quality campaign,
      complete compiler disposition, formatter/AST/idempotence and root postflight.
- [ ] Complete the separate HPIC composition campaign under its unchanged gates.
- [ ] Root current-main composition; push/CI/registration; update global tracker
      with unsupported normalized-dt future work and actual supported boundary.

## Historical initial SO implementation receipt

Candidate base: e92f326, isolated branch feature/scaled-fourier-sources-20261027;
version0.1.7 reserved by root. Scope here is the four shared numerical owners;
HPIC composition follows separately. For the new API, the input must be disjoint
from mutated FFT scratch; cache writes must be disjoint from readonly metadata.
A readonly input that happens to be metadata need not be rejected unless it
aliases mutated scratch. This clarifies the brief's compact alias wording.

- [x] Canonical sources, public contract, parent state and pre-edit evidence read.
- [x] Add canonical scalar loss predicate reuse and four scaled Fourier data owners.
- [x] Add independent8192-bit local multiplier/DFT/SBP, ordinary bits, range/alias/lifetime and actual FFT-count tests;12,098 assertions/runtime.
- [x] Complete legacy source AST comparison11/11 under exactly checked scalar predicate lift and new public/include additions.
- [x] Execute both runtime normal focused/full suites and capture raw allocation/compiler results for clean595.
- [x] Review all final resources, actual unrestricted/public captures and complete package reports; see the verification record.

Implementation evidence: `scaled-so-v017/focused-light-commands-v2.json` has
two exit0 commands with unchanged source freezes; local error <=0.352941 of
gamma+omission and DFT source error <=0.877973eps. `legacy-ast-v1` passes11/11.
Pinned JuliaFormatter1.0.62 plain whole-tree attempt failed because it corrupts
Julia `public` declarations; the AST guard caught it. The failed log and exact
interrupted snapshot are retained. Public-aware selected-file formatting keeps
those exact lines, preserves the full AST and passes two idempotent checks
(18 assertions); no production formatting exclusion was introduced.

Final local gates and root postflight are recorded in the linked verification record.
All four dependency CI jobs, registration, tag and release are complete; the
verification record binds their exact source. HPIC integration remains separate
and incomplete, with current progress owned by its capability plan.
