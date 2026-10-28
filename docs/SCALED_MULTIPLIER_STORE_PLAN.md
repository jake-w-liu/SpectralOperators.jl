# Normalized Fourier multiplier store boundary

Status (2026-10-29): complete both-runtime local gates and independent root review pass. The [verification record](SCALED_FOURIER_VERIFICATION.md) owns results and remaining publication, CI, registration and HPIC consumer gates.

`FourierGrid.cbuf` and the `CartesianIndices` created from that same array have the same rank. The private scalar store retains that relationship in its argument types, then delegates conversion and bounds checking to the original Base assignment. The multiplier keeps its Cartesian traversal and complete scalar arithmetic. There is no new public API, retained allocation, unchecked cast, bounds suppression or fallback.

The previous canonical candidate failed five unrestricted JET zero-report assertions. A vec-loop experiment cleared four entries but allocated 32 B in 324 profiles; direct and typed-OneTo loops failed four compiler assertions before later stages. These failures remain recorded. A separate 72-entry/144-vector IR observation shows unrestricted type erasure; it does not demonstrate a runtime Bool value or an allocation cause.

The new-name Cartesian-store experiment on Julia 1.12.7 clears all four unrestricted multiplier analyses and preserves all 1,200 output/state/storage profiles. All 18,000 individual allocation samples are 0 B; all 48 nonfinite error-state controls agree, with all 72 full IR entries retained. This is bounded experimental evidence. The canonical candidate subsequently passed the complete normal campaign on Julia 1.12.7 and 1.13.3; root postflight independently reviewed all results.

- [x] Read canonical owner, public callers and Base Array/CartesianIndex contracts.
- [x] Preserve failed prototypes and independently review all bounded outputs/resources.
- [x] Write the one canonical store path from the reviewed experiment.
- [x] Complete source inverse, formatter/idempotence and locked environment preparation.
- [x] Complete both-runtime provider/native-image and normal scientific/compiler/resource/package gates.
- [x] Review complete reports, exact final source and compatibility documentation.
- [ ] Commit/push, exact-SHA CI, Julia General registration and downstream integration.

[Zero-owner prerequisite](SCALED_ZERO_INFERENCE_PLAN.md) and [scaled-source contract](SCALED_FOURIER_SOURCE_PLAN.md) remain separate records. A successful primitive or package run alone does not accept HPIC hyperresistive composition.
