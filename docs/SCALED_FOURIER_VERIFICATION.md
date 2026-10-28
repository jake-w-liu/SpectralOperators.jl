# Scaled Fourier source verification

Checked 2026-10-29. The tested source is `816e0c7ff0a0156e781f14d45404c712c146d088`.
The publication composition preserves every tested non-Markdown byte and mode.
This record accepts bounded local package gates; exact-commit GitHub CI,
registration and HPIC consumer integration remain required.

The four borrowed operators retain coefficient-scaled source data through
intermediate range loss. Their supported behavior and lifetime are defined in
[the backend contract](FOURIER_BACKEND_CONTRACT.md#borrowed-scaled-fourier-data).
The floor, zero-buffer and Cartesian store repairs preserve supported numerical
arithmetic and checked Base conversion/bounds behavior.

| Gate | Actual result on each Julia 1.12.7 and 1.13.3 |
| --- | --- |
| Full package and Aqua | 26,886 assertions pass |
| Complete package compiler capture | Three reports in the original scaled-source candidate, one in the final candidate; the entire remaining constructor report matches under checked source locations |
| Concrete and unrestricted added captures | All 60 existing concrete, three floor/public and six zero-owner/multiplier captures are empty |
| Existing warmed API profiles | All 342 rows are exactly equal to the paired parent and allocate 0 Julia bytes |
| Floor helper profiles | All 380 allocate 0 Julia bytes |
| Complete output-bit profiles | All 432 match the paired parent |
| Zero-owner and inactive multiplier samples | All 840 and 320 individual samples allocate 0 Julia bytes; owned storage is unchanged |
| Mutation controls | All seven numerical/floor and five zero-owner mutants detected |
| Providers and native images | Complete graph/module/library comparisons, cache headers, embedded source and CRC proofs pass; 96 minimum / 95 latest negative controls reject |

Root postflight checked all 1,796 frozen inputs and 260 generated outputs,
all 20 serial normal-call receipts, and every provider/image/comparison side
receipt. It independently replayed the complete result and provider comparisons.
Both full compiler trees and all eight deltas per runtime remain retained;
no report, frame or message was filtered. The existing constructor report remains
an observed diagnostic, not a compiler-clean-package claim.

Evidence owner: `HybridPlasmaPIC-audit-evidence/2026-10-27/takeover/`.
The final receipt `so017-cartesian-normal20-actual-root-review-v1.json` has SHA-256
`d9fcc14a17d9c25a0211ec597efaf1a05854f5423523d9fa04b8f98b42600e9d`.
Earlier failed floor/zero/loop proposals and failed harnesses stay failed in their
immutable evidence namespaces. They are not results for the final candidate.

- [x] Canonical source/caller review, exact source inverse and formatter checks.
- [x] Both-runtime normal numerical, compatibility, resource and complete quality gates.
- [x] Independent complete root postflight.
- [ ] Publication commit and four exact-SHA Linux/macOS GitHub jobs.
- [ ] Julia General registration, tag and release.
- [ ] HPIC consumer composition and its separate physical/compiler/resource/CI gates.

Native/GMP peak accounting and exclusive latency remain unaccepted. The Julia
allocation results do not establish either. Full exponent-range coupled physics,
normalized-dt RK storage, distributed/GPU support and simulator completeness are
outside this package slice.
