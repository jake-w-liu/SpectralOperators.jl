# Scaled Fourier verification and release record

Sealed record, 2026-10-29. SpectralOperators 0.1.7 is registered, tagged and
released from commit `600cd262dd1d5fa71b378448713ab2f50c404bb3`, tree
`9074d7ef57054e4772f45134125eb4ec7125ef79`. Its 20 non-Markdown files preserve
every tested byte and mode of `816e0c7ff0a0156e781f14d45404c712c146d088`.
All four Linux/macOS minimum/latest jobs pass in the
[exact-revision CI run](https://github.com/jake-w-liu/SpectralOperators.jl/actions/runs/40784524878);
each full raw log records 26,886 passing assertions on Julia 1.12.7 or 1.13.3.
CI runs the package and Aqua suites. The complete compiler disposition below
comes from the source-identical paired local campaign.

Read [General's version/tree records](https://github.com/JuliaRegistries/General/blob/master/S/SpectralOperators/Versions.toml),
[GitHub releases](https://github.com/jake-w-liu/SpectralOperators.jl/releases) and
[CI](https://github.com/jake-w-liu/SpectralOperators.jl/actions) for later state.
This page records checked releases and tests; it does not maintain a live queue.

## Release provenance

| Release | Included change |
| --- | --- |
| [0.1.3](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.3) | Qualified reusable Fourier kernels; its missing annotated tag/release was restored against the exact already registered source |
| [0.1.4](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.4) | Infinite converted-strength exponential filter |
| [0.1.5](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.5) | Full-operation inverse-FFT headroom |
| [0.1.6](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.6) | Full-Laplacian Nyquist energy channel |
| [0.1.7](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.7) | Borrowed normalized source channels and canonical floor, zero-buffer and Cartesian store repairs |

The 0.1.7 [General request](https://github.com/JuliaRegistries/General/pull/173877)
merged on 2026-10-28 UTC. The annotated tag dereferences to the registered commit;
the public release and complete release notes were read back. Root receipt
`so017-registration-tag-release-actual-root-review-v1.json` has SHA-256
`5085411258e2d6038c7a3111425879c900f047a18853bf0f46b5c63aeaab3fb3`.

## Bounded local verification

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
- [x] Publication commit and four exact-SHA Linux/macOS GitHub jobs.
- [x] Julia General registration, tag and release.

HPIC consumer composition and its separate physical/compiler/resource/CI gates
remain open in the [consumer capability plan](https://github.com/jake-w-liu/HybridPlasmaPIC.jl/blob/main/docs/PLASMA_COMPLETENESS_IMPLEMENTATION_PLAN.md#cap-011--electron-transport-and-complete-energy-accounting).

Native/GMP peak accounting and exclusive latency remain unaccepted. The Julia
allocation results do not establish either. Full exponent-range coupled physics,
normalized-dt RK storage, distributed/GPU support and simulator completeness are
outside this package slice.
