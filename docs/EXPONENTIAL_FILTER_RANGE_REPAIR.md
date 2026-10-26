# Exponential filter at infinite converted strength

## Pre-Edit Brief (2026-10-27)

- Mode: repair the canonical shared filter multiplier.
- Requirement: accepted `α=Inf` removes every non-DC mode for every accepted positive `p`; finite fields must not become NaN merely because the mode power underflows. Preserve ordinary finite-strength arithmetic and zero-strength behavior.
- Baseline: registered 0.1.3, tree `036881a8627ae86e285825891bb427dd66508c08`, on Julia 1.13.3. Both serial and distributed `(8,4,4)` Float32/Float64 fields `3+cos(2πx/8)` return 128 NaNs at `α=Inf,p=10000`; `p=8` returns constant 3. Raw reproduction is retained in `HybridPlasmaPIC-audit-evidence/2026-10-27/distributed-takeover/filter-repro-v2-1133.log`.
- Root cause: `ratio2^p` underflows to zero for a positive non-DC mode; `Inf*0` in the multiplier creates NaN. The DC branch already preserves the zero mode separately.
- Contract: after conversion to the grid precision, positive infinite strength has transfer zero at every positive mode ratio and one at DC. This is the limit of `exp(-α r^(2p))` as α tends to infinity for fixed positive r and p. Positive finite strengths keep their existing arithmetic, including rounding behavior. No finite-output promise for arbitrary nonfinite raw kernel operands is added.
- Canonical owner: `_apply_exp_filter_multiplier!` in `src/filters.jl`, shared by serial FourierGrid and HPIC's distributed workspace. Update dependency version and tests; HPIC only adds a consumer regression.
- Compatibility: preserve the explicitly accepted infinite-strength API. Rejecting infinity would remove existing behavior and is unnecessary. No local consumer formula, parameter restriction, diagnostic suppression or tolerance change.
- Implementation: after the exact DC branch, choose zero transfer for positive infinite converted α; otherwise evaluate the existing expression unchanged. Release as 0.1.4 after review and CI.
- Verification: baseline/repaired serial and actual MPI 1/2/4/8, minimum Julia 1.12.7 and latest stable 1.13.3; both precisions, large/maximum accepted powers, finite strengths converting to infinity, DC/singletons/finite controls, independent transfer checks and failure rollback. Full SO tests, Aqua, concrete/compiler report comparison, exact ordinary allocations, full HPIC operator regressions. Root coordinates publication and exact-SHA CI.
- Wrong patches rejected: preserve DC only (non-DC underflow regression); reject infinite α (accepted API regression); special-case only p=10000 (different large/max powers); return an all-zero field (nonzero mean and singleton checks); rewrite finite arithmetic (bitwise control checks).

## Progress

- [x] Reproduce on registered dependency and trace canonical owner.
- [x] Implement and verify serial/shared-kernel regression: 7,674 package assertions on Julia 1.12.7 and 1.13.3; 128 finite-strength profiles per runtime reproduce baseline hashes and three zero-allocation samples exactly.
- [ ] Verify distributed consumer across all ranks and runtimes.
- [ ] Review compiler/resources and publish dependency plus consumer regression.
