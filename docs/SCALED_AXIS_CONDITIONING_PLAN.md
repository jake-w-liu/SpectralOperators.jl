# Directional scaled-source axis conditioning

Version 0.1.8 package verification, 2026-11-01. All five jobs pass in the
[exact-revision CI run](https://github.com/jake-w-liu/SpectralOperators.jl/actions/runs/41276540274)
at commit `1c3d0cc67278364cfe899bc20a4079bbd8f5a57e`, tree
`f768bbfa33c6b9e4eec6f03d3098b35be5c86928`. Each Linux/macOS package and Aqua job
records 26,986 passing assertions on actual Julia 1.12.7 or 1.13.3. The fifth
job passes JuliaFormatter 2.14.0 on Julia 1.12.7. The
[General registration](https://github.com/JuliaRegistries/General/pull/174315)
merged on 2026-11-01 UTC. The annotated tag and
[public release](https://github.com/jake-w-liu/SpectralOperators.jl/releases/tag/v0.1.8)
point to that exact tested commit and registered tree. The
[release record](SCALED_FOURIER_VERIFICATION.md) separates historical checks from
this update's acceptance and remaining limits.

## Source boundary

The cold normalized pipeline originally transformed the whole input before
applying its directional multiplier. A large background constant along the
consumed axis can therefore mix with a much smaller active source first.
The small D2 tests place a tiny varying slice beneath a normal transverse
background and swap the consumed axis. Their independent four-point exact
matrix comparisons pass in all four package jobs above.

## Bounded canonical change

For each line in axis `d`, let `r` be its first stored value. The field formed
by extending `r` along that line has Fourier index `m_d=0`. The centered first
derivative, axis Nyquist projector and centered derivative of the full
Laplacian all annihilate it. The full Laplacian alone can retain transverse
variation, so its ordinary and cold preparation remain unchanged.

Directional preparation may replace the input by `f-r` only when it can retain
the exact stored difference. Axis-constant pairs become zero before scaling;
their unrelated magnitude cannot force loss of an active tiny line. A common
binary downshift selected from active operands and existing unit-gain headroom
prevents subtraction overflow. Binary round trips certify every scaled operand.
The error-free TwoSum residual then certifies exact subtraction in nearest,
gradual Float32/64 arithmetic. Existing upward normalization retains its binary
exponent. If any certificate fails, the original preparation overwrites scratch
and keeps its existing scaling and rejection behavior.

TwoSum's exact decomposition, including gradual underflow, is established by
Algorithm 3.1 and Theorem 3.4 of
[Ogita, Rump and Oishi (2005)](https://ogilab.w.waseda.jp/ogita/math/doc/2005_OgRuOi.pdf).
The implementation uses ordinary ordered additions/subtractions and does not
introduce reassociation, higher precision or a solver dependency.

The same borrowed `g.cbuf`, scale/exponent interface, readonly metadata contract,
input disjointness checks and one-caller lifetime remain. There are no additional
FFT plans, persistent arrays or public options. Centered-null short axes and
odd Nyquist axes retain their existing early-return validation.

## Progress and limits

- [x] Inspect all four scaled kinds, numerical preparation, metadata and consumers.
- [x] Derive operator-null and exact scalar admission conditions before editing.
- [x] Implement one shared directional preparation and retain Laplacian behavior.
- [x] Include essential axis swaps, constants, exact-factor/range and alias controls.
- [x] Add a GitHub formatter check with a downloadable patch on failure.
- [x] Select JuliaFormatter 2.14.0, whose [upstream regression](https://github.com/JuliaEditorSupport/JuliaFormatter.jl/blob/v2.14.0/test/default_style.jl#L4984) covers Julia's
  `public` declarations. The initial v1 artifact joined these keywords to names;
  that patch was retained as evidence and was not applied.
- [x] Review and apply the v2 runner's formatting artifact, retaining the module's
  `public` declarations byte-for-byte.
- [x] Pass all four minimum/latest Linux/macOS package/Aqua jobs at the tested
  revision, plus its GitHub formatting job.
- [x] Verify General publication and the exact tested tag/release tree.
- [ ] Verify downstream consumer composition on GitHub.
- [ ] Check complete compiler impact and measure public costs before making new
  allocation or performance claims.

A certified null subtraction is not an exact-FFT theorem. Broad dynamic range
among genuinely active modes, arbitrary cancellation, non-admitted conditioning
cases and full hardware-mode pipelines remain separate boundaries. Fallback
retains the old method's limits; it does not claim to repair those cases.
