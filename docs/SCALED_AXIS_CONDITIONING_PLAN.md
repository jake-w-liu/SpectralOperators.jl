# Directional scaled-source axis conditioning

Version 0.1.8 candidate, 2026-11-01. Source and numerical design review are complete;
all Julia, formatter and runtime acceptance remain pending on GitHub. The
[registered 0.1.7 verification](SCALED_FOURIER_VERIFICATION.md) is historical
baseline evidence and does not certify this update.

## Source boundary

The cold normalized pipeline originally transformed the whole input before
applying its directional multiplier. A large background constant along the
consumed axis can therefore mix with a much smaller active source first.
The new small D2 tests place a tiny varying slice beneath a normal transverse
background and swap the consumed axis. They use independent four-point exact
matrices; execution will establish actual runtime behavior.

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
- [ ] Pass minimum/latest Julia package/Aqua and consumer checks on GitHub.
- [ ] Check formatting and complete compiler impact; measure public costs before
  making allocation or performance claims.

A certified null subtraction is not an exact-FFT theorem. Broad dynamic range
among genuinely active modes, arbitrary cancellation, non-admitted conditioning
cases and full hardware-mode pipelines remain separate boundaries. Fallback
retains the old method's limits; it does not claim to repair those cases.
