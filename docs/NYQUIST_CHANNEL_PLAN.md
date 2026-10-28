# Full-Laplacian Nyquist energy channel

This is the sealed operator design and local verification record. The channel
is registered, tagged and released in 0.1.6; see the
[release provenance](SCALED_FOURIER_VERIFICATION.md#release-provenance).
Its pre-edit base was `bdc5d3c919fe22df41dbb473ae107bc5a1d29e91` (0.1.5).
HPIC hyperresistive composition remains a separate acceptance gate.

`nyquist_channel!(out, f, g, d)` stores `N_d f = k_N P_N f`, where `P_N`
is the orthogonal projection onto the Nyquist mode along axis d. It preserves
all modes along the other axes. Odd axis lengths have no such mode and return
zero. Input/output may alias; malformed layout/axis and nonfinite input are
rejected before public output changes. The implementation reuses canonical
metadata, range-safe forward/inverse scaling, FFT plans and complex scratch.
It adds no persistent workspace and leaves existing operators unchanged.

The real first derivative D_d is skew-adjoint and omits axis Nyquist. The full
Laplacian retains it, hence

```
-L = sum_d(transpose(D_d) * D_d + k_N^2 P_N),
-<f,Lf> = sum_d(||D_d f||^2 + ||N_d f||^2).
```

This is a periodic inner-product identity. It does not assert a pointwise
product rule or a resolved sine derivative at Nyquist. The channel is the
positive discrete completion used by HPIC's hyperresistive electron heat.

Fixed pre-edit gates: independent direct Fourier/projector matrices; F32/F64
D1–3 odd/even/mixed/singleton grids; error ≤512eps(T) times operation scale;
in-place, shape/axis/finite-input errors; extreme finite amplitudes; zero
warmed allocations; minimum/latest Julia and the complete neighboring suite.

- [x] Source-supported pre-edit contract and factorization.
- [x] Canonical operator and 248 independent assertions per runtime.
- [x] Complete shared-operator suite and Aqua: 7986 assertions on Julia1.12.7
  and Julia1.13.3, with source at this candidate.
- [x] Dependency review, publication, registration and tagged release; see the release record.
