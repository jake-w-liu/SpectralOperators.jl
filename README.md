# SpectralOperators.jl

Reusable Julia operators for periodic Fourier grids and mixed SBP/Fourier
discretizations.

<!-- This package is the shared numerical-operator layer used by the sibling
`FourierPSTD.jl` and `HybridPlasmaPIC.jl` packages in this workspace. It has no dependency
on either solver package and does not use solver-global mutable state. -->

## Features

- Collocated periodic `FourierGrid` workspaces in 1D, 2D, and 3D.
- In-place first derivatives, gradient, divergence, curl, Laplacian, and
  divergence-free projection.
- Exponential spectral filtering and Orszag two-thirds dealiasing.
- FFTW wisdom import/export helper and FFT-friendly 7-smooth grid-size helper.
- Periodic binomial smoothing with the analytic transfer function exposed.
- A second-order SBP first derivative in `x` plus Fourier derivative in `y` for
  mixed non-periodic/periodic plasma-shock operators.

## Installation

For local development from this workspace:

```julia
using Pkg
Pkg.develop(path = "/path/to/SpectralOperators.jl")
```


## Quick Start

```julia
using SpectralOperators
using LinearAlgebra

n = (64,)
L = (2 * pi,)
g = FourierGrid(n, L)

x = [((i - 1) * g.dx[1]) for i in 1:n[1]]
f = sin.(3 .* x)

dfdx = deriv(f, g, 1)
exact = 3 .* cos.(3 .* x)

@assert norm(dfdx .- exact) / norm(exact) < 1e-12
```

For repeated use, prefer the in-place API:

```julia
out = similar(f)
deriv!(out, f, g, 1)
```

`FourierGrid` owns FFT plans and scratch arrays. Reuse it for repeated operator
calls on the same grid shape and element type. FFT-backed operators support
`Float32` and `Float64` fields.

## Vector Operators

```julia
g = FourierGrid((32, 32, 32), (2 * pi, 2 * pi, 2 * pi))
f = randn(Float64, g.n...)

grad = ntuple(_ -> similar(f), 3)
gradient!(grad, f, g)

A = ntuple(_ -> randn(Float64, g.n...), 3)
B = ntuple(_ -> similar(f), 3)
curl!(B, A, g)

divB = similar(f)
divergence!(divB, B, g)
```

`curl!` requires output component arrays that do not alias each other or the
input component arrays. This avoids silent overwrite of components that are
still needed for later curl components.

## Divergence-Free Projection

```julia
B = ntuple(_ -> randn(Float64, g.n...), 3)
project_divfree!(B, g)

divB = similar(B[1])
divergence!(divB, B, g)
```

The projection uses the same Nyquist-zeroed derivative wavenumbers as
`divergence!`: it keeps modes whose derivative wavenumber is zero unchanged and
removes the longitudinal component from derivative-resolved modes. Component
arrays must be distinct.

## Filters and Smoothing

```julia
field = randn(Float64, 64, 64)
g2 = FourierGrid((64, 64), (2 * pi, 2 * pi))

exp_filter!(field, g2; α = 36, p = 8)
dealias_two_thirds!(field, g2)
binomial_smooth!(field, g2; passes = 1)

smooth_work = BinomialSmoothWorkspace(g2)
binomial_smooth!(field, g2, smooth_work; passes = 1)

tf = smoothing_transfer(3.0, g2.dx[1]; passes = 2)
```

`exp_filter!`, `dealias_two_thirds!`, and `binomial_smooth!` preserve the mean
mode by construction. The workspace-backed smoothing overload avoids allocating
the line buffer in repeated calls.

## Mixed SBP/Fourier Operators

```julia
s = SBP1D(65, 1.0)
x = range(0.0, 1.0; length = s.n)
f = collect(sin.(pi .* x))

df = similar(f)
sbp_deriv!(df, f, s)

u = repeat(f, 1, 16)
dux = similar(u)
sbp_deriv_x!(dux, u, s)

duy = similar(u)
fourier_deriv_y!(duy, u, 2 * pi)

ywork = FourierDerivYWorkspace(u, 2 * pi)
fourier_deriv_y!(duy, u, ywork)
```

The SBP derivative is not an in-place transform: the output must not alias the
input. Reuse `FourierDerivYWorkspace` for repeated transverse Fourier
derivatives on the same matrix size.

## FFTW Wisdom

```julia
with_fftw_wisdom("fftw_wisdom.dat") do
    g = FourierGrid((128, 128), (2 * pi, 2 * pi))
    nothing
end
```

The helper imports wisdom if the file exists and attempts to export wisdom after
the block, even if the block throws. Import/export failures are non-fatal.

## Testing

From the package directory:

```julia
using Pkg
Pkg.test()
```

The test suite covers closed-form derivative oracles, spectral identities,
divergence-free projection, Nyquist handling, filters, dealiasing, FFTW wisdom,
binomial smoothing, mixed SBP/Fourier behavior, defensive shape and alias checks,
and steady-state zero allocation for the main periodic Fourier operators.

## Registration Readiness

The package has its own git repository and license file. Before Julia General
registration, run the normal release checks from a clean checkout and tag the
release intended for registration.
