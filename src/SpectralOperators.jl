module SpectralOperators

using FFTW

export FourierGrid
export deriv!, deriv, gradient!, divergence!, curl!, laplacian!, project_divfree!
export exp_filter!, dealias_two_thirds!, fft_friendly_size, with_fftw_wisdom
export smoothing_transfer, binomial_smooth!
export SBP1D, sbp_deriv!, sbp_deriv, sbp_deriv_x!, fourier_deriv_y!

include("periodic_fourier.jl")
include("differential_operators.jl")
include("divergence_control.jl")
include("filters.jl")
include("fdx_fourier_yz.jl")

end # module SpectralOperators
