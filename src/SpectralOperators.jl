module SpectralOperators

using FFTW

export FourierGrid
export deriv!, deriv, gradient!, divergence!, curl!, laplacian!, project_divfree!
export exp_filter!, dealias_two_thirds!, fft_friendly_size, with_fftw_wisdom
export smoothing_transfer, BinomialSmoothWorkspace, binomial_smooth!
export SBP1D, sbp_deriv!, sbp_deriv, sbp_deriv_x!, FourierDerivYWorkspace, fourier_deriv_y!

# Supported qualified backend contract; see docs/FOURIER_BACKEND_CONTRACT.md.
public _fourier_metadata, _apply_ik_store!, _apply_ik_store_scaled!, _apply_ik_accum!
public _curl_component_store!, _apply_laplacian_multiplier!, _laplacian_gain
public _project_fourier_modes!, _apply_exp_filter_multiplier!, _apply_two_thirds_mask!
public _fft_safe_magnitude, _fft_input_scale, _copy_fft_input!, _store_fft_output!

include("periodic_fourier.jl")
include("differential_operators.jl")
include("divergence_control.jl")
include("filters.jl")
include("fdx_fourier_yz.jl")

end # module SpectralOperators
