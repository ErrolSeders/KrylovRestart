module KrylovRestart

#using RationalFunctionApproximation
#using ComplexRegions
using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using SparseArrays
using KrylovKit
using Integrals

if Sys.ARCH == :aarch64
    using AppleAccelerate
end

include("./utils.jl")
include("./lanczos.jl")
include("./arnoldi.jl")
include("./rationalapprox.jl")
include("./quadrature_utils.jl")
include("./krylov_approximation.jl")

export RationalApproximation, bestapprox_expm_data
export arnoldi, lanczos
export StopCode, KrylovTrace, krylov_approx
export krylov_approx_quad, krylov_approx_quad2, message
export krylov_approx_chen_implicit, krylov_approx_chen_explicit

end # module KrylovRestart
