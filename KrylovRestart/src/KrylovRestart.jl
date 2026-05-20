module KrylovRestart

using FastGaussQuadrature
using LinearAlgebra
using SparseArrays
using KrylovKit
using Integrals

if Sys.ARCH == :aarch64
    using AppleAccelerate
end

include("./lanczos.jl")
include("./arnoldi.jl")
include("./rationalapprox.jl")
include("./krylov_approximation.jl")

export RationalApproximation, bestapprox_expm_data
export arnoldi, lanczos
export Params, StopCode, Results, krylov_approx
export krylov_approx_quad, krylov_approx_quad2, message

end # module KrylovRestart
