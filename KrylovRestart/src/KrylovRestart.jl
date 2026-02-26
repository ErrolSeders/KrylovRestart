module KrylovRestart

using LinearAlgebra
using SparseArrays
using KrylovKit

include("./lanczos.jl")
include("./arnoldi.jl")
include("./rationalapprox.jl")
include("./funm_krylov_restart.jl")

export RationalApproximation, bestapprox_expm_data
export arnoldi, lanczos
export Params, StopCode, Results, funm_krylov_restart_full, funm_krylov_restart_ra, message

end # module KrylovRestart
