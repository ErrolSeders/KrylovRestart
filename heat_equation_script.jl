"""
This script provides a sample symmetric problem for testing algorithms from `KrylovRestart`

Author: Errol Bucy
Date: 17.06.26
"""

using Pkg
Pkg.activate("KSMdev", shared = true)

using LinearAlgebra
using SparseArrays
using Tullio
using DataFrames
using KrylovKit
using UnicodePlots

push!(LOAD_PATH, joinpath(@__DIR__, "KrylovRestart"))
using KrylovRestart

"""
∂ₜu - Δu = 0 in Ω = (0,1)², t > 0
u(x,t) = 0 on Γ = ∂Ω,
t > 0 u(x,0) = u₀(x) in Ω

Heat equation discretized by the five-point stencil on a uniform grid involving
`n_1` interior grid points in each direction
"""

(⊗)(A, B) = kron(A, B)

n1 = 200

build_Lh(n1) = begin
    T = SymTridiagonal(fill(-2.0, n1), fill(1.0, n1 - 1))
    Id = I(n1) |> sparse
    h = 1 / (n1 + 1)
    Lh = (1 / h^2) * (T ⊗ Id + Id ⊗ T)
end

function eigenfuncs_initial_vec(n1::Int)
    h = 1 / (n1 + 1)

    I = collect(1:n1)
    Ip = collect(1:n1)

    # Precompute the sine matrices
    Sx = @. sin(I .* Ip' * π * h)   # n1 × n1
    Sy = @. sin(I .* Ip' * π * h)

    W = [1.0 / (ip + jp) for ip in Ip, jp in Ip]

    @tullio u[i, j] :=
        Sx[i, ip] * Sy[j, jp] * W[ip, jp]

    return reshape(u, :, 1)
end

u0 = eigenfuncs_initial_vec(n1)[:]

begin
    t = 0.2
    Lh = build_Lh(n1)
    tLh = t * Lh
    exact, _ = exponentiate(Lh, t, u0, issymmetric = true, maxiter = 2000, tol = 1.0e-16)
end

function test_alg()

end
