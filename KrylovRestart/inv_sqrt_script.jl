using Pkg
Pkg.activate("KSMdev", shared = true)

using LinearAlgebra
using SparseArrays
using MatrxDepot
using KrylovRestart
using JLD2

const N = 100

A = -(mdopen("poisson", N).A)

f = (z) -> z |> sqrt |> inv
