using RationalFunctionApproximation
using LinearAlgebra
import RationalFunctionApproximation: ContinuumApproximation


function (r::ContinuumApproximation)(A::AbstractArray, b::AbstractVector)::AbstractVector

    Z = r.fun.nodes
    W = r.fun.weights
    F = r.fun.values

    n = size(A, 1)

    numer = zeros(ComplexF64, n)
    denom = zeros(ComplexF64, n, n)

    for (z, w, f) in zip(Z, W, F)
        numer .+= w * f * ((A - z * I) \ b)
        denom .+= w * ((A - z * I) \ I)
    end

    return denom \ numer
end
