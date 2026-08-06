begin
    using Pkg
    Pkg.activate("KSMdev", shared = true)
end

begin
    using KrylovRestart
    using KrylovKit
    using LinearAlgebra
    using SparseArrays

    S = diagm(-1 => ones(39), 0 => fill(-2, 40), 1 => ones(39)) |> sparse
    I40 = I(40) |> sparse

    A = kron(I40, S) + kron(S, I40)

    b = rand(1600)

    exact, _ = exponentiate(A, 1, b)

    err_test(exact::AbstractVector, result::AbstractVector) = norm(exact - result) / norm(exact)
end
