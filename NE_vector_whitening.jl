using Pkg

begin
    Pkg.activate("KSMdev", shared = true)

    using LinearAlgebra
    using SparseArrays
    using MatrixDepot
    using KrylovRestart
end

retrieve_matrix() = begin
    mat = mdopen("HB/nos6")
    mat.A |> sparse
end

abs_err(exact, result) = norm(exact - result)
rel_err(exact, result) = abs_err(exact, result) / norm(exact)

function compute_exact(A, b)
    F = A |> Matrix |> Symmetric |> eigen
    λ = F.values
    @info minimum(real.(λ))
    @info maximum(real.(λ))
    V = F.vectors

    invsqrtλ = λ .|> sqrt .|> inv

    return V * (invsqrtλ .* (V' * b))
end

function test_alg(f::Function, args...; kwargs...)
    @info "Retreiving Matrix"
    A = retrieve_matrix()
    b = rand(size(A, 1))

    @info "Calculating Exact Solution"
    exact = compute_exact(A, b)

    tr = Trace(0, nothing, Dict{String, AbstractVector}())

    inv_sqrt = z -> z^(-1 / 2)

    @info "Calculating Krylov Approximation"
    res = f(inv_sqrt, A, b, args...; trace = tr, kwargs...)

    abserr = abs_err(exact, res)
    relerr = rel_err(exact, res)

    @info "Absolute Error: $abserr | Relative Error: $relerr"


    return tr
end
