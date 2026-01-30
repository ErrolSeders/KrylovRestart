"""
`arnoldi(A::AbstractArray,q₁::AbstractVector,m::Integer)`
Perform an Anoldi run for a Krylov Subspace ``\\cal{K}_m(A,q_1)`` of dimension `m` 
"""
function arnoldi(A, b, m)

    iter = ArnoldiIterator(A, b)
    fac = initialize(iter)

    for _ ∈ 1:m
        expand!(iter, fac)
    end

    Q = basis(fac) |> collect |> stack
    q_next = Q[:, m+1]
    Q = Q[:, 1:m]
    Hfull = rayleighquotient(fac)
    H = Hfull[1:m, 1:m]
    η = Hfull[m+1, m]

    return (Q, H, η, q_next)

end
