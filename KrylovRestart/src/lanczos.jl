"""
`lanczos(A,q₁,m)`
Perform an Lanczos run for a Krylov Subspace ``\\cal{K}_m(A,q_1)`` of dimension `m` 
"""
function lanczos(A, b, m)

    iter = LanczosIterator(A, b)
    fac = initialize(iter) # Generates first basis vector

    for _ ∈ 1:m
        expand!(iter, fac)
    end

    Q = fac |> basis |> collect |> stack
    Q = Q[:, 1:m]
    q_next = basis(fac)[end]
    Tfull = fac |> rayleighquotient
    T = Tfull[1:m, 1:m]
    η = Tfull[m+1, m]

    return (Q, T, η, q_next)

end
