"""
`arnoldi(A,q₁,m)`
Perform an Anoldi run for a Krylov Subspace ``\\cal{K}_m(A,q_1)`` of dimension `m`


Returns: \n
`Q`, `m` x `size(A,1)` orthonormal \n
`H`, `m` x `m` unreduced upper Hessenberg matrix \n
`η`, the final off diagonal term \n
`q_next`, the next basis vector that would be added to `Q`
"""
function arnoldi(A, b, m)

    iter = ArnoldiIterator(A, b)
    fac = initialize(iter)

    for _ in 1:m
        expand!(iter, fac)
    end

    Q = basis(fac) |> collect |> stack
    q_next = Q[:, m + 1]
    Q = Q[:, 1:m]
    Hfull = rayleighquotient(fac)
    H = Hfull[1:m, 1:m]
    η = Hfull[m + 1, m]

    return (Q, H, η, q_next)

end
