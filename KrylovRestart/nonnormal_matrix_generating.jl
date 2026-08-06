begin
    using Pkg
    Pkg.activate("KSMdev", shared = true)
end

begin
    using KrylovRestart
    using LinearAlgebra
    using SparseArrays
    using JLD2
end

function generate_normal()


end


"""
Norm of the commuter. Scales quadratically as matrices become nonnormal.
"""
commuter_metric(A) = norm(A * A' - A' * A)


"""
Henrici's departure from normality. A scalar index of the nonnormality of a matrix.
Given by the Frobenius norm of the strictly upper triangular component of the schur decomposition.

For a normal matrix = 0
"""
function henricis_departure(A)
    T, Z, _ = schur(A)
    R = triu(T, 1)
    return norm(R)
end
