"""
    krylov_approx(f, A, b; m=30, max_restarts=200, tol=nothing, bound=true, exact=nothing,
                 min_decay=0.95, trace=nothing, callback=nothing)

Restarted Krylov approximation of `f(A)b` (Alg. 1) or `r(A)b` for a `RationalApproximation` (Alg. 2).

The default API returns only the approximation vector and does not store per-restart history.
To collect convergence information for experiments, pass `trace=KrylovTrace()` and/or a
`callback(info)` function (e.g. to save iterates, plot errors, etc.).
"""

function update_alphas(α1, α2, H)

    if size(H) == (1, 1)
        μ = λ = H[1]
    else
        s, _ = eigen(H)
        μ = s .|> real |> minimum
        λ = s .|> real |> maximum
    end
    return min(μ, α1), max(λ, α2)
end

"""
Build the expanded Arnoldi decomposition matrix `Hhat`
The resulting matrix is block lower triangular and of size `km + m + l` x `km + m + l`
Note: `concat_block` is expected to be square
"""
function _build_Hhat(Hhat, H::AbstractMatrix, η_prev, k, m; η = nothing, concat_block = [])

    km = k * m
    l = size(concat_block, 1)
    Hhat_expand = zeros(eltype(Hhat), km + m + l, km + m + l)

    @views Hhat_expand[1:km, 1:km] .= Hhat[1:km, 1:km]
    @views Hhat_expand[(km + 1):(km + m), (km + 1):(km + m)] .= H

    Hhat_expand[((k - 1) * m) + 1, (k - 1) * m] = η_prev
    if !isnothing(η)
        Hhat_expand[km + m + 1, km + m] = η
    end
    Hhat_expand[(km + m + 1):(km + m + l), (km + m + 1):(km + m + l)] = concat_block

    return Hhat_expand
end

"""
Build the expanded Lanczos decomposition matrix `That`
The resulting matrix is `Tridiagonal` and of size `km + m + l` x `km + m + l` where `l`
is the size of `concat_block`.
Note: `concat_block` is expected to be square and preserve the Tridiagonal structure!
"""
function _build_Hhat(That, T::SymTridiagonal, η_prev, k = nothing, m = nothing; η = 0.0, concat_block = [])

    That_sup_diag, That_diag, That_sub_diag = isempty(concat_block) ? (
            [diag(That, 1); [0]; T.ev],
            [diag(That); T.dv],
            [diag(That, -1); [η_prev]; T.ev],
        ) : (
            [diag(That, 1); [0]; T.ev; [0]; diag(concat_block, 1)],
            [diag(That); T.dv ; diag(concat_block)],
            [diag(That, -1); [η_prev]; T.ev; [η]; diag(concat_block, -1)],
        )

    return Tridiagonal(That_sub_diag, That_diag, That_sup_diag)
end


krylov_approx(f, A::AbstractArray, b::AbstractVector, m::Integer; kwargs...) =
    krylov_approx(f, A, b; m = Int(m), kwargs...)

"""
    krylov_approx(f, A, b; m=30, max_restarts=200, tol=nothing, bound=true, exact=nothing,
                 min_decay=0.95, trace=nothing, callback=nothing)

Compute `f(A)b` in the manner corresponding to Alg. 1 in the paper.

If `trace::KrylovTrace` is provided, it is filled with per-restart statistics.
If `callback` is provided, it is called once per restart as `callback(info)` where
`info` is a named tuple containing the measured stopping quantity (and `stop` if triggered).
"""
function krylov_approx(
        f::Function,
        A::AbstractArray,
        b::AbstractVector;
        m::Int = 30,
        max_restarts::Int = 200,
        tol = nothing,
        bound::Bool = true,
        exact::Union{Nothing, AbstractVector} = nothing,
        min_decay = 0.95,
        trace::Union{Nothing, KrylovTrace} = nothing,
        callback = nothing
    )

    if !isnothing(exact) && bound
        @warn "Both `exact` provided and `bound=true`; stopping will use bounds."
    end
    m < 1 && throw(ArgumentError("m must be ≥ 1"))

    Tdefault = float(real(eltype(A)))
    tolT = tol === nothing ? eps(Tdefault) : float(tol)
    TT = typeof(tolT)
    tolT = TT(tolT)
    min_decayT = TT(min_decay)

    # state for the linear-convergence-based stopping rules
    up_decay = DecayTracker(TT)
    err_decay = DecayTracker(TT)

    α1, α2 = if bound
        prevfloat(Tdefault(Inf)), nextfloat(Tdefault(-Inf))
    else
        zero(Tdefault), zero(Tdefault)
    end

    fk = zeros(promote_type(eltype(A), eltype(b)), length(b))

    β = norm(b)
    qm = (1 / β) * b

    η_prev = NaN

    is_A_hermitian = ishermitian(A)

    if !is_A_hermitian
        Hhat = Array{eltype(H)}(undef, 0, 0)
    end

    for k in 1:max_restarts

        set_restart!(trace, k)

        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        if bound
            α1, α2 = update_alphas(α1, α2, H)
        end

        if k == 1
            # If I constantly have to do different things for Hermitian vs. Non-Hermitian then maybe I should dispatch a different algorithm for Hermitian.
            is_A_hermitian ? (
                    begin
                        Hhat = Tridiagonal(
                            [H.ev;[0.0, 0.0]],
                            [H.dv;[α1, α2]],
                            [H.ev;[η, 1.0]]
                        )
                    end
                ) : (
                    begin
                        Hhat = zeros(eltype(H), m + 2, m + 2)
                        Hhat[1:m, 1:m] = H
                        Hbar[m + 1, m] = η
                        Hhat[(m + 1):(m + 2), (m + 1):(m + 2)] = [α1 0.0; 1.0 α2]
                    end
                )
        else
            Hhat_expand = _build_Hhat(Hhat, H, η_prev, k, m, η = η, concat_block = [α1 0.0; 1.0 α2])
            Hhat = Hhat_expand
        end
        η_prev = η

        @views h = f(Hhat)[((k - 1) * m + 1):((k - 1) * m + m), 1]
        fk .+= β * (Q * h)

        stop = stop_conditions!(
            A, α1, β, h, qm, fk, m,
            bound, exact, tolT, min_decayT,
            trace, callback, k,
            up_decay, err_decay
        )

        if stop !== nothing
            set_stop!(trace, stop)
            return fk
        end
    end

    set_stop!(trace, MaxRestarts)
    return fk
end

krylov_approx(r::RationalApproximation, A::AbstractArray, b::AbstractVector, m::Integer; kwargs...) =
    krylov_approx(r, A, b; m = Int(m), kwargs...)

"""
    krylov_approx(r::RationalApproximation, A, b; m=30, max_restarts=200, tol=nothing, bound=true,
                 exact=nothing, min_decay=0.95, trace=nothing, callback=nothing)

Compute `r(A)b ≈ f(A)b` in the manner corresponding to Alg. 2 in the paper.
"""
function krylov_approx(
        r::RationalApproximation,
        A::AbstractArray,
        b::AbstractVector;
        m::Int = 30,
        max_restarts::Int = 200,
        tol = nothing,
        bound::Bool = true,
        exact::Union{Nothing, AbstractVector} = nothing,
        min_decay = 0.95,
        trace::Union{Nothing, KrylovTrace} = nothing,
        callback = nothing
    )

    if !isnothing(exact) && bound
        @warn "Both `exact` provided and `bound=true`; stopping will use bounds."
    end
    m < 1 && throw(ArgumentError("m must be ≥ 1"))

    Tdefault = float(real(eltype(A)))
    tolT = tol === nothing ? eps(Tdefault) : float(tol)
    TT = typeof(tolT)
    tolT = TT(tolT)
    min_decayT = TT(min_decay)

    up_decay = DecayTracker(TT)
    err_decay = DecayTracker(TT)

    α1, α2 = if bound
        # Set to 'virtually' Inf so that LU does not fail for solves with '\'
        prevfloat(Tdefault(Inf)), nextfloat(Tdefault(-Inf))
    else
        zero(Tdefault), zero(Tdefault)
    end

    real_res = (eltype(A) <: AbstractFloat) && (eltype(b) <: AbstractFloat)

    fk = real_res ? (r.absterm * b) : complex.(r.absterm * b)

    β = norm(b)
    qm = (1 / β) * b

    nr_single = length(r.single_poles)
    poles = [r.single_poles; r.conj_poles]
    coeff = [r.single_coeff; r.conj_coeff]

    Bbar = zeros(ComplexF64, m + 2, length(poles))
    Bbar[end - 2, :] .= 1.0 + 0.0im
    s = one(eltype(Bbar))

    e1 = zeros(eltype(A), m + 2)
    e1[1] = one(eltype(A))

    is_A_hermitian = ishermitian(A)

    Hbar = is_A_hermitian ? nothing : zeros(eltype(A), m + 2, m + 2)

    for k in 1:max_restarts

        set_restart!(trace, k)

        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        if bound
            α1, α2 = update_alphas(α1, α2, H)
        end

        # Hbar is overwritten in this block
        is_A_hermitian ? (
                begin
                    Hbar = Tridiagonal(
                        [H.ev;[0.0, 0.0]],
                        [H.dv;[α1, α2]],
                        [H.ev;[η, 1.0]]
                    )
                end
            ) : (
                begin

                    Hbar[1:m, 1:m] .= H
                    Hbar[(m + 1):(m + 2), (m + 1):(m + 2)] = [α1 0.0; 1.0 α2]
                    Hbar[m + 1, m] = η
                end
            )

        qbar = zeros(ComplexF64, m + 2)
        for p in eachindex(poles)
            qbar[1] = s * Bbar[m, p]
            Bbar[1:(m + 2), p] = (poles[p] * I - Hbar) \ qbar
        end

        # contribution from single poles
        @views h = nr_single > 0 ? Bbar[1:(m + 2), 1:nr_single] * coeff[1:nr_single] : zeros(ComplexF64, m + 2)

        if nr_single < length(poles) # add in the contrabution from the conjugate poles
            @views h .+= 2 * real(Bbar[1:(m + 2), (nr_single + 1):end] * coeff[(nr_single + 1):end])
        end

        s = η

        if real_res
            h = h .|> real
        end

        @views fk .+= β * (Q * h[1:m])

        stop = stop_conditions!(
            A, α1, β, h, qm, fk, m,
            bound, exact, tolT, min_decayT,
            trace, callback, k,
            up_decay, err_decay
        )

        if stop !== nothing
            set_stop!(trace, stop)
            return fk
        end
    end

    set_stop!(trace, MaxRestarts)
    return fk
end

"""
    Compute `f(A)b` by continually approximating `f` using the AAA algorithm.
"""
# function krylov_approx_barycentric(
#         f::Function,
#         A::AbstractArray,
#         b::AbstractVector;
#         m::Int = 30,
#         max_restarts::Int = 200,
#         tol = nothing,
#         exact::Union{Nothing, AbstractVector} = nothing,
#         min_decay = 0.95,
#         trace::Union{Nothing, KrylovTrace} = nothing,
#         callback = nothing
#     )

#     m < 1 && throw(ArgumentError("m must be ≥ 1"))

#     Tdefault = float(real(eltype(A)))
#     tolT = tol === nothing ? eps(Tdefault) : float(tol)
#     TT = typeof(tolT)
#     tolT = TT(tolT)
#     min_decayT = TT(min_decay)

#     up_decay = DecayTracker(TT)
#     err_decay = DecayTracker(TT)

#     real_res = (eltype(A) <: AbstractFloat) && (eltype(b) <: AbstractFloat)

#     β = norm(b)
#     qm = (1 / β) * b

#     e1 = zeros(ComplexF64, m + 2)
#     e1[1] = 1.0 + 0.0im

#     # H is allocated once
#     Hbar = zeros(eltype(A), m + 2, m + 2)

#     is_A_hermitian = ishermitian(A)

#     α1, α2 = Inf, -Inf

#     Tvec = promote_type(eltype(A), eltype(b))
#     fk = zeros(real_res ? Tvec : ComplexF64, size(A, 1))

#     for k in 1:max_restarts

#         set_restart!(trace, k)

#         (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

#         vals = eigvals(H)

#         α1 = min(α1, vals .|> real |> minimum)
#         α2 = max(α2, vals .|> real |> maximum)

#         # Include the two extra spectral points that appear in Hbar.
#         vals_aug = vcat(vals, α1, α2)
#         R, c = fetch_contour_shape(vals_aug, 1.1)

#         contour = Circle(c, R)

#         r = approximate(f, contour)

#         Z = r.fun.nodes
#         W = r.fun.weights
#         WF = r.fun.w_times_f

#         # H is overwritten in this block
#         Hbar[1:m, 1:m] .= H
#         Hbar[(m + 1):(m + 2), (m + 1):(m + 2)] = [α1 0.0; 1.0 α2]
#         Hbar[m + 1, m] = η

#         # Evaluate the AAA barycentric rational on the projected matrix:
#         #   r(M) = (∑_j w_j f_j (z_j I - M)^{-1}) (∑_j w_j (z_j I - M)^{-1})^{-1}
#         # applied to the restart RHS q = e1.
#         q = e1
#         numer = zeros(ComplexF64, m + 2)
#         denom = zeros(ComplexF64, m + 2, m + 2)
#         for j in eachindex(Z)
#             M = Z[j] * I - Hbar
#             Fm = lu(M)
#             numer .+= WF[j] .* (Fm \ q)
#             denom .+= W[j] .* (Fm \ I)
#         end
#         h = denom \ numer

#         if real_res
#             h = real.(h)
#         end

#         @views fk .+= β * (Q * h[1:m])

#         stop = stop_conditions!(
#             A, α1, β, h, qm, fk, m,
#             true, exact, tolT, min_decayT,
#             trace, callback, k,
#             up_decay, err_decay
#         )

#         if stop !== nothing
#             set_stop!(trace, stop)
#             return fk
#         end
#     end

#     set_stop!(trace, MaxRestarts)
#     return fk
# end

"""
    integral_error_correction(f,H,Hs,ηs,order,contour_safety)

    evaluate the integral formulation for the error term at the transformed quadrature points

"""
function integral_error_correction(
        f::Function,
        H::AbstractArray,
        Hs,
        ηs::AbstractVector,
        order::Int,
        R::Float64,
        c,
        e1
    )

    m = size(H, 1)
    #@info "solving quad of order $order"

    x, w = contour_quad_nodes(order, R, c)

    S = zeros(ComplexF64, m)


    for (t, ω) in zip(x, w)
        y = (t * I - H) \ e1

        ϕ = prod(
            ηs[j] * ((t * I - Hs[j]) \ e1)[end] for j in eachindex(ηs)
        )

        S += ω * f(t) * ϕ * y
    end
    return S
end

"""
    krylov_approx_quad(f,A,b,)

Approximate ``f(A)b`` by Krylov subsapce method with a quadrature formulation of the error.
Here we assume that the function `f` is holomorphic in a neighborhood around the spectrum of `A`. In this case
we are guarenteed that the interal exists and the quadrature based error formula actually represents the error.

"""
function krylov_approx_quad(
        f::Function,
        A::AbstractArray,
        b::AbstractVector,
        m::Int
        ;
        tol = 1.0e-16,
        max_restarts = 200,
        contour_safety = 1.1
    )

    β = norm(b)

    qm = (1 / β) * b
    is_A_hermitian = ishermitian(A)

    (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

    Hs = Vector([H])
    ηs = [η]

    h = f(H)[1, :]
    fk = β * Q * h

    e1 = zeros(eltype(H), m)
    e1[1] = one(eltype(H))

    ```
       For the integral we throw a circular contour around the ritz values
    ```
    R, c = fetch_contour_shape(H, contour_safety)

    quad1::Int64 = 8
    quad2::Int64 = round(sqrt(2) * quad1)

    for _k in 2:max_restarts
        #@info "Iteration $_k"
        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        accurate = false
        refined = false

        h2 = zeros(m)
        err_prev = 1.0e16
        while !accurate

            R, c = fetch_contour_shape(H, contour_safety)
            h1 = integral_error_correction(f, H, Hs, ηs, quad1, R, c, e1)
            h2 = integral_error_correction(f, H, Hs, ηs, quad2, R, c, e1)

            err = norm(h2 - h1)

            #@info "Quadrature error $(err)"
            if err < tol
                #@info "accurate!"
                accurate = true
            else
                #@info "refining"
                quad1 = quad2
                quad2 = round(sqrt(2) * quad1)
                refined = true
            end
        end

        update = β * Q * h2

        #@info "norm of update $(norm(update))"

        fk += update
        if norm(update) < tol
            return fk
        end

        push!(Hs, H)
        push!(ηs, η)

        if !refined
            quad2 = quad1
            quad1 = round(quad2 / sqrt(2))
        end
    end
    return
end

function make_integral_function(f, H, Hs, ηs, R, c, e1)

    m = size(H, 1)

    function F(x, p = nothing)

        # Coordinate transform from [-1,1] to circle contour
        θ = π * x
        expθ = exp(im * θ)
        t = c + R * expθ

        y = (t * I - H) \ e1

        ϕ = one(eltype(H))

        res = Vector{ComplexF64}(undef, m)

        for j in eachindex(ηs)
            ldiv!(res, (t * I - Hs[j]), e1)
            ϕ *= ηs[j] * res[end]
        end

        return (0.5 * R * expθ * f(t) * ϕ) * y
    end

    return F
end

"""
    krylov_approx_quad2(f,A,b,m)

    Approximate ``f(A)b`` using Gauss-Kronrod quadrature
    directly assuming that we can attain an accurate quadrature based
    error correction in 1 shot using Julia's quadrature libraries

    !! Right now this appears to give incorrect answers for non-symmetric matrices !!
"""
function krylov_approx_quad2(
        f::Function,
        A::AbstractArray,
        b::AbstractVector,
        m::Int
        ;
        tol = 1.0e-16,
        max_restarts = 200,
        contour_safety = 1.1,
        alg = QuadGKJL()
    )

    real_res = false

    if (eltype(A) <: AbstractFloat) && (eltype(b) <: AbstractFloat)
        real_res = true
    end

    β = norm(b)

    qm = (1 / β) * b
    is_A_hermitian = ishermitian(A)

    (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

    Hs = [UpperHessenberg(H)]
    ηs::Vector{eltype(H)} = [η]

    h = f(H)[1, :]
    fk = β * Q * h

    e1::Vector{eltype(H)} = zeros(eltype(H), m)
    e1[1] = one(eltype(H))

    for _k in 2:max_restarts
        @info "Iteration $_k"
        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        H = UpperHessenberg(H)

        R, c = fetch_contour_shape(H, contour_safety)
        F = make_integral_function(f, H, Hs, ηs, R, c, e1)
        prob = IntegralProblem(F, (-1.0, 1.0))
        sol = solve(prob, alg, reltol = 1.0e-16, abstol = 1.0e-16)

        update = β * Q * sol.u

        @info "norm of update $(norm(update))"

        fk += update
        if norm(update) < tol
            if real_res

                fk = fk .|> real
            end
            return fk
        end

        push!(Hs, H)
        push!(ηs, η)
    end
    if real_res
        fk = fk .|> real
    end
    return fk
end

"""
Implementation of ``f(A)b`` for HPD matrices with an error indicator given
from "2-Norm Error bounds and estimates for Lanczos approximations to
linear systems and rational matrix functions".

For this algorithm we need `f` to be a rational function
(or able to be accurately approximated by a rational function)
and furthermore `A` to be HPD.

"""
function krylov_approx_2norm(f, A, b, m)

end

"""
Chen dissertation: Lemma 7.8
"""
function norm_of_hwz(z, w::Real, a, b)
    z_re = real(z)
    z_im = imag(z)
    left = abs((a - w) / (a - z))
    right = abs((b - w) / (b - z))
    x = (z_re^2 + z_im^2 - z_re * w) / (z_re - w)
    if x >= a && x <= b && z_im != 0
        middle = abs((1.0 / z_im) * (z - w))
    else
        middle = 0
    end

    @assert middle != NaN

    return max(left, middle, right)
end

function build_contour_chen(f, T, contour_safety, order, w)
    λs, _ = eigen(T)

    λmax = maximum(λs)
    λmin = minimum(λs)

    R, c = fetch_contour_shape(λs, contour_safety)
    nodes, weights = contour_quad_nodes(order, R, c)

    fval = nodes .|> f
    absfval = fval .|> abs

    norm_curry = z -> norm_of_hwz(z, w, λmin, λmax)
    norm_val = nodes .|> norm_curry

    integral_constant = (R * im) / (2 * π)

    return R, c, nodes, weights, fval, absfval, norm_val, integral_constant
end

"""
Implementation of an implicit restarted KSM using an error indicator adapted from Chen et al. and an implicit restarting scheme derived from the cauchy integral form of the KSM approximation

For this algorithm we require that `A` be Hermitian as the error bound is only valid for the Lanczos decomposition.
The adaptation for the restarted scheme was derived assuming that we have a decomposition as given by Eiermann and Ernst (2006).

This algorithm is rather sensitive to its hyperparameters, we ideally want to choose `w` less than the minimum eigenvalue of `A`. Furthermore, in this algorithm we try determine a countour the fully encloses the eigenvalues of `A` and all subsequent `Tk` in one fell swoop based upon the initial ritz values extracted from `T1`. The ritz values of `T1` should be an okay approximation to the extreme values of `A`. Therefore, we have the variable `contour_safety` which will scale the radius of the circular contour beyond the initial maximum ritz value. A sufficiently large number of quadrature nodes must be choosen as well with `order`.

Two key variables that need to be controlled for this algorithm to work are the value `w` which is a single complex point that we need to choose
to be disjoint from both the spectrum of `A` and the ritz-values of each `T`.
We also will need to choose a contour for evaluating an integral such that every value of
the contour is far away from the spectrum of `A` and ritz-values of each `T`.
"""
function krylov_approx_chen_implicit(
        f, A, b, m, w;
        tol = 1.0e-16,
        max_restarts = 200,
        contour_safety = 2.0,
        order = 80
    )

    @assert ishermitian(A)
    n = size(A, 1)

    err = 0.0

    β = norm(b)
    qm = (1 / β) * b

    (Q, T, η, qm) = lanczos(A, qm, m)

    T = SymTridiagonal(T)

    #The contour is set once based on the initial Ritz values
    R, c, nodes, weights, fval, absfval, norm_val, error_integral_constant = build_contour_chen(f, T, contour_safety, order, w)
    update_integral_constant = (R / (2 * π * im))

    h = f(T)[1, :]
    fk = β * Q * h

    # Dets are computed in log-space for numerical stability
    w_det_accum = (T - w * I) |> logabsdet |> first
    dets_accum = [(T - node * I) |> logabsdet |> first for node in nodes]

    μs = Vector{ComplexF64}[]
    μ_coeff = iseven(m) ? 1.0 : -1.0 #(-1)^m

    Q_herm_b = Q' * b
    Q_herm_bs = [Q_herm_b]

    ηs = [η]

    e1 = zeros(m); e1[1] = 1.0
    em = zeros(m); em[end] = 1.0

    # eₘᵀ(T - zI)⁻¹ at each quad node
    right_vecs = [mapreduce(node -> (T - node * I)' \ em, hcat, nodes)]

    #Preallocation - Continually write over the same arrays
    left_vec = Array{ComplexF64}(undef, m, order)
    right_vec = Array{ComplexF64}(undef, m, order)
    cauchy_integral = Array{ComplexF64}(undef, n, order)
    dets = Array{Tuple{Float64, ComplexF64}}(undef, order)
    recurrence_terms = Array{ComplexF64}(undef, size(A, 1), order)

    for k in 2:max_restarts

        (Q, T, η, qm) = lanczos(A, qm, m)

        Q_herm_b = Q' * b
        push!(Q_herm_bs, Q_herm_b)

        T::SymTridiagonal = SymTridiagonal(T)

        sub_diag = diag(T, -1)
        sub_diag_prod = prod(sub_diag)

        w_det_accum += (T - w * I) |> logabsdet |> first

        #Compute the quantities in the integral at each quadrature node
        for (j, node) in enumerate(nodes)

            T_shift = T - node * I
            F = factorize(T_shift)

            left_vec[:, j] = F \ e1
            right_vec[:, j] = F' \ em
            cauchy_solve = F \ Q_herm_b
            cauchy_integral[:, j] = weights[j] * Q * cauchy_solve

            dets[j] = logabsdet(F)
            dets_accum[j] += dets[j] |> first

        end

        push!(right_vecs, deepcopy(right_vec))

        det_sign = [det[1] for det in dets]
        inv_logabsdet_val = [-det[2] for det in dets]
        inv_det_val = (μ_coeff * det_sign) .* exp.(inv_logabsdet_val)

        push!(μs, sub_diag_prod .* inv_det_val)

        fill!(recurrence_terms, 0.0 + 0.0im)

        η_suffix = one(η)
        μ_suffix = ones(eltype(μs[1]), order)
        sgn = one(η)

        for i in (k - 1):-1:1
            coeff = sgn * η_suffix * η

            @inbounds for j in 1:order
                recurrence_terms[:, j] += weights[j] * fval[j] * coeff * Q * μ_suffix[j] * left_vec[:, j] * dot(right_vecs[i][:, j], Q_herm_bs[i])
            end

            η_suffix *= ηs[i]
            μ_suffix .*= μs[i]
            sgn = -sgn
        end

        h = (update_integral_constant * (sum(cauchy_integral, dims = 2) + sum(recurrence_terms, dims = 2))) .|> real
        fk .+= h

        push!(ηs, η)

        det_prod = exp.(-dets_accum .+ w_det_accum)
        # print("det_prod: "); display(det_prod)
        error_indicator = error_integral_constant * dot(weights, absfval .* det_prod .* norm_val)

        if isnan(error_indicator)
            throw(OverflowError("Error bound is NaN, choose a different value for `w`"))
        end

        if abs(real(error_indicator)) < tol
            @info "Error bound met"
            @info k
            return fk
        end
    end
    @info "Max Restarts"
    return fk
end

"""
Implementation of an explicit restarted KSM using an error indicator adapted from Chen et al.

Here we contruct the extended block diagonal matrix `That` like in `krylov_approx`. We evaluate the error using quadrature. Because we have access to the full matrix `That` we are able to rebuild our contour every restart based on the current Ritz values of `That`. This allows us to use a lower `order`.
For this algorithm we require that `A` be Hermitian as the error bound is only valid for the Lanczos decomposition.
The adaptation for the restarted scheme was derived assuming that we have a decomposition as given by Eiermann and Ernst (2006).
"""
function krylov_approx_chen_explicit(
        f, A, b, m, w;
        tol = 1.0e-16,
        max_restarts = 200,
        contour_safety = 2.0,
        order = 20,
        rebuild_contour = true
    )

    @assert ishermitian(A)

    β = norm(b)
    qm = (1 / β) * b

    (Q, T, η_prev, qm) = lanczos(A, qm, m)

    That = T

    R, c, nodes, weights, fval, absfval, norm_val, error_integral_constant = build_contour_chen(f, T, contour_safety, order, w)

    fk = β * Q * f(That)[:, 1]

    for k in 2:max_restarts

        (Q, T, η, qm) = lanczos(A, qm, m)

        That = _build_Hhat(That, T, η_prev)

        η_prev = η

        @views h = f(That)[((k - 1) * m + 1):((k - 1) * m + m), 1]

        fk .+= β * (Q * h)

        if rebuild_contour
            R, c, nodes, weights, fval, absfval, norm_val, error_integral_constant = build_contour_chen(f, That, contour_safety, order, w)
        end

        norm_hwz_T = Array{Float64}(undef, order)
        det_w = first(logabsdet(That - w * I))
        for (i, node) in enumerate(nodes)
            norm_hwz_T[i] = det_w - first(logabsdet(That - node * I))
        end

        error_indicator = (error_integral_constant * dot(weights, absfval .* exp.(norm_hwz_T) .* norm_val))

        if isnan(error_indicator)
            throw(OverflowError("Error bound is NaN, choose a different value for `w`"))
        end

        error_indicator = error_indicator |> real |> abs

        if error_indicator < tol
            @info k, error_indicator
            return fk
        end
    end
    @info "Max Restarts"
    return fk
end
