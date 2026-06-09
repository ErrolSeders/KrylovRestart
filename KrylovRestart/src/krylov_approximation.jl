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

    Hhat = Matrix{eltype(A)}(undef, 0, 0)
    η_prev = NaN

    is_A_hermitian = ishermitian(A)

    for k in 1:max_restarts

        set_restart!(trace, k)

        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        if bound
            α1, α2 = update_alphas(α1, α2, H)
        end

        if k == 1
            Hhat = zeros(eltype(H), m + 2, m + 2)
            Hhat[1:m, 1:m] = H
            Hhat[(m + 1):(m + 2), (m + 1):(m + 2)] = [α1 0.0; 1.0 α2]
        else
            km = size(Hhat, 1) - 2
            Hhat_expand = zeros(eltype(Hhat), km + m + 2, km + m + 2)
            Hhat_expand[1:km, 1:km] .= Hhat[1:km, 1:km]
            Hhat_expand[(km + 1):(km + m), (km + 1):(km + m)] .= H
            Hhat_expand[((k - 1) * m) + 1, (k - 1) * m] = η_prev
            Hhat_expand[km + m + 1, km + m] = η
            Hhat_expand[(km + m + 1):(km + m + 2), (km + m + 1):(km + m + 2)] = [α1 0.0; 1.0 α2]
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

    # H is allocated once
    Hbar = zeros(eltype(A), m + 2, m + 2)

    is_A_hermitian = ishermitian(A)

    for k in 1:max_restarts

        set_restart!(trace, k)

        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        if bound
            α1, α2 = update_alphas(α1, α2, H)
        end

        # H is overwritten in this block
        Hbar[1:m, 1:m] .= H
        Hbar[(m + 1):(m + 2), (m + 1):(m + 2)] = [α1 0.0; 1.0 α2]
        Hbar[m + 1, m] = η

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
function krylov_approx_barycentric(
        f::Function,
        A::AbstractArray,
        b::AbstractVector;
        m::Int = 30,
        max_restarts::Int = 200,
        tol = nothing,
        exact::Union{Nothing, AbstractVector} = nothing,
        min_decay = 0.95,
        trace::Union{Nothing, KrylovTrace} = nothing,
        callback = nothing
    )

    m < 1 && throw(ArgumentError("m must be ≥ 1"))

    Tdefault = float(real(eltype(A)))
    tolT = tol === nothing ? eps(Tdefault) : float(tol)
    TT = typeof(tolT)
    tolT = TT(tolT)
    min_decayT = TT(min_decay)

    up_decay = DecayTracker(TT)
    err_decay = DecayTracker(TT)

    real_res = (eltype(A) <: AbstractFloat) && (eltype(b) <: AbstractFloat)

    β = norm(b)
    qm = (1 / β) * b

    e1 = zeros(ComplexF64, m + 2)
    e1[1] = 1.0 + 0.0im

    # H is allocated once
    Hbar = zeros(eltype(A), m + 2, m + 2)

    is_A_hermitian = ishermitian(A)

    α1, α2 = Inf, -Inf

    Tvec = promote_type(eltype(A), eltype(b))
    fk = zeros(real_res ? Tvec : ComplexF64, size(A, 1))

    for k in 1:max_restarts

        set_restart!(trace, k)

        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        vals = eigvals(H)

        α1 = min(α1, vals .|> real |> minimum)
        α2 = max(α2, vals .|> real |> maximum)

        # Include the two extra spectral points that appear in Hbar.
        vals_aug = vcat(vals, α1, α2)
        R, c = fetch_contour_shape(vals_aug, 1.1)

        contour = Circle(c, R)

        r = approximate(f, contour)

        Z = r.fun.nodes
        W = r.fun.weights
        WF = r.fun.w_times_f

        # H is overwritten in this block
        Hbar[1:m, 1:m] .= H
        Hbar[(m + 1):(m + 2), (m + 1):(m + 2)] = [α1 0.0; 1.0 α2]
        Hbar[m + 1, m] = η

        # Evaluate the AAA barycentric rational on the projected matrix:
        #   r(M) = (∑_j w_j f_j (z_j I - M)^{-1}) (∑_j w_j (z_j I - M)^{-1})^{-1}
        # applied to the restart RHS q = e1.
        q = e1
        numer = zeros(ComplexF64, m + 2)
        denom = zeros(ComplexF64, m + 2, m + 2)
        for j in eachindex(Z)
            M = Z[j] * I - Hbar
            Fm = lu(M)
            numer .+= WF[j] .* (Fm \ q)
            denom .+= W[j] .* (Fm \ I)
        end
        h = denom \ numer

        if real_res
            h = real.(h)
        end

        @views fk .+= β * (Q * h[1:m])

        stop = stop_conditions!(
            A, α1, β, h, qm, fk, m,
            true, exact, tolT, min_decayT,
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
            ηs[j] * ((t * I - Hs[j]) \ e1)[end] for j in 1:length(ηs)
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

        @inbounds for j in eachindex(ηs)
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
        A::AbstractMatrix,
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
