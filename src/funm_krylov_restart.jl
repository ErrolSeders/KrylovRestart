"""
	Params(restart_length, max_restarts, stop_accuracy, bound, exact, min_decay)
using Base: update_stackframes_callback
Parameters for the restarted Krylov subspace methods:  
`restart_length`: Number of iterates in each arnoldi/lanczos run  
`max_restarts`: Maximum number of restarts  
`stop_accuracy`: Accuracy tolerance used to determine stopping, what this parameter means depends on stopping condition.  
`bound`: Compute upper and lower error bounds and stop according to their value.  
`exact`: Provided exact value of f(A)b. If provided then absolute error between iterates and exact value is used as the stopping condition  
`min_decay`: Stop if the linear convergence rate of the upper error bound or absolute error exceeds this number.  

If `bound` is `false` and `exact` is `nothing` then the algorithm is stopped when the norm of iterate updates is less than `stop_accuracy`
"""
struct Params{T<:AbstractFloat,S<:Int}
    restart_length
    max_restarts
    stop_accuracy
    bound
    exact
    min_decay
    function Params{T,S}(rl, mr, sa, bnd, ex, md) where {T<:AbstractFloat,S<:Int}
        if !isnothing(ex) && bnd
            @warn "Both exact solution provided and bound set true!"
        end

        new(rl, mr, sa, bnd, ex, md)
    end
end

Params(
    rl::S,
    mr::S,
    sa::T,
    bnd::Bool,
    ex::Union{Nothing,AbstractVector},
    md::T
) where {T<:AbstractFloat,S<:Int} =
    Params{T,S}(rl, mr, sa, bnd, ex, md)


@enum StopCode begin
    MaxRestarts #0
    AbsErrAcc #1
    UpBndAcc #2
    AbsErrLinConv #3
    UpBndLinConv #4
    UpdateAcc #5
end

mutable struct Results{T<:AbstractFloat}
    fk::AbstractVector
    num_restarts::Int64
    lowbnd::AbstractVector{T}
    upbnd::AbstractVector{T}
    errs::AbstractVector{T}
end

Results{T}(fk_init::AbstractVector) where {T<:AbstractFloat} = Results(
    fk_init,
    0,
    Vector{T}(undef, 0),
    Vector{T}(undef, 0),
    Vector{T}(undef, 0)
)

message(::Val{MaxRestarts}) = "Maximum number of restarts reached."
message(::Val{AbsErrAcc}) = "Absolute error below stopping accuracy."
message(::Val{UpBndAcc}) = "Upper bound below stopping accuracy."
message(::Val{AbsErrLinConv}) = "Linear convergence rate of absolute error greater than minimum decay."
message(::Val{UpBndLinConv}) = "Lienar convergence rate upper bound greater than minimum decay."
message(::Val{UpdateAcc}) = "Norm of updates below stopping accuracy."
message(::Val{M}) where {M} = "Invalid stop code encountered! This should never happen!"
message(c::StopCode) = message(Val(c))

"""
Checks whether the entries of v are of the form 1,...,1,0,...0,1,1 (i.e., last two entries are true/1, but not all entries are true/1)
"""
function check_stop_condition(v)
    if length(v) < 2
        return false
    end

    if all(v[end-1:end]) && !all(v)
        return true
    end

    return false
end

"""
Check whether to stop according to the provided parameters.
If `exact` is provided then absolute error is used.
Then, if `bound` is true the upper error bound estimate is used.
Lastly, if neither is set, the norm of the iteration updates is used.
"""
function stop_conditions(A, α1, β, h, qm, res::Results, p::Params)::Union{StopCode,Nothing}

    not_using_exact = isnothing(p.exact) || isempty(p.exact)

    if p.bound

        append!(res.lowbnd, β * abs(h[end-1]))
        w = A * qm
        append!(res.upbnd, β * norm((h[end-1] - α1 * h[end]) * qm + h[end] * w))

        if length(res.upbnd) > 2 && check_stop_condition(
            res.upbnd[2:end] ./ res.upbnd[1:end-1] .> p.min_decay
        )
            return UpBndLinConv
        end
        if res.upbnd[end] < p.stop_accuracy
            return UpBndAcc
        end

    elseif !not_using_exact
        append!(res.errs, norm(res.fk - p.exact))

        if length(res.errs) > 2 && check_stop_condition(res.errs[2:end] ./ res.errs[1:end-1] .> p.min_decay)
            return AbsErrLinConv
        end

        if res.errs[end] < p.stop_accuracy
            return AbsErrAcc
        end
    elseif β * norm(h) < p.stop_accuracy
        return UpdateAcc
    end

    return nothing
end

function update_alphas(α1, α2, H)

    if size(H) == (1, 1)
        μ = λ = H[1]
    else
        s,_ = eigen(H)
        μ = s .|> real |> minimum
        λ = s .|> real |> maximum
    end
    return min(real(μ), α1), max(real(λ), α2)
end

"""
	funm_krylov_restart_full(f, A, b, p)


Compute ``f(A)b`` in the manner corresponding to Alg. 1 in the paper \\

The function `f` must be callable for Matrix arguments.
"""
function funm_krylov_restart_full(f, A::AbstractArray, b::AbstractVector, p::Params)

    if p.bound
        α1 = prevfloat(Inf)
        α2 = nextfloat(-Inf)
    else
        α1 = 0.0
        α2 = 0.0
    end

    res = Results{Float64}(zeros(eltype(A), length(b)))

    m = p.restart_length

    β = norm(b)
    qm = (1 / β) * b

    Hhat = Matrix{eltype(A)}(undef, m, m)
    η_prev = NaN64

    is_A_hermitian = ishermitian(A)

    for k ∈ 1:p.max_restarts

        res.num_restarts = k
        
        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        if p.bound
            α1, α2 = update_alphas(α1, α2, H)
        end

        if k == 1
            Hhat = H
        else
            km = size(Hhat)[1]
            Hhat_expand = zeros(eltype(Hhat), km + m, km + m)
            Hhat_expand[1:km, 1:km] .= Hhat
            Hhat_expand[km+1:km+m, km+1:km+m] .= H
            Hhat_expand[((k-1)*m)+1, (k-1)*m] = η_prev
            Hhat = Hhat_expand
        end
        η_prev = η

        @views h = f(Hhat)[(k-1)*m+1:k*m, 1]
        res.fk .+= β * (Q * h)

        stop = stop_conditions(A, α1, β, h, qm, res, p)

        if stop !== nothing
            return res.fk, stop, res
        end
    end
    res.fk, MaxRestarts, res
end

"""
	funm_krylov_restart_ra(r, A, b, p)


Compute ``r(A)b ≈ f(A)b`` in the manner corresponding to Alg. 2 in the paper \\
"""
function funm_krylov_restart_ra(r::RationalApproximation, A::AbstractArray, b::AbstractVector, p::Params)

    real_res = false
    
    if (eltype(A) <: AbstractFloat) && (eltype(b) <: AbstractFloat)
        real_res = true        
    end
    
    m = p.restart_length

    if p.bound
        # Set to 'virtually' Inf so that LU does not fail for '\'
        α1 = prevfloat(Inf)
        α2 = nextfloat(-Inf)
    else
        α1 = 0.0
        α2 = 0.0
    end

    if real_res
        fk_init = r.absterm * b
    else
        fk_init = complex.(r.absterm *b)
    end
    
    res = Results{Float64}(fk_init)


    β = norm(b)
    qm = (1 / β) * b

    nr_single = length(r.single_poles)
    poles = [r.single_poles; r.conj_poles]
    coeff = [r.single_coeff; r.conj_coeff]
    Bbar = zeros(ComplexF64, m + 2, length(poles))
    Bbar[end-2, :] .= 1.0 + 0.0im
    s = 1.0

    e1 = zeros(eltype(A), m + 2)
    e1[1] = 1.0

    # H is allocated once
    Hbar = zeros(eltype(A), m + 2, m + 2)

    is_A_hermitian = ishermitian(A)

    for k ∈ 1:p.max_restarts

        res.num_restarts = k
        
        (Q, H, η, qm) = is_A_hermitian ? lanczos(A, qm, m) : arnoldi(A, qm, m)

        if p.bound
            α1, α2 = update_alphas(α1, α2, H)
        end

        # H is overwritten in this block
        Hbar[1:m, 1:m] .= H
        Hbar[m+1:m+2, m+1:m+2] = [α1 0.0; 1.0 α2]
        Hbar[m+1, m] = η

        qbar = zeros(ComplexF64,m+2)
        for p ∈ eachindex(poles)
            qbar[1] = s * Bbar[m, p]
            Bbar[1:m+2, p] = (poles[p] * I - Hbar) \ qbar
        end
        
        # contribution from single poles
        @views h = nr_single > 0 ? Bbar[1:m+2, 1:nr_single] * coeff[1:nr_single] : zeros(m+2)
        
        if nr_single < length(poles) #add in the contrabution from the conjugate poles
            @views h .+= 2 * real(Bbar[1:m+2, nr_single+1:end] * coeff[nr_single+1:end])
        end
         
        s = η

        if real_res
            h .|> real
        end
    
        @views res.fk .+= β * (Q * h[1:m])

        stop = stop_conditions(A, α1, β, h, qm, res, p)

        if stop != nothing
            return res.fk, stop, res
        end
    end
    return res.fk, MaxRestarts, res
end

