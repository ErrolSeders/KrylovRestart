@enum StopCode begin
    MaxRestarts #0
    AbsErrAcc #1
    UpBndAcc #2
    AbsErrLinConv #3
    UpBndLinConv #4
    UpdateAcc #5
end

message(::Val{MaxRestarts}) = "Maximum number of restarts reached."
message(::Val{AbsErrAcc}) = "Absolute error below stopping accuracy."
message(::Val{UpBndAcc}) = "Upper bound below stopping accuracy."
message(::Val{AbsErrLinConv}) = "Linear convergence rate of absolute error greater than minimum decay."
message(::Val{UpBndLinConv}) = "Linear convergence rate upper bound greater than minimum decay."
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

    if all(v[(end - 1):end]) && !all(v)
        return true
    end

    return false
end

"""Tracks the 'last-two-true but not all true' decay pattern without storing full history."""
mutable struct DecayTracker{T <: AbstractFloat}
    prev_metric::T
    have_prev::Bool
    last_ratio_gt::Bool
    prev_ratio_gt::Bool
    seen_ratio_false::Bool
end

DecayTracker{T}() where {T <: AbstractFloat} = DecayTracker(T(NaN), false, false, false, false)
DecayTracker(::Type{T} = Float64) where {T <: AbstractFloat} = DecayTracker{T}()

@inline function update_decay!(dt::DecayTracker{T}, metric::T, min_decay::T) where {T <: AbstractFloat}
    if !dt.have_prev
        dt.prev_metric = metric
        dt.have_prev = true
        return false
    end

    ratio_gt = (metric / dt.prev_metric) > min_decay
    dt.prev_metric = metric

    dt.prev_ratio_gt = dt.last_ratio_gt
    dt.last_ratio_gt = ratio_gt

    if !ratio_gt
        dt.seen_ratio_false = true
    end

    return dt.prev_ratio_gt && dt.last_ratio_gt && dt.seen_ratio_false
end

"""Internal stopping logic (and optional logging)."""
function stop_conditions!(
        A,
        α1,
        β,
        h,
        qm,
        fk,
        m::Int,
        bound::Bool,
        exact,
        tol::T,
        min_decay::T,
        trace,
        callback,
        k::Int,
        up_decay::DecayTracker{T},
        err_decay::DecayTracker{T}
    )::Union{StopCode, Nothing} where {T <: AbstractFloat}

    not_using_exact = exact === nothing || isempty(exact)

    if bound
        low = T(β * abs(h[end - 1]))
        w = A * qm
        up = T(β * norm((h[end - 1] - α1 * h[end]) * qm + h[end] * w))

        log_lowbnd!(trace, low)
        log_upbnd!(trace, up)

        stop = nothing
        if update_decay!(up_decay, up, min_decay)
            stop = UpBndLinConv
        elseif up < tol
            stop = UpBndAcc
        end

        maybe_callback(callback, (; k, lowbnd = low, upbnd = up, stop))
        return stop

    elseif !not_using_exact
        err = T(norm(fk - exact))
        log_err!(trace, err)

        stop = nothing
        if update_decay!(err_decay, err, min_decay)
            stop = AbsErrLinConv
        elseif err < tol
            stop = AbsErrAcc
        end

        maybe_callback(callback, (; k, err, stop))
        return stop

    else
        update_norm = T(β * norm(@view h[1:m]))
        log_update_norm!(trace, update_norm)

        stop = update_norm < tol ? UpdateAcc : nothing
        maybe_callback(callback, (; k, update_norm, stop))
        return stop
    end
end

"""Optional collector for convergence statistics (one entry per restart)."""
mutable struct KrylovTrace{T <: AbstractFloat}
    num_restarts::Int
    stop::Union{Nothing, StopCode}
    lowbnd::Vector{T}
    upbnd::Vector{T}
    errs::Vector{T}
    update_norms::Vector{T}
end

KrylovTrace{T}() where {T <: AbstractFloat} = KrylovTrace(0, nothing, T[], T[], T[], T[])
KrylovTrace(::Type{T} = Float64) where {T <: AbstractFloat} = KrylovTrace{T}()

function reset!(tr::KrylovTrace)
    tr.num_restarts = 0
    tr.stop = nothing
    empty!(tr.lowbnd)
    empty!(tr.upbnd)
    empty!(tr.errs)
    empty!(tr.update_norms)
    return tr
end

set_restart!(::Nothing, _k) = nothing
set_restart!(tr::KrylovTrace, k) = (tr.num_restarts = k)

set_stop!(::Nothing, _stop::StopCode) = nothing
set_stop!(tr::KrylovTrace, stop::StopCode) = (tr.stop = stop)

log_lowbnd!(::Nothing, _x) = nothing
log_lowbnd!(tr::KrylovTrace, x) = push!(tr.lowbnd, x)

log_upbnd!(::Nothing, _x) = nothing
log_upbnd!(tr::KrylovTrace, x) = push!(tr.upbnd, x)

log_err!(::Nothing, _x) = nothing
log_err!(tr::KrylovTrace, x) = push!(tr.errs, x)

log_update_norm!(::Nothing, _x) = nothing
log_update_norm!(tr::KrylovTrace, x) = push!(tr.update_norms, x)

maybe_callback(::Nothing, _info) = nothing
maybe_callback(cb, info) = cb(info)
