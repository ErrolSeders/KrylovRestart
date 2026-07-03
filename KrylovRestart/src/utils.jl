function unit_vector(T::Type, size::Int, pos::Int)
    unit = zeros(T, size)
    unit[pos] = one(T)
    return unit
end

@enum StopCode begin
    MaxRestarts #0
    AbsErrAcc #1
    UpBndAcc #2
    AbsErrLinConv #3
    UpBndLinConv #4
    UpdateAcc #5
    IndicatorAcc #6
end

message(::Val{MaxRestarts}) = "Maximum number of restarts reached."
message(::Val{AbsErrAcc}) = "Absolute error below stopping accuracy."
message(::Val{UpBndAcc}) = "Upper bound below stopping accuracy."
message(::Val{AbsErrLinConv}) = "Linear convergence rate of absolute error greater than minimum decay."
message(::Val{UpBndLinConv}) = "Linear convergence rate upper bound greater than minimum decay."
message(::Val{UpdateAcc}) = "Norm of updates below stopping accuracy."
message(::Val{IndicatorAcc}) = "Error Indicator below stopping accuracy"
message(::Val{M}) where {M} = "Invalid stop code encountered! This should never happen!"
message(c::StopCode) = message(Val(c))

"""Optional collector for convergence statistics (one entry per restart)."""
mutable struct Trace
    num_restarts::Int
    stop::Union{Nothing, StopCode}
    data::Dict{String, AbstractVector}
end

function Base.show(io::IO, tr::Trace)
    print(io, "Trace(restarts=$(tr.num_restarts), stop=$(isnothing(tr.stop) ? "nothing" : tr.stop), metrics=$(length(tr.data)))")
end

function Base.show(io::IO, ::MIME"text/plain", tr::Trace)
    println(io, "Trace")
    println(io, "  restarts: ", tr.num_restarts)
    if isnothing(tr.stop)
        println(io, "  stop: not set")
    else
        println(io, "  stop: ", tr.stop, " — ", message(tr.stop))
    end

    if isempty(tr.data)
        print(io, "  metrics: (none)")
        return
    end

    println(io, "  metrics:")
    for key in sort!(collect(keys(tr.data)))
        values = tr.data[key]
        n = length(values)
        print(io, "    - ", key, ": ", n, " entr", n == 1 ? "y" : "ies")
        if n > 0
            print(io, " (last=", values[end], ")")
        end
        println(io)
    end
end

function reset!(tr::Trace)
    tr.num_restarts = 0
    tr.stop = nothing
    empty!(tr.data)
    return tr
end

set_restart!(::Nothing, _k) = nothing
set_restart!(tr::Trace, k) = (tr.num_restarts = k)

set_stop!(::Nothing, _stop::StopCode) = nothing
set_stop!(tr::Trace, stop::StopCode) = (tr.stop = stop)

log_item!(::Nothing, item_name, _x) = nothing
log_item!(tr::Trace, item_name, value) = push!(get!(tr.data, item_name, typeof(value)[]), value)

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
mutable struct DecayTracker
    prev_metric
    have_prev::Bool
    last_ratio_gt::Bool
    prev_ratio_gt::Bool
    seen_ratio_false::Bool
end

DecayTracker() = DecayTracker(NaN, false, false, false, false)

@inline function update_decay!(dt::DecayTracker, metric, min_decay)
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
        tol,
        min_decay,
        trace::Union{Trace, Nothing},
        k::Int,
        up_decay::DecayTracker,
        err_decay::DecayTracker
    )::Union{StopCode, Nothing}

    not_using_exact = exact === nothing || isempty(exact)

    if bound
        low = β * abs(h[end - 1])
        w = A * qm
        up = β * norm((h[end - 1] - α1 * h[end]) * qm + h[end] * w)

        log_item!(trace, "low_bnd", low)
        log_item!(trace, "up_bnd", up)

        stop = nothing
        if update_decay!(up_decay, up, min_decay)
            stop = UpBndLinConv
        elseif up < tol
            stop = UpBndAcc
        end

        return stop

    elseif !not_using_exact
        err = norm(fk - exact)
        log_item!(trace, "abs_err", err)

        stop = nothing
        if update_decay!(err_decay, err, min_decay)
            stop = AbsErrLinConv
        elseif err < tol
            stop = AbsErrAcc
        end

        return stop

    else
        update_norm = β * norm(h[1:m])
        log_update_norm!(trace, "update_norm", update_norm)

        stop = update_norm < tol ? UpdateAcc : nothing

        return stop
    end
end
