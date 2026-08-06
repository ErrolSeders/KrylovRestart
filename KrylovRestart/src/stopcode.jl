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
