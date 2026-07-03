"""
This script provides a sample symmetric problem for testing algorithms from `KrylovRestart`

Author: Errol Bucy
Date: 17.06.26
"""

using Pkg
Pkg.activate("KSMdev", shared = true)

using LinearAlgebra
using SparseArrays
using Tullio
using DataFrames
using KrylovKit
using UnicodePlots

using KrylovRestart

@info "Done loading packages"

"""
∂ₜu - Δu = 0 in Ω = (0,1)², t > 0
u(x,t) = 0 on Γ = ∂Ω,
t > 0 u(x,0) = u₀(x) in Ω

Heat equation discretized by the five-point stencil on a uniform grid involving
`n_1` interior grid points in each direction
"""

(⊗)(A, B) = kron(A, B)

n1 = 50

build_Lh(n1) = begin
    T = SymTridiagonal(fill(-2.0, n1), fill(1.0, n1 - 1))
    Id = I(n1) |> sparse

    #Lh = 1/h^2 (T⊗In⊗In + In⊗T⊗In + In⊗In⊗T) ∈ R^(n^3 × n^3).

    (⊗)(A, B) = kron(A, B)

    h = 1 / (n1 + 1)

    Lh = (1 / h^2) * (T ⊗ Id ⊗ Id + Id ⊗ T ⊗ Id + Id ⊗ Id ⊗ T)
end

abs_err_test(exact::AbstractVector, result::AbstractVector) = norm(exact - result)
rel_err_test(exact::AbstractVector, result::AbstractVector) = abs_err_test(exact, result) / norm(exact)

function eigenfuncs_initial_vec(n1::Int)
    h = 1 / (n1 + 1)

    I = collect(1:n1)
    Ip = collect(1:n1)

    # Precompute the sine matrices
    Sx = @. sin(I .* Ip' * π * h)   # n1 × n1
    Sy = @. sin(I .* Ip' * π * h)
    Sz = @. sin(I .* Ip' * π * h)

    W = [1.0 / (ip + jp + kp) for ip in Ip, jp in Ip, kp in Ip]

    @tullio u[i, j, k] :=
        Sx[i, ip] * Sy[j, jp] * Sz[k, kp] * W[ip, jp, kp]

    return reshape(u, :, 1)
end

@info "Building initial vector"

u0 = eigenfuncs_initial_vec(n1)[:]

@info "Building Exact Answer"

begin
    Lh = build_Lh(n1)
    t = 0.2
    tLh = t * Lh
    exact, conv_info = exponentiate(Lh, t, u0, issymmetric = true, maxiter = 2000, tol = 1.0e-16)
end

@info "KrylovKit exponentiate: $conv_info"

make_trace() = KrylovRestart.Trace(0, nothing, Dict{String, AbstractVector}())

function final_trace_error(trace::KrylovRestart.Trace)
    for key in ("error_indicator", "abs_err", "up_bnd", "update_norm", "low_bnd")
        if haskey(trace.data, key) && !isempty(trace.data[key])
            return (metric = key, value = trace.data[key][end])
        end
    end
    return (metric = nothing, value = nothing)
end

function test_krylov(name::String, f, args...; kwargs...)
    trace = make_trace()
    res = try
        f(args...; kwargs..., trace = trace)
    catch err
        println(name)
        Base.showerror(stdout, err)
        println()
        Base.show(stdout, MIME"text/plain"(), trace)
        println()
        rethrow()
    end
    errors = (
        rel_err = rel_err_test(exact, res),
        abs_err = abs_err_test(exact, res),
        final_trace_error = final_trace_error(trace),
    )

    println(name)
    Base.show(stdout, MIME"text/plain"(), errors)
    println()
    Base.show(stdout, MIME"text/plain"(), trace)
    println()

    return (errors = errors, trace = trace, result = res)
end

function test_algs()
    r = bestapprox_expm_data(16)

    test_krylov("krylov_approx", krylov_approx, exp, tLh, u0, 40)
    test_krylov("krylov_approx_rational", krylov_approx, r, tLh, u0, 40)
    test_krylov(
        "krylov_approx_chen_explicit",
        krylov_approx_chen_explicit,
        exp,
        tLh,
        u0,
        40,
        -10_000.0;
        contour_safety = 2.0
    )
    # test_krylov("krylov_approx_chen_implicit", krylov_approx_chen_implicit, exp, tLh, u0, 40, -6400.0)

    return nothing
end
