using KrylovRestart
using LinearAlgebra
using Random
using Test

const KR = KrylovRestart

run_krylov(args...; kwargs...) = begin
    tr = KR.KrylovTrace()
    fk = KR.krylov_approx(args...; trace = tr, kwargs...)
    return fk, tr
end

@testset "KrylovRestart" begin
    Random.seed!(1234)

    @testset "Arnoldi invariants" begin
        n, m = 10, 5
        A = 0.1 .* randn(n, n)
        b = randn(n)

        Q, H, η, qnext = KR.arnoldi(A, b, m)

        @test size(Q) == (n, m)
        @test size(H) == (m, m)
        @test isapprox(Q' * Q, I(m); rtol = 0, atol = 1.0e-10)
        @test norm(Q' * qnext) ≤ 1.0e-10

        emT = zeros(eltype(A), 1, m)
        emT[end] = one(eltype(A))
        R = A * Q - Q * H - (η .* qnext) * emT
        @test norm(R) / norm(A * Q) ≤ 1.0e-10
    end

    @testset "Lanczos invariants" begin
        n, m = 12, 6
        M = randn(n, n)
        A = 0.1 .* (M + M') / 2
        b = randn(n)

        Q, T, η, qnext = KR.lanczos(A, b, m)

        @test size(Q) == (n, m)
        @test size(T) == (m, m)
        @test isapprox(Q' * Q, I(m); rtol = 0, atol = 1.0e-10)
        @test norm(Q' * qnext) ≤ 1.0e-10
        @test isapprox(T, T'; rtol = 0, atol = 1.0e-12)

        T_off = copy(T)
        for i in 1:m, j in 1:m
            if abs(i - j) ≤ 1
                T_off[i, j] = 0
            end
        end
        @test norm(T_off) ≤ 1.0e-10

        emT = zeros(eltype(A), 1, m)
        emT[end] = one(eltype(A))
        R = A * Q - Q * T - (η .* qnext) * emT
        @test norm(R) / norm(A * Q) ≤ 1.0e-10
    end

    @testset "krylov_approx(exp, A, b) accuracy" begin
        n = 10

        # Hermitian path (Lanczos)
        M = randn(n, n)
        Aherm = 0.05 .* (M + M') / 2
        b = randn(n)
        fk, tr = run_krylov(exp, Aherm, b, 5; max_restarts = 60, tol = 1.0e-12, bound = false)
        exact = exp(Matrix(Aherm)) * b
        @test tr.stop != KR.MaxRestarts
        @test norm(fk - exact) / norm(exact) ≤ 1.0e-8

        # Non-Hermitian path (Arnoldi)
        Anon = 0.05 .* randn(n, n)
        c = randn(n)
        fk2, tr2 = run_krylov(exp, Anon, c, 6; max_restarts = 80, tol = 1.0e-12, bound = false)
        exact2 = exp(Matrix(Anon)) * c
        @test tr2.stop != KR.MaxRestarts
        @test norm(fk2 - exact2) / norm(exact2) ≤ 1.0e-8
    end

    @testset "krylov_approx other matrix functions" begin
        n = 10

        M = randn(n, n)
        Aherm = 0.05 .* (M + M') / 2
        b = randn(n)

        Anon = 0.05 .* randn(n, n)
        c = randn(n)

        S = randn(n, n)
        Aspd = I + 0.05 .* (S' * S)
        d = randn(n)

        @testset "sin" begin
            fkH, trH = run_krylov(sin, Aherm, b, 6; max_restarts = 100, tol = 1.0e-12, bound = false)
            exactH = sin(Matrix(Aherm)) * b
            @test trH.stop != KR.MaxRestarts
            @test norm(fkH - exactH) / norm(exactH) ≤ 1.0e-8

            fkN, trN = run_krylov(sin, Anon, c, 6; max_restarts = 120, tol = 1.0e-12, bound = false)
            exactN = sin(Matrix(Anon)) * c
            @test trN.stop != KR.MaxRestarts
            @test norm(fkN - exactN) / norm(exactN) ≤ 1.0e-8
        end

        @testset "cos" begin
            fkH, trH = run_krylov(cos, Aherm, b, 6; max_restarts = 100, tol = 1.0e-12, bound = false)
            exactH = cos(Matrix(Aherm)) * b
            @test trH.stop != KR.MaxRestarts
            @test norm(fkH - exactH) / norm(exactH) ≤ 1.0e-8

            fkN, trN = run_krylov(cos, Anon, c, 6; max_restarts = 120, tol = 1.0e-12, bound = false)
            exactN = cos(Matrix(Anon)) * c
            @test trN.stop != KR.MaxRestarts
            @test norm(fkN - exactN) / norm(exactN) ≤ 1.0e-8
        end

        @testset "polynomial" begin
            f(A) = A^3 + 0.25 * A + I

            fkH, trH = run_krylov(f, Aherm, b, 6; max_restarts = 100, tol = 1.0e-12, bound = false)
            exactH = f(Matrix(Aherm)) * b
            @test trH.stop != KR.MaxRestarts
            @test norm(fkH - exactH) / norm(exactH) ≤ 1.0e-8

            fkN, trN = run_krylov(f, Anon, c, 6; max_restarts = 120, tol = 1.0e-12, bound = false)
            exactN = f(Matrix(Anon)) * c
            @test trN.stop != KR.MaxRestarts
            @test norm(fkN - exactN) / norm(exactN) ≤ 1.0e-8
        end

        @testset "sqrt (SPD)" begin
            fk, tr = run_krylov(sqrt, Aspd, d, 6; max_restarts = 120, tol = 1.0e-12, bound = false)
            exact = sqrt(Matrix(Aspd)) * d
            @test tr.stop != KR.MaxRestarts
            @test norm(fk - exact) / norm(exact) ≤ 1.0e-8
        end

        @testset "log (SPD)" begin
            fk, tr = run_krylov(log, Aspd, d, 6; max_restarts = 120, tol = 1.0e-12, bound = false)
            exact = log(Matrix(Aspd)) * d
            @test tr.stop != KR.MaxRestarts
            @test norm(fk - exact) / norm(exact) ≤ 1.0e-8
        end
    end

    @testset "Stopping codes smoke tests" begin
        n, m = 8, 4
        A = 0.1 .* randn(n, n)
        b = randn(n)

        # Update-based stopping should trigger immediately for huge tol
        _fk, tr = run_krylov(exp, A, b, m; max_restarts = 50, tol = 1.0e9, bound = false)
        @test tr.stop == KR.UpdateAcc

        # Bound-based stopping should trigger immediately for huge tol (needs m≥2)
        _fk2, tr2 = run_krylov(exp, A, b, 3; max_restarts = 50, tol = 1.0e9, bound = true)
        @test tr2.stop == KR.UpBndAcc

        # Exact-based stopping criterion
        exact = exp(Matrix(A)) * b
        tr3 = KR.KrylovTrace()
        _fk3 = KR.krylov_approx(exp, A, b, m; max_restarts = 50, tol = 1.0e9, bound = false, exact = exact, trace = tr3)
        @test tr3.stop == KR.AbsErrAcc
        @test !isempty(tr3.errs)
    end

    @testset "RationalApproximation" begin
        r = KR.bestapprox_expm_data(16)

        # Scalar sanity: approximation for exp(x) on (-∞, 0]
        t = 5.0
        approx = r([-t 0; 0 -t], [1.0, 1.0])[1]
        @test abs(approx - exp(-t)) ≤ 1.0e-12

        @testset "krylov_approx(r, A, b) Hermitian" begin
            n, m = 10, 6
            vals = collect(range(0.1, 2.0; length = n))
            A = -Diagonal(vals) # Hermitian with spectrum in (-∞, 0]
            @test maximum(diag(A)) ≤ 0

            b = randn(n)
            approx = r(Matrix(A), b)
            exact = exp(A) * b

            fk, tr = run_krylov(r, Matrix(A), b, m; max_restarts = 200, tol = 1.0e-12, bound = false)

            @test tr.stop != KR.MaxRestarts
            @test norm(fk - approx) / norm(approx) ≤ 1.0e-10
            @test norm(fk - exact) / norm(exact) ≤ 1.0e-10
        end

        @testset "krylov_approx(r, A, b) non-Hermitian" begin
            n, m = 10, 6
            vals = collect(range(0.1, 2.0; length = n))

            A = Matrix(-Diagonal(vals))
            A .+= 0.02 .* triu(randn(n, n), 1) # upper triangular; eigenvalues = diag(A)
            @test maximum(diag(A)) ≤ 0
            @test !ishermitian(A)

            b = randn(n)
            approx = r(A, b)
            exact = exp(A) * b

            fk, tr = run_krylov(r, A, b, m; max_restarts = 250, tol = 1.0e-12, bound = false)

            @test tr.stop != KR.MaxRestarts
            @test norm(fk - approx) / norm(approx) ≤ 1.0e-10
            @test norm(fk - exact) / norm(exact) ≤ 1.0e-10
        end
    end
end
