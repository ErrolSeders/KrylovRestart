struct RationalApproximation{T <: AbstractFloat, S <: Int}
    degree::S
    absterm::T
    conj_poles::Vector{Complex{T}}
    conj_coeff::Vector{Complex{T}}
    single_poles::Vector{Complex{T}}
    single_coeff::Vector{Complex{T}}
end

"""
Compute ``r(A)`` for a given rational approximation: ``r(A) ≈ f(A)``
"""
function (r::RationalApproximation)(A::AbstractArray)

    if issparse(A)
        A = A |> Matrix
    end

    n = first(size(A))
    out = fill(r.absterm, n) |> spdiagm

    for j in eachindex(r.single_poles)
        out = out + r.single_coeff[j] * inv(r.single_poles[j] * I - A)
    end

    if isreal(A)
        for j in eachindex(r.conj_poles)
            out = out + 2 * real(r.conj_coeff[j] * inv(r.conj_poles[j] * I - A))
        end
    else
        for j in eachindex(r.conj_poles)
            out = out + r.conj_coeff[j] * inv(r.conj_poles[j] * I - A)
            out = out + conj(r.conj_coeff[j]) * inv(conj(r.conj_poles[j]) * I - A)
        end
    end
    return out
end

"""
Compute ``r(A)b`` for a given rational approximation: ``r(A)b ≈ f(A)b``
"""
function (r::RationalApproximation)(A::AbstractArray, b::AbstractVector)

    out = r.absterm * b

    for j in eachindex(r.single_poles)
        out = out + r.single_coeff[j] * ((r.single_poles[j] * I - A) \ b)
    end
    if isreal(A)
        for j in eachindex(r.conj_poles)
            out = out + 2 * real(r.conj_coeff[j] * (((r.conj_poles[j] * I) - A) \ b))
        end
    else
        for j in eachindex(r.conj_poles)
            out = out + r.conj_coeff[j] * ((r.conj_poles[j] * I - A) \ b)
            out = out + conj(r.conj_coeff[j]) * ((conj(r.conj_poles[j]) * I - A) \ b)
        end
    end
    return out
end

"""
	bestapprox_expm_data(deg)

Provides the best uniform rational approximation to ``e^A`` on ``(-∞, 0]`` of degree `deg`
derived from data provided by Richard Varga et. al.

Note that the approximation can be shifted by exp(z)
"""
function bestapprox_expm_data(deg::Integer)
    bestapproximations = Dict(
        2 => (
            z = [
                -5.8479863263849063e-1 + 1im * -1.1855377812415784e+0
                -5.8479863263849063e-1 + 1im * 1.1855377812415784e+0
            ],
            kappa = [
                -1.6885872734712942e-1 + 1im * 8.0945019772077931e-1
                -1.6885872734712942e-1 + 1im * -8.0945019772077931e-1
            ],
            absterm = 7.3586701695805279e-3,
        ), 3 => (
            z = [
                -1.982084787507408e-1 + 1im * 2.4107321029593383e+0
                -1.982084787507408e-1 + 1im * -2.4107321029593383e+0
                -1.3688493230338539e+0 + 1im * 0.0e+0
            ],
            kappa = [
                -6.9116392297431128e-1 + 1im * 4.3182218073431179e-2
                -6.9116392297431128e-1 + 1im * -4.3182218073431179e-2
                1.4838485523029701e+0 + 1im * 0.0e+0
            ],
            absterm = -7.9938063633568784e-4,
        ), 16 => (
            z = [
                1.0843917078696986e+1 + 1im * 1.9277446167181651e+1
                1.0843917078696986e+1 + 1im * -1.9277446167181651e+1
                5.2649713434426522e+0 + 1im * 1.6220221473167935e+1
                5.2649713434426522e+0 + 1im * -1.6220221473167935e+1
                1.4139284624888904e+0 + 1im * 1.3497725698892721e+1
                1.4139284624888904e+0 + 1im * -1.3497725698892721e+1
                -1.4193758971857036e+0 + 1im * 1.0925363484496751e+1
                -1.4193758971857036e+0 + 1im * -1.0925363484496751e+1
                -3.5091036084148506e+0 + 1im * 8.436198985884376e+0
                -3.5091036084148506e+0 + 1im * -8.436198985884376e+0
                -4.9931747377180455e+0 + 1im * 5.9968817136039023e+0
                -4.9931747377180455e+0 + 1im * -5.9968817136039023e+0
                -5.9481522689511781e+0 + 1im * 3.5874573620183634e+0
                -5.9481522689511781e+0 + 1im * -3.5874573620183634e+0
                -6.4161776990994213e+0 + 1im * 1.1941223933701262e+0
                -6.4161776990994213e+0 + 1im * -1.1941223933701262e+0
            ],
            kappa = [
                5.0901521865241194e-7 + 1im * -2.4220017652852255e-5
                5.0901521865241194e-7 + 1im * 2.4220017652852255e-5
                -2.115174218246386e-4 + 1im * 4.3892969647380299e-3
                -2.115174218246386e-4 + 1im * -4.3892969647380299e-3
                -4.1023136835406988e-2 + 1im * -1.5743466173455728e-1
                -4.1023136835406988e-2 + 1im * 1.5743466173455728e-1
                1.479300711355878e+0 + 1im * 1.7686588323782464e+0
                1.479300711355878e+0 + 1im * -1.7686588323782464e+0
                -1.5059585270022891e+1 + 1im * -5.7514052776425739e+0
                -1.5059585270022891e+1 + 1im * 5.7514052776425739e+0
                6.2518392463209082e+1 + 1im * -1.1190391094283951e+1
                6.2518392463209082e+1 + 1im * 1.1190391094283951e+1
                -1.1339775178483859e+2 + 1im * 1.0194721704215897e+2
                -1.1339775178483859e+2 + 1im * -1.0194721704215897e+2
                6.4500878025537119e+1 + 1im * -2.245944076265209e+2
                6.4500878025537119e+1 + 1im * 2.245944076265209e+2
            ],
            absterm = 2.1248537104952236e-16,
        )
    )

    (z, kappa, absterm) = bestapproximations[deg]

    conj_poles = -z[2:2:end]
    conj_coeff = kappa[2:2:end]

    if deg % 2 == 1
        single_poles = [-real(z[end]) + 0im]
        single_coeff = [real(kappa[end]) + 0im]
    else
        single_poles = Vector{ComplexF64}(undef, 0)
        single_coeff = Vector{ComplexF64}(undef, 0)
    end

    return RationalApproximation(deg, absterm, conj_poles, conj_coeff, single_poles, single_coeff)
end
