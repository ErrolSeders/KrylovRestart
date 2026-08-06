"""
     fetch_contour_shape(H, contour_safety)

Find a radius `R` and center `c` around the Ritz values of `H`
using Base: Ordered
with a loose tolerance defined by `contour_safety`.

"""
function fetch_contour_circle(H, contour_safety)
    vals = eigvals(H)
    c = sum(vals) / length(vals)
    dist = maximum(abs.(c .- vals))
    R = contour_safety * dist
    R = max(contour_safety, R)
    return R, c
end

function fetch_contour_circle(vals::AbstractVector, contour_safety)
    c = sum(vals) / length(vals)
    dist = maximum(abs.(c .- vals))
    R = contour_safety * dist
    R = max(contour_safety, R)
    return R, c
end

"""
Construct quadrature nodes and weights for a circular contour integral for an evenly spaced trapezoidal rule.
"""

function contour_quad_nodes_trapezoid(order, R, c = 0.0 + 0.0im)

    h = 2π / (order - 1)
    nodes = [i for i in 0:h:2π]

    cisnodes = cis.(nodes)

    x = c .+ (R .* cisnodes)
    w = (R / order) .* cisnodes

    return x, w
end

"""
    contour_quad_nodes(order,R,[c=0.0+0.0im])

Construct Gauss-Legendre quadrature nodes and weights for
a circular controur intregral

From Γ := t(θ) = c + R*exp(im*θ)
To x ∈ [-1,1]

"""
function contour_quad_nodes(order, R, c = 0.0 + 0.0im)

    ξ, ω = gausslegendre(order)

    θ = π .* ξ
    expθ = θ .|> cis

    x = c .+ R .* expθ
    w = 0.5 .* ω .* R .* expθ

    return x, w
end
