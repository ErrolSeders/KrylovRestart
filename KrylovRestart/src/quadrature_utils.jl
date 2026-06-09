"""
     fetch_contour_shape(H, contour_safety)

Find a radius `R` and center `c` around the Ritz values of `H`
with a loose tolerance defined by `contour_safety`.

"""
function fetch_contour_shape(H, contour_safety)
    vals = eigvals(H)
    c = sum(vals) / length(vals)
    dist = maximum(abs.(c .- vals))
    R = contour_safety * dist
    R = max(contour_safety, R)
    return R, c
end

function fetch_contour_shape(vals::AbstractVector, contour_safety)
    c = sum(vals) / length(vals)
    dist = maximum(abs.(c .- vals))
    R = contour_safety * dist
    R = max(contour_safety, R)
    return R, c
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
    expθ = exp.(im .* θ)

    x = c .+ R .* expθ
    w = 0.5 .* ω .* R .* expθ

    return x, w
end
