abstract type FunctionProperty end


"""
Mark the function as Stieltjes, this allows an integral from [-∞,0] for
algorithms dependent on quadrature. This contour is independent of the spectrum
of the matrix. 
"""
struct Stieltjes <: FunctionProperty end
