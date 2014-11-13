from Scientific.Functions.LeastSquares import *
from Numeric import log

def _cash(model, parameters, data):
    n_param = len(parameters)
    chi_sq = 0.
    alpha = Numeric.zeros((n_param, n_param))
    for point in data:
	n_i = (point[1]/point[2])**2
	scale = point[1]/n_i
	if len(point) == 3:
	    sigma = point[2]
	f = model(parameters, point[0])
	m_i = f/scale
	chi_sq = chi_sq + (2*m_i-2*n_i*log(m_i))  #calculate Cash stat
        #calculate change in cash for each derivative
        d = Numeric.array((2*m_i-2*n_i*log(m_i))[1])/scale
	alpha = alpha + d[:,Numeric.NewAxis]*d
    print chi_sq[0],parameters
    return chi_sq, alpha

def CashFit(model, parameters, data, max_iterations=2000,
                    stopping_limit = 0.000):
    """As LeastSquaresFit, but using Cash statistics for Poissonian
    noise. All bins must be positive (no background subtraction).
    """
    if len(data[0]) != 3:
        raise TypeError("Require error column in input data for Cash fitting""")
    n_param = len(parameters)
    p = ()
    i = 0
    for param in parameters:
	p = p + (DerivVar(param, i),)
	i = i + 1
    id = Numeric.identity(n_param)
    l = 0.001
    chi_sq, alpha = _cash(model, p, data)
    niter = 0
    while 1:
	delta = LA.solve_linear_equations(alpha+l*Numeric.diagonal(alpha)*id,
					  -0.5*Numeric.array(chi_sq[1]))
	next_p = map(lambda a,b: a+b, p, delta)
	next_chi_sq, next_alpha = _cash(model, next_p, data)
	if next_chi_sq > chi_sq:
	    l = 10.*l
	else:
	    l = 0.1*l
	    #if chi_sq[0] - next_chi_sq[0] < stopping_limit: break
	    p = next_p
	    chi_sq = next_chi_sq
	    alpha = next_alpha
        niter = niter + 1
        if max_iterations is not None and niter == max_iterations:
            print "Iterations exceeded"
            break
    return map(lambda p: p[0], next_p), next_chi_sq[0]

