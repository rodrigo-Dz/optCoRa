# Antithetic feedback (v01)
#   with inactive W in complex form

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

    # ODE system with feedback
    function FB(du, u, p, t)
        Y, U, W, C = u
        g, mY, gY, mU, gU, mW, gW, e0, eP, eM, mUs = p
    
        du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
        du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
        du[3] = mW       - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
        du[4] =          - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
    end
    
	# ODE system without feedback
    function nFB(du, u, p, t)
        Y, U, W, C = u
        g, mY, gY, mU, gU, mW, gW, e0, eP, eM, mUs = p
    
        du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
        du[2] = (mUs)    - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
        du[3] = mW       - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
        du[4] =          - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
    end

	# Define system's output (total Y):
	function out_FB(SS)
        return SS[1]
	end

	function out_nFB(SS)
		return SS[1]
	end

	# Define locally analogous system:
	function localNF(p, SS)
		p[:mUs] = p[:mU] * SS[1]
	end

end
