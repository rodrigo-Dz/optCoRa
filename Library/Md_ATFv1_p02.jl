# Antithetic feedback (v01)
#   with inactive W in complex form

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

    # ODE system with feedback
    function FB(du, u, p, t)
        Y, U, W, C , Y0, Y1= u
        g, mY, kD, gY, mU, gU, mW, gW, e0, eP, eM, m0, m1, k1, mUs = p
    
        du[1] = (mY * (kD/(Y1 + kD))) - ((g + gY) * Y) # dY/dt
        du[2] = (mU * Y)              - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C) # dU/dt
        du[3] =    mW                 - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C) # dW/dt
        du[4] =                       - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C) # dC/dt
        du[5] = (m0 * W)              - ((g + gY) * Y0)
        du[6] = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)
    end
    
	# ODE system without feedback
    function nFB(du, u, p, t)
        Y, U, W, C , Y0, Y1= u
        g, mY, kD, gY, mU, gU, mW, gW, e0, eP, eM, m0, m1, k1, mUs = p
    
        du[1] = (mY * (kD/(Y1 + kD))) - ((g + gY) * Y) # dY/dt
        du[2] =   (mUs)               - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
        du[3] =    mW                 - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C) # dW/dt
        du[4] =                       - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C) # dC/dt
        du[5] = (m0 * W)              - ((g + gY) * Y0)
        du[6] = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)
    end

	# Define system's output (total Y):
	function out_FB(SS)
        return SS[1]
	end

	function out_nFB(SS)
		return SS[1];
	end

	# Define locally analogous system:
	function localNF(p, SS)
		p[:mUs] = p[:mU] * SS[1];
	end

end
