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
        g, mY, gY, mU, gU, mW, gW, e0, eP, eM, m0, m1, mP, kP, mUs = p
    
        du[1] = (mY * W)  + (mP * (Y/(Y + kP))) - ((g + gY) * Y)
        du[2]  = (mU * Y)                          - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
        du[3]  =    mW                             - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
        du[4]  =                                   - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
    end
    
	# ODE system without feedback
    function nFB(du, u, p, t)
        Y, U, W, C = u
        g, mY, gY, mU, gU, mW, gW, e0, eP, eM, m0, m1, mP, kP, mUs = p
    
        du[1] = (mY * W)  + (mP * (Y/(Y + kP)))  - ((g + gY) * Y)
        du[2]  = (mUs)                          - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
        du[3]  =    mW                             - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
        du[4]  =                                   - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
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
