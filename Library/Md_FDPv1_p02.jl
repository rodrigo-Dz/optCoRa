# Brink motif feedback (v01)
#   with Y inducing B synthesis,

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

    # ODE system with feedback
    function FB(du, u, p, t)
        Y, U, W, C, Y0, Y1= u
        g, mY, gY, mU, gU, mW,kD, gW, e0, eP, eM, m0, m1, k1, mUs = p
    
		du[1] = (mY * Y1)            - ((g + gY) * Y)
		du[2] = (mU * Y)            - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		du[3] = (mW * (W/(W + kD))) - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		du[4] =                     - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
		du[5] = (m0 * W)              - ((g + gY) * Y0)
		du[6] = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)

	end

    # ODE system without feedback
    function nFB(du, u, p, t)
        Y, U, W, C, Y0, Y1 = u
        g, mY, gY, mU, gU, mW,kD, gW, e0, eP, eM,m0, m1, k1, mUs = p
    
		du[1]  = (mY * Y1)            - ((g + gY) * Y)
		du[2]  =    mUs              - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		du[3]  = (mW * (W/(W + kD))) - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		du[4]  =                     - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
		du[5] = (m0 * W)              - ((g + gY) * Y0)
		du[6] = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)
    end
    

	# Define system's output (total Y):
	function out_FB(ss)
		return ss[1];
	end;
	function out_NF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mUs] = p[:mU] * ss[1];
	end;
end