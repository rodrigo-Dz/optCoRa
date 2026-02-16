# Brink motif feedback (v01)
#   with Y inducing B synthesis,

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

    # ODE system with feedback
    function FB(du, u, p, t)
        Y, A, B, C, U, Us = u
        g, mY, gY, mA, kD, mB, mU, e0, eP, bA, bI, mAs = p
    
		du[1]  = (mY * U)              - ((g + gY) * Y)
		du[2]  = (mA * (kD/(kD + Y)))  - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		du[3]  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		du[4]  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		du[5]  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		du[6]  =                       - (g * Us) - (bA * A * Us) + (bI * B * U)
    end

    # ODE system without feedback
    function nFB(du, u, p, t)
        Y, A, B, C, U, Us = u
        g, mY, gY, mA, kD, mB, mU, e0, eP, bA, bI, mAs = p
    
		du[1]  = (mY * U)              - ((g + gY) * Y)
		du[2]  =  mAs                  - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		du[3]  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		du[4]  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		du[5]  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		du[6] =                       - (g * Us) - (bA * A * Us) + (bI * B * U)
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
		p[:mAs] = p[:mA] * (p[:kD]/(p[:kD] + ss[1]));
	end;
end