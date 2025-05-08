using DifferentialEquations, BifurcationKit
using NLsolve

# Define the ODE system
definition!(du, u, p, t) = begin
    Y, U, W, C = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM = p
    
    du[1] = (mY * W) - ((g + gY) * Y)
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
    du[3] =    mW    - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
    du[4] =          - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
end

# Parameter values
params = [0.01, 1.0, 0.1, 0.0152, 0.05, 0.10, 0.0001, 0.0001, 0.0375, 0.5]  # g, mY, gY, ...

# Initial conditions
u0 = [0.0, 0.0, 0.0, 0.0]

# Time span
tspan = (0.0, 1000.0)

# Solve ODEs
prob = ODEProblem(definition!, u0, tspan, params)
sol = DifferentialEquations.solve(prob, Tsit5())

# Function to find steady states
function find_equilibrium(params)
    f(u) = definition!(zeros(4), u, params, 0.0);  # Residual function
    return nlsolve(f, [1.0, 1.0, 1.0, 1.0]).zero
end

equilibrium = find_equilibrium(params)

# Two-parameter bifurcation analysis
br = bifurcationdiagram(prob, (@lens _[2]), 
    (@lens _[4]), range(0.5, 2.0, length=50), range(0.005, 0.03, length=50))  # Varying mY (0.5 to 2.0) and mU (0.005 to 0.03)

plot(br)
