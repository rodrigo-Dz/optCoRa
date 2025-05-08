using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve

# Define the ODE system
function ode_system!(du, u, p, t)
    Y, U, W, C = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM = p

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)  # dU/dt
    du[3] = mW - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] = -((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end

# Parameters
params = [
    0.01,  # g
    1.0,   # mY (se variará)
    0.1,   # gY
    1, # mU (se variará)
    0.05,  # gU
    400,  # mW
    0.0001, # gW
    0.0001, # e0
    100, # eP
    0.5   # eM
]

# Initial conditions
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Time span
tspan = (0.0, 1000.0)  # Long time span to reach steady state

# Find equilibrium points
function find_equilibrium(p)
    function equilibrium_condition(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
    return result.zero
end

# Compute the Jacobian at an equilibrium point
function compute_jacobian(u_eq, p)
    # Define a function that returns the derivatives as an array
    function system_derivatives(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    ForwardDiff.jacobian(system_derivatives, u_eq)
end

# Bifurcation analysis: Vary mY and mU
mY_range = 10 .^ range(-3, 3, length=100)  # Range of mY values
mU_range = 10 .^ range(-3, 3, length=100) # Range of mU values
bifurcation_data = []

for mY in mY_range, mU in mU_range
    params[2] = mY  # Update mY
    params[6] = mU  # Update mU
    u_eq = find_equilibrium(params)
    J = compute_jacobian(u_eq, params)
    eigenvalues = eigen(J).values
    println(eigenvalues)
    max_real_eigenvalue = minimum(real.(eigenvalues))  # Track the maximum real part
    push!(bifurcation_data, (mY, mU, max_real_eigenvalue))
end

# Prepare data for heatmap
mY_values = [data[1] for data in bifurcation_data]
mU_values = [data[2] for data in bifurcation_data]
max_real_eigenvalues = [data[3] for data in bifurcation_data]

# Reshape data for heatmap
mY_grid = reshape(mY_values, (length(mY_range), length(mU_range)))
mU_grid = reshape(mU_values, (length(mY_range), length(mU_range)))
eigenvalue_grid = reshape(max_real_eigenvalues, (length(mY_range), length(mU_range)))

# Plot 2D bifurcation diagram
heatmap(mY_range, mU_range, eigenvalue_grid, xlabel="mY", ylabel="mW", xscale=:log10, yscale=:log10,
        title="Bifurcation Diagram (Max Real Part of Eigenvalues)", color=:viridis)