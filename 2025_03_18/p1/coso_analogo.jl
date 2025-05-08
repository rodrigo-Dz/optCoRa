
using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve

# Define el sistema de ODEs

function ode_system!(du, u, p, t)
    Y, U, W, C, Ys = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM, mYs = p

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
    du[3] = mW - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] = -((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
    du[5] = (mYs * W) - ((g + gY) * Ys)  # dY/dt

end

# Parámetros
params = [
        0.0001,  # g
        0.125,   # mY (se variará)
        1,       # gY
        0.125,   # mU 
        0.0001,  # gU
        0.1,     # mW
        0.0001,  # gW
        0.0001,  # e0
        0.0375,  # eP
        0.5,      # eM
        0.01     #mYp
]

# Condiciones iniciales
u0 = [0.0, 0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Rango de tiempo
tspan = (0.0, 10000.0)  # Tiempo largo para alcanzar el estado estacionario

# Encontrar puntos de equilibrio
function find_equilibrium(p)
    function equilibrium_condition(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
    return result.zero
end

# Calcular el Jacobiano en un punto de equilibrio
function compute_jacobian(u_eq, p)
    function system_derivatives(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    ForwardDiff.jacobian(system_derivatives, u_eq)
end

# Bifurcation analysis: Variar mY y mW en escala logarítmica
mY_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mY: 0.001 a 1000
mW_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mW: 0.001 a 1000
oscillation_data = []  # Almacenar datos de oscilaciones

for mY in mY_range, mW in mW_range
    params[2] = mY  # Actualizar mY
    params[6] = mW  # Actualizar mW
    params[11] = params[2] * 1.5
    u_eq = find_equilibrium(params)

    steady_state_Y = log10(u_eq[5])  # Almacenar el valor de Y en el estado estacionario

    J = compute_jacobian(u_eq, params)
    eigenvalues = eigvals(J)  # Calcular autovalores

    # Verificar si hay autovalores con parte imaginaria no nula
    has_oscillations = any(abs.(imag.(eigenvalues)) .> 1e-6)

    push!(oscillation_data, (mY, mW, steady_state_Y, has_oscillations))
end

# Preparar datos para el heatmap
mY_values = [data[1] for data in oscillation_data]
mW_values = [data[2] for data in oscillation_data]
steady_state_Y = [data[3] for data in oscillation_data]
has_oscillations = [data[4] for data in oscillation_data]

# Reshape data for heatmap
mY_grid = reshape(mY_values, (length(mY_range), length(mW_range)))
mW_grid = reshape(mW_values, (length(mY_range), length(mW_range)))
steady_state_Y_grid = reshape(steady_state_Y, (length(mY_range), length(mW_range)))
oscillation_grid = reshape(has_oscillations, (length(mY_range), length(mW_range)))

# Gráfica de regiones con oscilaciones
heatmap(mY_range, mW_range, oscillation_grid, xlabel="mY", ylabel="mW",
        title="Regiones con oscilaciones", xscale=:log10, yscale=:log10,
        color=:coolwarm, colorbar_title="Oscilaciones (1 = Sí, 0 = No)")
savefig("./Results_bifurcation/Oscillations_an.png")

heatmap(mY_range, mW_range, steady_state_Y_grid, xlabel="mY", ylabel="mW",
        title="Regiones con oscilaciones", xscale=:log10, yscale=:log10,
        color=:coolwarm, colorbar_title="SS Y")
savefig("./Results_bifurcation/SS_an.png")