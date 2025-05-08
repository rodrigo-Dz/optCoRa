using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve

# Define el sistema de ODEs
function ode_system!(du, u, p, t)
    Y, U, W, C = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM = p

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)  # dU/dt
    du[3] = mW - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] = -((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end

# Parámetros
params = [
    0.01,  # g
    1.0,   # mY (se variará)
    0.1,   # gY
    0.001, # mU (se variará)
    0.1,  # gU
    400,  # mW
    0.0001, # gW
    0.001, # e0
    1, # eP
    100    # eM
]

# Condiciones iniciales
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Rango de tiempo
tspan = (0.0, 1000.0)  # Tiempo largo para alcanzar el estado estacionario

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
    # Función que retorna las derivadas como un arreglo
    function system_derivatives(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    ForwardDiff.jacobian(system_derivatives, u_eq)
end

# Bifurcation analysis: Variar mY y mU en escala logarítmica
mY_range = 10 .^ range(-3, 3, length=100)  # Rango logarítmico para mY: 0.01 a 100
mU_range = 10 .^ range(-3, 3, length=100)  # Rango logarítmico para mU: 0.01 a 100
bifurcation_data = []

for mY in mY_range, mU in mU_range
    params[2] = mY  # Actualizar mY
    params[6] = mU  # Actualizar mU
    u_eq = find_equilibrium(params)
    J = compute_jacobian(u_eq, params)
    eigenvalues = eigen(J).values

    # Determinar el color basado en las condiciones
    has_complex_positive = any(imag(λ) != 0 && real(λ) > 0 for λ in eigenvalues)
    if has_complex_positive
        color = :orange  # Eigenvalues complejos con parte real positiva
    else
        color = :grey  # Todo lo demás
    end

    push!(bifurcation_data, (mY, mU, color))
end

# Preparar datos para el heatmap
mY_values = [data[1] for data in bifurcation_data]
mU_values = [data[2] for data in bifurcation_data]
colors = [data[3] for data in bifurcation_data]

# Reshape data for heatmap
mY_grid = reshape(mY_values, (length(mY_range), length(mU_range)))
mU_grid = reshape(mU_values, (length(mY_range), length(mU_range)))
color_grid = reshape(colors, (length(mY_range), length(mU_range)))

# Crear una paleta de colores personalizada
custom_palette = [:grey, :orange]

# Graficar el heatmap
heatmap(mY_range, mU_range, color_grid, xlabel="mY", ylabel="mW",
        title="Diagrama de Bifurcación", xscale=:log10, yscale=:log10,
        color=custom_palette, clim=(1, 2), colorbar=false)