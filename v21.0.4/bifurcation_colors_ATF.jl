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
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
    du[3] = mW - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] = -((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end

# Parámetros
params = [
    0.0001,  # g
    1.0,   # mY (se variará)
    0.1,   # gY
    0.01, # mU (se variará)
    0.0001,  # gU
    4,  # mW
    0.0001, # gW
    0.0001, # e0
    0.001, # eP
    0.5    # eM
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
mY_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mY: 0.01 a 100
mU_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mU: 0.01 a 100
bifurcation_data = []

for mY in mY_range, mU in mU_range
    params[2] = mY  # Actualizar mY
    params[6] = mU  # Actualizar mU
    u_eq = find_equilibrium(params)
    J = compute_jacobian(u_eq, params)
    eigenvalues = eigen(J).values

    # Parte real positiva máxima (usando `init` para evitar errores con colecciones vacías)
    #max_real_positive = maximum(real(λ) for λ in eigenvalues if real(λ) > 0; init=0.0)
    max_real_positive = log10(abs(maximum(real(λ) for λ in eigenvalues)))
    min_real = log10(abs(minimum(real(λ) for λ in eigenvalues)))
    # Parte imaginaria máxima (usando `init` para evitar errores con colecciones vacías)
    max_imaginary = maximum(abs(imag(λ)) for λ in eigenvalues; init=0.0)

    push!(bifurcation_data, (mY, mU, max_real_positive, min_real, max_imaginary))
end

# Preparar datos para los heatmaps
mY_values = [data[1] for data in bifurcation_data]
mU_values = [data[2] for data in bifurcation_data]
max_real_positive_values = [data[3] for data in bifurcation_data]
min_real_values = [data[4] for data in bifurcation_data]
max_imaginary_values = [data[5] for data in bifurcation_data]

# Reshape data for heatmaps
mY_grid = reshape(mY_values, (length(mY_range), length(mU_range)))
mU_grid = reshape(mU_values, (length(mY_range), length(mU_range)))
max_real_positive_grid = reshape(max_real_positive_values, (length(mY_range), length(mU_range)))
min_real_grid = reshape(min_real_values, (length(mY_range), length(mU_range)))
max_imaginary_grid = reshape(max_imaginary_values, (length(mY_range), length(mU_range)))

# Gráfica 1: Parte real positiva de los eigenvalues
heatmap(mY_range, mU_range, max_real_positive_grid, xlabel="mY", ylabel="mW",
        title="Parte Real Maxima de Eigenvalues", xscale=:log10, yscale=:log10,
        color=:heat, colorbar_title="Máx Re(λ) > 0")
        savefig("Bif_max_real.png")

heatmap(mY_range, mU_range, min_real_grid, xlabel="mY", ylabel="mW",
    title="Parte Real Minima de Eigenvalues", xscale=:log10, yscale=:log10,
    color=:viridis, colorbar_title="Min Re(λ) > 0")
    savefig("Bif_min_real.png")
# Gráfica 2: Parte imaginaria de los eigenvalues
heatmap(mY_range, mU_range, max_imaginary_grid, xlabel="mY", ylabel="mW",
        title="Parte Imaginaria de Eigenvalues", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Máx |Im(λ)|")
        savefig("Bif_imag.png")