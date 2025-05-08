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
    1, # mU (se variará)
    0.05,  # gU
    400,  # mW
    0.0001, # gW
    0.0001, # e0
    100, # eP
    0.5  # eM
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

    # Extraer partes reales e imaginarias de todos los eigenvalues
    real_parts = real.(eigenvalues)
    imaginary_parts = imag.(eigenvalues)

    push!(bifurcation_data, (mY, mU, real_parts, imaginary_parts))
end

# Preparar datos para los heatmaps
mY_values = [data[1] for data in bifurcation_data]
mU_values = [data[2] for data in bifurcation_data]
real_parts_all = vcat([data[3] for data in bifurcation_data]...)  # Aplanar la lista
imaginary_parts_all = vcat([data[4] for data in bifurcation_data]...)  # Aplanar la lista

# Normalizar los valores de color
normalize_color(vals) = (vals .- minimum(vals)) ./ (maximum(vals) - minimum(vals))

real_parts_normalized = normalize_color(real_parts_all)
imaginary_parts_normalized = normalize_color(imaginary_parts_all)

# Crear una paleta de colores
palette = cgrad(:viridis)  # Puedes cambiar :viridis por otra paleta

# Gráfica 1: Parte real de los eigenvalues
scatter(
    repeat(mY_values, inner=4),  # Repetir mY para cada eigenvalue
    repeat(mU_values, inner=4),  # Repetir mU para cada eigenvalue
    color=[palette[val] for val in real_parts_normalized],  # Mapear valores a la paleta
    xlabel="mY", ylabel="mU",
    title="Parte Real de Eigenvalues", xscale=:log10, yscale=:log10,
    markersize=4, markerstrokewidth=0
)
savefig("p_real.png")


# Gráfica 2: Parte imaginaria de los eigenvalues
scatter(
    repeat(mY_values, inner=4),  # Repetir mY para cada eigenvalue
    repeat(mU_values, inner=4),  # Repetir mU para cada eigenvalue
    color=[palette[val] for val in imaginary_parts_normalized],  # Mapear valores a la paleta
    xlabel="mY", ylabel="mU",
    title="Parte Imaginaria de Eigenvalues", xscale=:log10, yscale=:log10,
    markersize=4, markerstrokewidth=0
)
savefig("p_imag.png")