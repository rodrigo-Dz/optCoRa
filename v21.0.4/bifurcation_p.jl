using DifferentialEquations
using Plots
using NLsolve

# Define el sistema de ODEs
function logistic_growth!(du, u, p, t)
    x = u[1]
    r, K, h = p
    du[1] = r * x * (1 - x/K) - h
end

# Parámetros
r = 1.0  # Tasa de crecimiento intrínseco
K = 10.0 # Capacidad de carga
h_range = range(0.0, 3.0, length=100)  # Rango de h para la bifurcación

# Condiciones iniciales
u0 = [1.0]  # Población inicial

# Tiempo de simulación
tspan = (0.0, 10.0)

# Función para encontrar puntos de equilibrio
function find_equilibrium(p)
    r, K, h = p
    # Resolver la ecuación cuadrática: r*x*(1 - x/K) - h = 0
    discriminant = r^2 - 4*r*h/K
    if discriminant >= 0
        x1 = (r + sqrt(discriminant)) * K / (2*r)
        x2 = (r - sqrt(discriminant)) * K / (2*r)
        return [x1, x2]
    else
        return []  # No hay equilibrios reales
    end
end

# Bifurcation analysis
bifurcation_data = []

for h in h_range
    p = [r, K, h]
    equilibria = find_equilibrium(p)
    if !isempty(equilibria)
        push!(bifurcation_data, (h, equilibria))
    end
end

# Preparar datos para el diagrama de bifurcación
h_values = [data[1] for data in bifurcation_data]
x1_values = [data[2][1] for data in bifurcation_data]
x2_values = [data[2][2] for data in bifurcation_data]

# Graficar el diagrama de bifurcación
plot(h_values, x1_values, label="Equilibrio Estable", lw=2, color=:blue)
plot!(h_values, x2_values, label="Equilibrio Inestable", lw=2, linestyle=:dash, color=:red)
xlabel!("h (Tasa de extracción)")
ylabel!("Población (x)")
title!("Diagrama de Bifurcación Saddle-Node")
vline!([r*K/4], linestyle=:dash, color=:black, label="Bifurcación (h = rK/4)")