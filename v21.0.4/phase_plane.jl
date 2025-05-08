using DifferentialEquations
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
    1.0,   # mY
    0.1,   # gY
    0.0152, # mU
    0.05,  # gU
    0.10,  # mW
    0.0001, # gW
    0.0001, # e0
    0.0375, # eP
    0.5    # eM
]

# Condiciones iniciales
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Rango de tiempo
tspan = (0.0, 100.0)  # Tiempo de simulación

# Resolver el sistema de ODEs
prob = ODEProblem(ode_system!, u0, tspan, params)
sol = DifferentialEquations.solve(prob, Tsit5())  # Calificar solve con DifferentialEquations

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

# Calcular el punto de equilibrio
u_eq = find_equilibrium(params)

# Graficar el plano fase (Y vs U)
plot(sol, idxs=(1, 2), xlabel="Y", ylabel="U", label="Trayectoria", lw=2)
scatter!([u_eq[1]], [u_eq[2]], label="Punto de Equilibrio", color=:red, markersize=8)
title!("Plano Fase: Y vs U")