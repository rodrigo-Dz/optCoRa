using DifferentialEquations
using Plots

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
    0.01,  # g
    1,   # mY
    0.1,   # gY
    1, # mU
    0.05,  # gU
    100,  # mW
    0.0001, # gW
    0.0001, # e0
    100, #0.0375, # eP
    0.5    # eM
]

# Condiciones iniciales
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Rango de tiempo
tspan = (0.0, 500.0)  # Simular desde t=0 hasta t=100

# Crear el problema de ODE
prob = ODEProblem(ode_system!, u0, tspan, params)

# Resolver el sistema de ODEs
sol = DifferentialEquations.solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Graficar la dinámica del sistema
plot(sol, xlabel="Tiempo (t)", ylabel="Concentración", label=["Y" "U" "W" "C"],
     title="Dinámica del Sistema", lw=2)