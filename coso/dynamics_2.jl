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

# Segunda ecuación diferencial para Y sola
function n!(du, u, p, t)
    Y = u[1]
    g, mY, gY = p
    du[1] = mY - ((g + gY) * Y)  # dY/dt
end

# Parámetros
params = [0.0001, 
0.125, #mY
1,  #gY
0.125, 
0.0001, 
0.1, 
0.0001, 
0.0001, 
0.0375, 
0.5]

p2 = [0.0001, 0.125, 1]

p3 = [0.0001, 
0.125, #mY
1,  #gY
5.6526, #mU
0.0001, 
0.7415, #mW 
0.0001, 
0.0001, 
9.1512, #eP 
0.5]

p4 = [0.0001, 
0.125, #mY
1,  #gY
1.2729, #mU
0.0001, 
0.7765, #mW 
0.0001, 
0.0001, 
0.1024, #eP 
0.5]

p5 = [0.0001, 
0.125, #mY
1,  #gY
3.9436, #mU
0.0001, 
0.04511, #mW 
0.0001, 
0.0001, 
117.9416, #eP 
0.5]

p6 = [0.0001, 
0.125, #mY
1,  #gY
7.1322, #mU
0.0001, 
0.5231, #mW 
0.0001, 
0.0001, 
18.7132, #eP 
0.5]
# Condiciones iniciales
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C
u20 = [0.0]  # Solo Y

# Rango de tiempo
tspan = (0.0, 500.0)

# Crear y resolver los problemas de ODE
prob = ODEProblem(ode_system!, u0, tspan, params)
prob2 = ODEProblem(n!, u20, tspan, p2)
prob3 = ODEProblem(ode_system!, u0, tspan, p3)
prob4 = ODEProblem(ode_system!, u0, tspan, p4)
prob5 = ODEProblem(ode_system!, u0, tspan, p5)
prob6 = ODEProblem(ode_system!, u0, tspan, p6)

sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
sol2 = solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8)
sol3 = solve(prob3, Tsit5(), reltol=1e-8, abstol=1e-8)
sol4 = solve(prob4, Tsit5(), reltol=1e-8, abstol=1e-8)
sol5 = solve(prob5, Tsit5(), reltol=1e-8, abstol=1e-8)
sol6 = solve(prob6, Tsit5(), reltol=1e-8, abstol=1e-8)

# Graficar los resultados de Y en ambas soluciones
plot(sol.t, sol[1, :], label="start", xlabel="Tiempo", ylabel="Concentración", lw=2)
plot!(sol2.t, sol2[1, :], label="Modelo Reducido", lw=2, linestyle=:dash)
plot!(sol3.t, sol3[1, :], label="op1", lw=2)
plot!(sol4.t, sol4[1, :], label="op2", lw=2)
plot!(sol5.t, sol5[1, :], label="op3", lw=2)
plot!(sol6.t, sol6[1, :], label="op4", lw=2)
title!("Comparación de Y en ambos modelos")
