using DifferentialEquations
using Plots
using NLsolve

# Define el sistema de ODEs completo
function ode_system!(du, u, p, t)
    Y, U, W, C = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM = p

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)  # dU/dt
    du[3] = mW - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] = -((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end

# Parámetros del sistema completo
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
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C (sistema completo)

# Rango de tiempo
tspan = (0.0, 500.0)

# Resolver los sistemas de ODEs
prob = ODEProblem(ode_system!, u0, tspan, params)
prob3 = ODEProblem(ode_system!, u0, tspan, p3)
prob4 = ODEProblem(ode_system!, u0, tspan, p4)
prob5 = ODEProblem(ode_system!, u0, tspan, p5)
prob6 = ODEProblem(ode_system!, u0, tspan, p6)
sol = solve(prob, Tsit5())
sol3 = solve(prob3, Tsit5())
sol4 = solve(prob4, Tsit5())
sol5 = solve(prob5, Tsit5())
sol6 = solve(prob6, Tsit5())


# Encontrar punto de equilibrio del sistema completo
function find_equilibrium(p)
    function equilibrium_condition(u)
        du = similar(u)
        ode_system!(du, u, p, 0.0)
        return du
    end
    result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
    return result.zero
end

u_eq = find_equilibrium(params)
u_eq3 = find_equilibrium(p3)
u_eq4 = find_equilibrium(p4)
u_eq5 = find_equilibrium(p5)
u_eq6 = find_equilibrium(p6)


# Graficar el plano fase (Y vs U) para ambos modelos
#plot(sol, idxs=(1, 2), label="Trayectoria (Completo)", lw=2, xlabel="Y", ylabel="U")
plot(sol3, idxs=(1, 2), label="opt1", lw=2)
plot!(sol4, idxs=(1, 2), label="opt2", lw=2)
plot!(sol5, idxs=(1, 2), label="opt3", lw=2)
plot!(sol6, idxs=(1, 2), label="opt4", lw=2)
#scatter!([u_eq[1]], [u_eq[2]], label="Punto de Equilibrio", color=:red, markersize=8)
scatter!([u_eq3[1]], [u_eq3[2]], label="Eq Op 1", markersize=4)
scatter!([u_eq4[1]], [u_eq4[2]], label="Eq Op 2",  markersize=4)
scatter!([u_eq5[1]], [u_eq5[2]], label="Eq Op 3", markersize=4)
scatter!([u_eq6[1]], [u_eq6[2]], label="Eq Op 4", markersize=4)
title!("Plano Fase: Y vs U (Comparación Completo vs. Reducido)")
