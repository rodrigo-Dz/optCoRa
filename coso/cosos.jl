using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve

# Define el primer sistema de ODEs
function ode_system!(du, u, p, t)
    Y, U, W, C = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM, mUs = p

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mU * Y) - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
    du[3] = mW       - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] =          - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end


function ode_systemNF!(du, u, p, t)
    Y, U, W, C = u
    g, mY, gY, mU, gU, mW, gW, e0, eP, eM, mUs = p

    du[1] = (mY * W) - ((g + gY) * Y)  # dY/dt
    du[2] = (mUs)    - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)  # dU/dt
    du[3] = mW       - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] =          - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end





# Parámetros para el primer sistema
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
    NaN
]

p_copy = params

# Condiciones iniciales para el primer sistema
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C


# Encontrar puntos de equilibrio para el primer sistema
function find_equilibrium(p, u0, system!)
    function equilibrium_condition(u)
        du = similar(u)
        system!(du, u, p, 0.0)
        return du
    end
    result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
    return result.zero
end

# Calcular el Jacobiano en un punto de equilibrio
function compute_jacobian(u_eq, p, system!)
    function system_derivatives(u)
        du = similar(u)
        system!(du, u, p, 0.0)
        return du
    end
    ForwardDiff.jacobian(system_derivatives, u_eq)
end

# Bifurcation analysis: Variar mY y mW en escala logarítmica
mY_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mY: 0.001 a 1000
mW_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mW: 0.001 a 1000

# Función para realizar el análisis de bifurcación
function bifurcation_analysis(mY_range, mW_range, params, u0, system!, system_name, grid)
    oscillation_data = []  # Almacenar datos de oscilaciones
    for i in 1:length(mY_range), j in 1:length(mW_range)
        params[2] = mY_range[i]  # Actualizar mY
        params[6] = mW_range[j]  # Actualizar mW

        if grid != 0
            params[11] = 0.125 * grid[j, i]
            
        end

        #if pert != 0
        #    params[pert] = mY * 2
        #end 

        u_eq = find_equilibrium(params, u0, system!)
        steady_state_Y = (u_eq[1])  # Almacenar el valor de Y en el estado estacionario

        J = compute_jacobian(u_eq, params, system!)
        eigenvalues = eigvals(J)  # Calcular autovalores

        # Verificar si hay autovalores con parte imaginaria no nula
        has_oscillations = any(imag(l) != 0 && real(l) > 0 for l in eigenvalues)
        
        push!(oscillation_data, (mY_range, mW_range, steady_state_Y, has_oscillations))
    end

    # Preparar datos para el heatmap
    mY_values = [data[1] for data in oscillation_data]
    mW_values = [data[2] for data in oscillation_data]
    steady_state_Y = [data[3] for data in oscillation_data]
    has_oscillations = [data[4] for data in oscillation_data]

    println(steady_state_Y)
    # Reshape data for heatmap
    mY_grid = reshape(mY_values, (length(mY_range), length(mW_range)))
    mW_grid = reshape(mW_values, (length(mY_range), length(mW_range)))
    steady_state_Y_grid = reshape(steady_state_Y, (length(mY_range), length(mW_range)))
    oscillation_grid = reshape(has_oscillations, (length(mY_range), length(mW_range)))
    println(steady_state_Y_grid)
    # Gráfica de regiones con oscilaciones
    heatmap(mY_range, mW_range, log10.(oscillation_grid), xlabel="mY", ylabel="mW",
            title="Regiones con oscilaciones ($system_name)", xscale=:log10, yscale=:log10,
            color=:coolwarm, colorbar_title="Oscilaciones (1 = Sí, 0 = No)")
    savefig("./Results_bifurcation/Oscillations_$system_name.png")

    heatmap(mY_range, mW_range, log10.(steady_state_Y_grid), xlabel="mY", ylabel="mW",
            title="Estado estacionario de Y ($system_name)", xscale=:log10, yscale=:log10,
            color=:coolwarm, colorbar_title="SS Y")
    savefig("./Results_bifurcation/SS_Y_$system_name.png")

    return steady_state_Y_grid
end

# Realizar el análisis de bifurcación para el primer sistema
ss_fb = bifurcation_analysis(mY_range, mW_range, params, u0, ode_system!, "Sistema_1", 0)

params = p_copy
# Realizar el análisis de bifurcación para el segundo sistema
bifurcation_analysis(mY_range, mW_range, params, u0, ode_systemNF!, "Sistema_2", ss_fb)


#bifurcation_analysis(mY_range, mW_range, params, u0, ode_system!, "Sistema_1_pert", 2, 1)

# Realizar el análisis de bifurcación para el segundo sistema
#bifurcation_analysis(mY_range, mW_range, params2, u02, ode_system2!, "Sistema_2_pert", 11, 5)