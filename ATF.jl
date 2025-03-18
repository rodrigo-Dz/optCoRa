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
    du[2] = (mUs)    - ((g + gU) * U) - (eP * U * W) + ((e0 + gW ) * C)  # dU/dt
    du[3] = mW       - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)  # dW/dt
    du[4] =          - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)  # dC/dt
end

# Parámetros para el primer sistema
params = [
    0.0001,  #1 g    # Dilution rate (e.g. [0.01,0.24] 1/min)
    0.125,   #2 mY  # Y synthesis rate dependent of W (nM/min)
    0.1,     #3 gY    # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    1.14,    #4 mU   U synthesis rate dependent of Y (nM/min)
    0.0001,  #5 gU   U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    43.75,     #6 mW   W constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    0.0001,  #7 gW   W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    0.0001,  #8 e0   U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    38.45,   #9 eP   U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    0.5,     #10 eM   U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)
    NaN   
]

p_copy = params

# Condiciones iniciales para el primer sistema
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Rango de tiempo (aumentado para asegurar convergencia al SS)
tspan = (0.0, 1e8)  # Tiempo largo para alcanzar el estado estacionario

# Función para resolver el sistema de ODEs y verificar el estado estacionario
function solve_to_steady_state(system!, u0, p, tspan)
    prob = ODEProblem(system!, u0, tspan, p)
    sol = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)  # Usar un solver robusto

    # Verificar si el sistema alcanzó el estado estacionario
    du = similar(u0)
    system!(du, sol.u[end], p, sol.t[end])
    if maximum(abs.(du)) < 1e-6
        #println("El sistema alcanzó el estado estacionario.")
        i = 1
    else
        i = 0
        println("El sistema no alcanzó el estado estacionario.")
    end

    return sol
end

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
    function equilibrium_condition(u)
        du = similar(u)
        system!(du, u, p, 0.0)
        return du
    end
    ForwardDiff.jacobian(equilibrium_condition, u_eq)
end

# Bifurcation analysis: Variar mY y mW en escala logarítmica
mY_range = 10 .^ range(-3, 3, length=30)  # Rango logarítmico para mY: 0.001 a 1000
mW_range = 10 .^ range(-3, 3, length=30)  # Rango logarítmico para mW: 0.001 a 1000

# Función para realizar el análisis de bifurcación
function bifurcation_analysis(mY_range, mW_range, params, u0, system!, system_name, grid, pert)
    bifurcation_data = []  # Almacenar datos de oscilaciones
    for i in 1:length(mY_range), j in 1:length(mW_range)
        params[2] = mY_range[i]  # Actualizar mY
        params[10] = mW_range[j]  # Actualizar mW

        if grid != 0  # definir sintesis fija de U:  mUs = mU*yss
            params[11] = params[4] * grid[j, i]
        end
        
        if pert == 1  # perturbacion al parametro mY
            params[2] = params[2] * 1.05
        end 

        ## Resolver el sistema hasta el estado estacionario
        sol = solve_to_steady_state(system!, u0, params, tspan)
        u_eq = sol.u[end]  # Usar el último punto de la solución como estado estacionario
        
        #u_eq = find_equilibrium(params, u0, system!)

        steady_state_Y = (u_eq[1])  # Almacenar el valor de Y en el estado estacionario

        J = compute_jacobian(u_eq, params, system!)
        eigenvalues = eigvals(J)  # Calcular autovalores

        # Verificar si hay autovalores con parte imaginaria no nula
        has_oscilations = any(imag(l) != 0 && real(l) > 0 for l in eigenvalues)
        if has_oscilations == true
            color = 1  # Eigenvalues complejos con parte real positiva
        else
            color = 0  # Todo lo demás
        end

        

        push!(bifurcation_data, (mY_range[i], mW_range[j], steady_state_Y, color))
    end

    # Preparar datos para el heatmap
    mY_values = [data[1] for data in bifurcation_data]
    mW_values = [data[2] for data in bifurcation_data]
    steady_state_Y = [data[3] for data in bifurcation_data]
    colors = [data[4] for data in bifurcation_data]

    # Reshape data for heatmap
    mY_grid = reshape(mY_values, (length(mY_range), length(mW_range)))
    mW_grid = reshape(mW_values, (length(mY_range), length(mW_range)))
    steady_state_Y_grid = reshape(steady_state_Y, (length(mY_range), length(mW_range)))
    color_grid = reshape(colors, (length(mY_range), length(mW_range)))

    #custom_palette = [:grey, :yellow]
    #heatmap(mY_range, mW_range, color_grid, xlabel="mY", ylabel="mW",
    #    title="Diagrama de Bifurcación", xscale=:log10, yscale=:log10,
    #    color=custom_palette, clim=(1, 2), colorbar=false)
    #savefig("./Results_bifurcation/Oscillations_$system_name.png")    
    
    # Gráfica de regiones con oscilaciones
    #heatmap(mY_range, mW_range, oscillation_grid, xlabel="mY", ylabel="mW",
    #    title="Diagrama de Bifurcación", xscale=:log10, yscale=:log10,
    #    color=custom_palette, clim=(1, 2), colorbar=false)
    #savefig("./Results_bifurcation/Oscillations_$system_name.png")

    #heatmap(mY_range, mW_range, log10.(oscillation_grid), xlabel="mY", ylabel="mW",
     #       title="Regiones con oscilaciones ($system_name)", xscale=:log10, yscale=:log10,
      #      color=:coolwarm, colorbar_title="Oscilaciones (1 = Sí, 0 = No)")

    heatmap(mY_range, mW_range, log10.(steady_state_Y_grid), xlabel="mY", ylabel="mW",
            title="Estado estacionario de Y ($system_name)", xscale=:log10, yscale=:log10,
            color=:coolwarm, colorbar_title="SS Y")
    savefig("./Results_bifurcation/SS_Y_$system_name.png")

    return steady_state_Y_grid, color_grid
end

# Realizar el análisis de bifurcación para el primer sistema
ss_fb, oscilations = bifurcation_analysis(mY_range, mW_range, params, u0, ode_system!, "Sistema_1", 0, 0)

params = p_copy
# Realizar el análisis de bifurcación para el segundo sistema
ss_nf, x = bifurcation_analysis(mY_range, mW_range, params, u0, ode_systemNF!, "Sistema_2", ss_fb, 0)

params = p_copy
ss_fb_p, x = bifurcation_analysis(mY_range, mW_range, params, u0, ode_system!, "Sistema_1_p", 0, 1)

same = (ss_fb == ss_nf)
println(same)

params = p_copy
# Realizar el análisis de bifurcación para el segundo sistema
ss_nf_p, x = bifurcation_analysis(mY_range, mW_range, params, u0, ode_systemNF!, "Sistema_2_p", ss_fb, 1)



CoRa = log10.(ss_fb_p./ss_fb) ./ log10.(ss_nf_p./ss_nf)
using Plots

# Graficar con la escala de colores fija en el rango [0,1]
heatmap(mY_range, mW_range, CoRa, xlabel="mY", ylabel="eM",
        title="CoRas", xscale=:log10, yscale=:log10,
        color=:viridis, colorbar_title="Coritas",
        clim=(0,1))  # Fija la escala de colores entre 0 y 1

# Guardar la imagen
savefig("./Results_bifurcation/CoRas.png")

#savefig("./Results_bifurcation/CoRas_2.png")
CoRa[oscilations .== 1] .= 1


heatmap(mY_range, mW_range, CoRa, xlabel="mY", ylabel="mW",
 title="CoRas", xscale=:log10, yscale=:log10,
 color=:viridis, colorbar_title="Coritas",
 clim=(0,1))  # Fija la escala de colores entre 0 y 1
savefig("./Results_bifurcation/Oscillations.png")    
   
