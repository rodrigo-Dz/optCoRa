using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles;
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
p = Dict([
    :g   => 0.0001,    # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,     # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,         # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 0.01,     # U synthesis rate dependent of Y (nM/min)
    :gU  => 0.0001,    # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mW  => 100,       # W constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gW  => 0.0001,    # W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :e0  => 0.0001,    # U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    :eP  => 0.001,    # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :eM  => 0.5,       # U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)
    :mUs => NaN,       # LOCAL: Ys constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
]);


p_copy = p

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


# Función para realizar el análisis de bifurcación
function bifurcation_analysis(n_points, a, b, n_params, u0, system_FB, system_nFB, tspan)
    results = []
    # Crea el generador de secuencias Sobol en 10 dimensiones
    sobol_p = SobolSeq(n_params)
    # Genera los puntos y aplica la transformación logarítmica
    sobol_p = [10.0 .^ (a .+ (b - a) .* next!(sobol_p)) for _ in 1:n_points]

    for i in 1:n_points
        
        p = copy(sobol_p[i])  # Crear una copia de los parámetros
        copy_p = copy(p)  

        #u_eq = find_equilibrium(params, u0, system!)
        p[1] = 0.0001
        p[3] = 0.1
        p[5] = 0.0001
        p[7] = 0.0001
        p[8] = 0.0001
        p[10] = 0.5
        sol = solve_to_steady_state(system_FB, u0, p, tspan)
        SS = sol.u[end]  # Usar el último punto de la solución como estado estacionario
        FB = (SS[1])  # Almacenar el valor de Y en el estado estacionario

        p[2] = p[2] * 1.05
        p[1] = 0.0001
        p[3] = 0.1
        p[5] = 0.0001
        p[7] = 0.0001
        p[8] = 0.0001
        p[10] = 0.5
        sol = solve_to_steady_state(system_FB, SS, p, tspan)
        ss = sol.u[end]
        FB_p = (ss[1])


        p = copy(copy_p)
        p[1] = 0.0001
        p[1] = 0.0001
        p[3] = 0.1

        p[5] = 0.0001
        p[7] = 0.0001
        p[8] = 0.0001
        p[10] = 0.5

        p[11] = p[4] * FB
        sol = solve_to_steady_state(system_nFB, SS, p, tspan)
        ss = sol.u[end]
        nFB = (ss[1])


        p = copy(copy_p)
        p[1] = 0.0001
        p[1] = 0.0001
        p[3] = 0.1

        p[5] = 0.0001
        p[7] = 0.0001
        p[8] = 0.0001
        p[10] = 0.5
        p[2] = p[2] * 1.05
        p[11] = p[4] * FB

        sol = solve_to_steady_state(system_nFB, SS, p, tspan)
        ss = sol.u[end]
        nFB_p = (ss[1])

        cora = log10(FB_p/FB) / log10.(nFB_p/nFB)
        push!(results,  vcat(copy_p, cora))
    end

    # Convertir results a una matriz
    results_matrix = reduce(vcat, transpose.(results))
    # Guardar los resultados en un archivo CSV
    writedlm("bifurcation_results_3p_7.csv", results_matrix, ',')

    return results_matrix

end

# Condiciones iniciales para el primer sistema
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C

# Rango de tiempo (aumentado para asegurar convergencia al SS)
tspan = (0.0, 1e8)  # Tiempo largo para alcanzar el estado estacionario

#writedlm(outfile1, [vcat(string("Row"), r)],'\t');
a = -3.0  # límite inferior para la escala logarítmica, e.g., 10^1
b = 3.0  # límite superior para la escala logarítmica, e.g., 10^3
n_points = 4096 # cantidad de puntos (recomendable usar una potencia de 2)
n_params = 11


r = bifurcation_analysis(n_points, a, b, n_params, u0, ode_system!, ode_systemNF!, tspan)
println(r)


