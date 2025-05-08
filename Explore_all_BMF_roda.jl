using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles;
using Statistics

# Define feedback system
# before run, check that the irden is the same as when you print the dictionary
# mW, eP, gU, gY, g, mY, mU, gW, eM, mUs, e0)

function ode_system!(du, u, p, t)
    Y, A, B, C, U, Us = u
    g, mY, gY, mA, mB, mU, eP, e0, bA, bI, mBs = p

    du[1] = (mY * U) - ((g + gY) * Y)
	du[2] =  mA      - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
	du[3] = (mB * Y) - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
	du[4] =          - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
	du[5] =  mU      - (g * U) + (bA * A * Us) - (bI * B * U)
	du[6] =          - (g * Us) - (bA * A * Us) + (bI * B * U)

end

function ode_systemNF!(du, u, p, t)
    Y, A, B, C, U, Us = u
    g, mY, gY, mA, mB, mU, eP, e0, bA, bI, mBs = p

    du[1] = (mY * U) - ((g + gY) * Y)
	du[2] =  mA      - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
	du[3] = (mBs) - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
	du[4] =          - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
	du[5] =  mU      - (g * U) + (bA * A * Us) - (bI * B * U)
	du[6] =          - (g * Us) - (bA * A * Us) + (bI * B * U)
end

using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.01,    
    :mY  => 1.0,     
    :gY  => 0.1,       
    :mA  => 0.0338,
    :mB  => 0.0125,      
    :mU  => 0.1,
    :eP  => 0.0001,
    :e0  => 0.05,
    :bA  => 0.5,
    :bI  => 100,
    :mBs => NaN
)


# Solve 
function solve_to_steady_state(system!, u0, p, tspan)
    p_values = collect(values(p))
    prob = ODEProblem(system!, u0, tspan, p_values)

    sol = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)  # Usar un solver robusto

    # Verificar si el sistema alcanzó el estado estacionario
    du = similar(u0)
    system!(du, sol.u[end], p_values, sol.t[end])
    if maximum(abs.(du)) < 1e-6
        #println("El sistema alcanzó el estado estacionario.")
        return sol
    else
        println(" no ss")
        return 0
    end

    
end


# Encontrar puntos de equilibrio para el primer sistema

function find_equilibrium(p, u0, system!)
    p_values = collect(values(p))
    function equilibrium_condition(u)
        du = similar(u)
        system!(du, u, p_values, 0.0)
        return du
    end
    result = nlsolve(equilibrium_condition, u0, autodiff=:forward)
    return result.zero
end




# Calcular el Jacobiano en un punto de equilibrio
function compute_jacobian(u_eq, p, system!)
    p_values = collect(values(p))
    function equilibrium_condition(u)
        du = similar(u)
        system!(du, u, p_values, 0.0)
        return du
    end
    ForwardDiff.jacobian(equilibrium_condition, u_eq)
end




function Check(ssR, soR, rtol)
    if isnan(ssR) || isnan(soR) || (abs(ssR) - abs(soR)) > 1e-4
        rtol *= 1e-3
        if(rtol < 1e-24)
            println("ERROR: Check NF system (reltol=",rtol,").")
            #println(vcat(pert.p,i,[p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps],mm.outFB(ssR),mm.outNF(soR)))
            if(abs(ssR - soR)/ssR > 0.01)
                ssR, soR = Restart(ssR, soR)
                println("Error too large. SS results excluded!")
            end
        end
        return ssR, soR, rtol, "Insufficient"
    else
        return ssR, soR, rtol, "Sufficient"
    end
end;


# Bifurcation analysis
function bifurcation_analysis(n_points, sobol_p, u0, system_FB, system_nFB, tspan, pOp, p, mY_range)
    results = []
    p_ori = copy(p)

    for i in 1:n_points
        println(i)
        p = copy(p_ori)
        for j in 1:length(pOp)
			#p[mrw.pOp[i]] = 10.0 .^ (rand(Uniform(mrw.pMin[i], mrw.pMax[i])));
            p[pOp[j]] = sobol_p[i][j];

		end

        curve = zeros(length(mY_range))
        SSs = zeros(length(mY_range))
        next_ci = u0
        for j in 1:length(mY_range)
            p[:mY] = mY_range[j]
            copy_p = copy(p)  

            #sol = find_equilibrium(p, next_ci, system_FB)
            #SS = sol

            sol = solve_to_steady_state(system_FB, next_ci, p, tspan)

            if sol == 0
                n_sol = find_equilibrium(p, next_ci, system_FB)
                J = compute_jacobian(n_sol, p, system_FB)
                eigenvalues = eigvals(J)  # Calcular autovalores
                has_oscilations = any(imag(l) != 0 && real(l) > 0 for l in eigenvalues)
                if has_oscilations == true
                    curve[j] = 1
                    SSs[j] = NaN # Eigenvalues complejos con parte real positiva
                    println("oscila")
                else
                    println("No oscila pero pasa algo raro, intenta con otro solver")  # Todo lo demás
                    curve[j] = 2
                    SSs[j] = NaN
                end     
            else
                #SS = sol
                SS = sol.u[end]  # Usar el último punto de la solución como estado estacionario

                FB = (SS[1])  # Almacenar el valor de Y en el estado estacionario
                next_ci = SS

                p = copy(copy_p)
                p[:mBs] = p[:mB] * FB
    
                #sol_ = find_equilibrium(p, SS, system_nFB)
                #ss_= sol_[1]
    
                sol_ = solve_to_steady_state(system_nFB, SS, p, tspan)
                ss_ = sol_.u[end]
    
                nFB = (ss_[1])


                if abs(FB - nFB) < 0.0001

                    p = copy(copy_p)
                    p[:mY] = p[:mY] * 1.05
    
                    #sol_ = find_equilibrium(p, SS, system_FB)
                    #ss_= sol_[1]
    
                    sol_ = solve_to_steady_state(system_FB, SS, p, tspan)
                    ss_ = sol_.u[end]
    
                    FB_p = (ss_[1])
    
                    p = copy(copy_p)
                    p[:mY] = p[:mY] * 1.05
                    p[:mBs] = p[:mB] * FB
    
                    #sol_ = find_equilibrium(p, SS, system_nFB)
                    #ss_= sol_[1]
    
                    sol_ = solve_to_steady_state(system_nFB, SS, p, tspan)
                    ss_ = sol_.u[end]
    
                    nFB_p = (ss_[1])
    
                    curve[j] = log10(FB_p/FB) / log10(nFB_p/nFB)
                    SSs[j] = FB
                else
                    println("no son iguales y no oscila")
                    curve[j] = 3
                    SSs[j] = NaN
                end  
            end
            
            



          
        end  
        r_01 = count(x -> x < 0.1, curve)
        r_02 = count(x -> x < 0.2, curve)
        r_03 = count(x -> x < 0.3, curve)
        os = count(x -> x == 1, curve)
        raro = count(x -> x == 2, curve)
        noeq = count(x -> x == 3, curve)

        indices = curve .< 0.1  # Condición para filtrar valores en 'a' menores a 0.1
        # Obtener los valores filtrados de 'b'
        filtered_SSs = SSs[indices]
        # Calcular el promedio de los valores filtrados de 'b'
        ss = mean(filtered_SSs)
        push!(results,  vcat(collect(values(p)), r_01, r_02, r_03, os, ss, raro, noeq))
    end

    # Convertir results a una matriz
    results_matrix = reduce(vcat, transpose.(results))
    # Guardar los resultados en un archivo CSV
    writedlm("bifurcation_results_BMF_roda.csv", results_matrix, ',')

    return results_matrix

end

# Condiciones iniciales para el primer sistema
u0 = [0.0, 0.0, 0.0, 0.0]  # Y, U, W, C
mY_range = 10 .^ range(-3, 3, length=25)  # Rango logarítmico para mY: 0.001 a 1000
# Rango de tiempo (aumentado para asegurar convergencia al SS)
tspan = (0.0, 1e8)  # Tiempo largo para alcanzar el estado estacionario
pOp  = [:mU,:mW,:eP]	# Parameters to optimize
#pMin = [-3,-3,-3],		# Minimum parameter value to explore (log10)
#pMax = [3,3,3]
#writedlm(outfile1, [vcat(string("Row"), r)],'\t');
a = -3.0  # límite inferior para la escala logarítmica, e.g., 10^1
b = 3.0  # límite superior para la escala logarítmica, e.g., 10^3
n_points = 2048 # cantidad de puntos (recomendable usar una potencia de 2)
n_params = 3

# Set of parameters    
sobol_p = SobolSeq(n_params)
sobol_p = [10.0 .^ (a .+ (b - a) .* next!(sobol_p)) for _ in 1:n_points]

r = bifurcation_analysis(n_points, sobol_p, u0, ode_system!, ode_systemNF!, tspan, pOp, p, mY_range)

#=
using Base.Threads

# Dividir la tabla de Sobol en dos partes
n_half = div(n_points, 2)
sobol_p1 = sobol_p[1:n_half]
sobol_p2 = sobol_p[n_half+1:end]

# Ejecutar en paralelo
task1 = Threads.@spawn bifurcation_analysis(n_half, sobol_p1, u0, ode_system!, ode_systemNF!, tspan, pOp, p, mY_range)
task2 = Threads.@spawn bifurcation_analysis(n_half, sobol_p2, u0, ode_system!, ode_systemNF!, tspan, pOp, p, mY_range)

# Esperar los resultados
r1 = fetch(task1)
r2 = fetch(task2)

# Combinar los resultados
r = vcat(r1, r2)

writedlm("bifurcation_results_FAD_roda.csv", results_matrix, ',')

=#