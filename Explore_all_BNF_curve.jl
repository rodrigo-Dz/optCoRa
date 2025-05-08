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
    Y, U, Us= u
    g, mY, gY, mU, kD, gU, b, bs, mUs = p

    du[1] = (mY * U)             - ((g + gY) * Y)
    du[2]  = (mU * (kD/(Y + kD))) - ((g + gU) * U)  - (b * U) + (bs * Us)
    du[3] =                      - ((g + gU) * Us) + (b * U) - (bs * Us)
end

# Define anologous system
function ode_systemNF!(du, u, p, t)
    Y, U, Us= u

    g, mY, gY, mU, kD, gU, b, bs, mUs = p

    du[1] = (mY * U)             - ((g + gY) * Y)
    du[2]  = (mUs) - ((g + gU) * U)  - (b * U) + (bs * Us)
    du[3] =                      - ((g + gU) * Us) + (b * U) - (bs * Us)
end

using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 1.0,       # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 31.622777,       # U maximum synthesis rate repressed by Y (nM/min)
    :kD  => 1.0,       # Y repression KD (nM)
    :gU  => 0.0001,    # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :b   => 0.02,    # U to Us transition rate (1/min)
    :bs  => 0.0316,    # Us to U transition rate (1/min)
    :mUs => NaN,  
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
        return 0.0
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
function bifurcation_analysis(n_points, a, b, n_params, u0, system_FB, system_nFB, tspan, pOp, p, mY_range)
    results = []
    curve = zeros(length(mY_range))
    SSs = zeros(length(mY_range))
    next_ci = u0
    for j in 1:length(mY_range)
        p[:mY] = mY_range[j]
        copy_p = copy(p)  
        sol = find_equilibrium(p, next_ci, system_FB)
        #sol = solve_to_steady_state(system_FB, next_ci, p, tspan)
        SS = sol
        #SS = sol.u[end]  # Usar el último punto de la solución como estado estacionario
        FB = (SS[1])  # Almacenar el valor de Y en el estado estacionario
        next_ci = SS

        p = copy(copy_p)
        p[:mUs] = p[:mU] * (p[:kD]/(FB + p[:kD]))
        sol_ = find_equilibrium(p, SS, system_nFB)
        ss_= sol_[1]
        #sol = solve_to_steady_state(system_nFB, SS, p, tspan)
        #ss_ = sol.u[end]
        nFB = (ss_[1])

        if abs(FB - nFB) < 0.0001

            p = copy(copy_p)
            p[:mY] = p[:mY] * 1.05
            sol_ = find_equilibrium(p, SS, system_FB)
            ss_= sol_[1]
            #sol = solve_to_steady_state(system_FB, SS, p, tspan)
            #ss_ = sol.u[end]
            FB_p = (ss_[1])

            p = copy(copy_p)
            p[:mY] = p[:mY] * 1.05
            p[:mUs] = p[:mU] * (p[:kD]/(FB + p[:kD]))
            sol_ = find_equilibrium(p, SS, system_nFB)
            ss_= sol_[1]
            #sol = solve_to_steady_state(system_nFB, SS, p, tspan)
            #ss_ = sol.u[end]
            nFB_p = (ss_[1])

            if (FB_p < 0) || (FB<0) || (nFB_p<0) || (nFB<0)
                println("error en uno de estos:", FB, FB_p, nFB_p, nFB)
                curve[j] = NaN
                SSs[j] = NaN
            else
                curve[j] = log10(FB_p/FB) / log10(nFB_p/nFB)
                SSs[j] = FB
            end
        else
            J = compute_jacobian(sol, p, system_FB)
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
            println("no son iguales")
            #curve[j] = NaN
            #SSs[j] = NaN
        end            
    end  
    return curve

end

# Condiciones iniciales para el primer sistema
u0 = [0.0, 0.0, 0.0]  # Y, U, W, C
mY_range = 10 .^ range(-3, 3, length=25)  # Rango logarítmico para mY: 0.001 a 1000
# Rango de tiempo (aumentado para asegurar convergencia al SS)
tspan = (0.0, 1e8)  # Tiempo largo para alcanzar el estado estacionario
pOp  =  [:mU,:b,:bs]	# Parameters to optimize
#pMin = [-3,-3,-3],		# Minimum parameter value to explore (log10)
#pMax = [3,3,3]
#writedlm(outfile1, [vcat(string("Row"), r)],'\t');
a = -3.0  # límite inferior para la escala logarítmica, e.g., 10^1
b = 3.0  # límite superior para la escala logarítmica, e.g., 10^3
n_points = 2 # cantidad de puntos (recomendable usar una potencia de 2)
n_params = 3


r = bifurcation_analysis(n_points, a, b, n_params, u0, ode_system!, ode_systemNF!, tspan, pOp, p, mY_range)
println(r)

r_01 = count(x -> x < 0.1, r) 
println(r_01)
plot(mY_range, r, seriestype = :line, xscale = :log10,
     xlabel = "mY", ylabel = "CoRa(mY)", ylims=[0,1],
     legend = :topright, grid = true)