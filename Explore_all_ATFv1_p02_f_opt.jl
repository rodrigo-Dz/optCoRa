using DifferentialEquations
using ForwardDiff
using LinearAlgebra
using Plots
using NLsolve
using Sobol
using DelimitedFiles;
using Statistics
using Distributions
# Define feedback system
# before run, check that the irden is the same as when you print the dictionary
# mW, eP, gU, gY, g, mY, mU, gW, eM, mUs, e0)

function ode_system!(du, u, p, t)
    Y, U, W, C , Y0, Y1= u
    g, mY, kD, gY, mU, gU, mW, gW, e0, eP, eM, m0, m1, k1, mUs = p

    du[1] = (mY * (kD/(Y1 + kD))) - ((g + gY) * Y) # dY/dt
    du[2] = (mU * Y)              - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C) # dU/dt
    du[3] =    mW                 - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C) # dW/dt
    du[4] =                       - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C) # dC/dt
    du[5] = (m0 * W)              - ((g + gY) * Y0)
    du[6] = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)
end

# Define anologous system
function ode_systemNF!(du, u, p, t)
    Y, U, W, C , Y0, Y1= u
    g, mY, kD, gY, mU, gU, mW, gW, e0, eP, eM, m0, m1, k1, mUs = p

    du[1] = (mY * (kD/(Y1 + kD))) - ((g + gY) * Y) # dY/dt
    du[2] =   (mUs)   - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
    du[3] =    mW                 - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C) # dW/dt
    du[4] =                       - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C) # dC/dt
    du[5] = (m0 * W)              - ((g + gY) * Y0)
    du[6] = (m1 * (k1/(Y0 + k1))) - ((g + gY) * Y1)
end

using DataStructures

# Diccionario ordenado
p = OrderedDict(
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,     # Y maximum synthesis rate dependent of Y1 (nM/min)
	:kD  => 1.0,       # Activation threshold for Y synthesis dependent of Y1 (nM)
    :gY  => 1.0,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 0.01,     # U synthesis rate dependent of Y (nM/min)
    :gU  => 0.0001,    # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mW  => 0.01,       # W constitutive synthesis rate (e.g. [0.005,0.02] nM/min)
    :gW  => 0.0001,    # W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :e0  => 0.0001,    # U:W dissociation (unbinding) rate (e.g. [0.05,140] 1/min)
    :eP  => 0.01,    # U:W association (binding) rate (e.g. [0.0012,2000] 1/nM 1/min)
    :eM  => 0.5,       # U:W (C) mutual anhilation/co-degradation (e.g. [0.03,2] 1/nM 1/min)
    :m0  => 1.25,      # Y0 synthesis rate dependent of W (nM/min)
    :m1  => 12.5,      # Y1 maximum synthesis rate dependent of Y0 (nM/min)
	:k1  => 1.0,       # Activation threshold for Y1 synthesis dependent of Y0 (nM)
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


function CC(mY_range, u0, p, system_FB, system_nFB)
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
        p[:mUs] = p[:mU] * FB
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
            p[:mUs] = p[:mU] * FB
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



function mrwPar(pOp, cov, M, pMin, pMax, p)
    r0 = zeros(length(pOp));			# Vector of parameters to optimize in reference system
    rI = rand(MvNormal(zeros(length(pOp)), zeros(length(pOp)) .+ cov)); # Random values to update parameters
    for pI in 1:length(pOp)
        r0[pI] = p[pOp[pI]];			 # Save previous value
        p[pOp[pI]] *= (M .^ rI[pI]); # Update value
        p[:mUs] = NaN
    #p[mrw.pOp[pI]] *= (M .^ rI[pI]); # Update value
    # Exclude values outside regime of exploration:
        if p[pOp[pI]] < (10.0 ^ pMin[pI])
            p[pOp[pI]] = (10.0 ^ pMin[pI])
        elseif p[pOp[pI]] > (10.0 ^ pMax[pI])
            p[pOp[pI]] = (10.0 ^ pMax[pI])
        end
    end 
    p1 = copy(p)
    return p1
end

# Bifurcation analysis
function optimization(p, mY_range, u0, iter, pOp, cov, M, pMin, pMax, system_FB, system_nFB)
    results = []
    p_ori = copy(p)

    # initialize
    s0 = CC(mY_range, u0, p, system_FB, system_nFB)
    m0 = count(x -> x < 0.1, s0)
    min0 = minimum(s0)
    p0 = copy(p)
    println(p0)

    i = 1
    while (i != iter)
        p1 = copy(p0)
        r0 = zeros(length(pOp));			# Vector of parameters to optimize in reference system
        rI = rand(MvNormal(zeros(length(pOp)), zeros(length(pOp)) .+ cov)); # Random values to update parameters
        for pI in 1:length(pOp)
            r0[pI] = p1[pOp[pI]];			 # Save previous value
            p1[pOp[pI]] *= (M .^ rI[pI]); # Update value
            p1[:mUs] = NaN
        #p[mrw.pOp[pI]] *= (M .^ rI[pI]); # Update value
        # Exclude values outside regime of exploration:
            if p1[pOp[pI]] < (10.0 ^ pMin[pI])
                p1[pOp[pI]] = (10.0 ^ pMin[pI])
            elseif p1[pOp[pI]] > (10.0 ^ pMax[pI])
                p1[pOp[pI]] = (10.0 ^ pMax[pI])
            end
        end 

        i = i + 1
       
        s1 = CC(mY_range, u0, p1, system_FB, system_nFB)
        m1 = count(x -> x < 0.1, s1)
        os = count(x -> x == 1, s1)
        min1 = minimum(s1)
        println(m0, ", ", m1)
        println(min0,", ", min1)
        println(p0[:mU], ", ", p1[:mU])


		## If CoRa>=eps for some conditions, evaluate the |CoRa<=eps| for both sets:
		### NOTE: As op0,op1=[0,1], correct exponential with the expected variance of ~U(0,1)
#=
        xiC = (min0 ^ 2) / (2 * 0.083);
		xiP = (min1 ^ 2) / (2 * 0.083);
		c2 = (minimum([m0,m1])==1.0) && (rand() < exp(xiC - xiP));
        
        xiC = ((1 - m0)^2) / (2 * 0.083);
        xiP = ((1 - m1)^2) / (2 * 0.083);
        c3 = rand() < exp(xiC - xiP);

        if c3 || c2
            s0 = s1
            m0 = m1
            min0 = min1
            p0 = copy(p1)
        else
                # If not, revert to previous parameter values
            p0 = p0
        end  
=#
		# Evaluate if accept new parameter values or not:
		## Only accept in the regime of interest, i.e. CoRa>=0:
		## If CoRa>eps for all conditions, evaluate the min(CoRa) for both sets:
		### NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)

		## If CoRa>=eps for some conditions, evaluate the |CoRa<=eps| for both sets:
		### NOTE: As op0,op1=[0,1], correct exponential with the expected variance of ~U(0,1)
 
        if min1 < min0 && min0 > 0.1 && os == 0.0 && min1>0.0
           s0 = s1
            m0 = m1
            min0 = min1
            p0 = copy(p1)
        end
        if m1 >= m0  && os == 0.0 && min0 > 0.0 && min1 > 0.0
                # If yes, update "reference" system
            s0 = s1
            m0 = m1
            min0 = min1
            p0 = copy(p1)
        else
                # If not, revert to previous parameter values
            p0 = p0
        end  
       # println(p1)
        push!(results,  vcat(collect(values(p0)), min0, m0))

    end

    results_matrix = reduce(vcat, transpose.(results))
    # Guardar los resultados en un archivo CSV
    writedlm("optimization_ATFv1_p03_f1_3.csv", results_matrix, ',')

    return results_matrix
    println(p0)
end



# Condiciones iniciales para el primer sistema
u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Y, U, W, C
mY_range = 10.0 .^ range(-3, 3, length=50)  # Rango logarítmico para mY: 0.001 a 1000
# Rango de tiempo (aumentado para asegurar convergencia al SS)
tspan = (0.0, 1e8)  # Tiempo largo para alcanzar el estado estacionario
pOp  = [:mU,:mW,:eP]	# Parameters to optimize
#pMin = [-3,-3,-3],		# Minimum parameter value to explore (log10)
#pMax = [3,3,3]
#writedlm(outfile1, [vcat(string("Row"), r)],'\t');
a = -3.0  # límite inferior para la escala logarítmica, e.g., 10^1
b = 3.0  # límite superior para la escala logarítmica, e.g., 10^3
n_points = 4192 # cantidad de puntos (recomendable usar una potencia de 2)
n_params = 3

pOp  = [:mU,:mW,:eP]	# Parameters to optimize
pMin = [-3.0,-3.0,-3.0]		# Minimum parameter value to explore (log10)
pMax = [3.0, 3.0, 3.0]			# Maximum parameter value to explore (log10)
iter = 2000.0				# Number of iterations per optimization run
cov  = [0.1, 0.1, 0.1]	# Covariance to calculate parameter random walk
M    = 5				# "Mutation step size" for multiplicative random walk
println("hi")



r = optimization(p, mY_range, u0, iter, pOp, cov, M, pMin, pMax, ode_system!, ode_systemNF!)

plot(r[:,5], label="mU", lw=2, ylims=[0.001,1000])
plot!(r[:,7], label="mW", lw=2)
plot!(r[:,10], label="eP", lw=2)
yaxis!(:log10)
ylabel!("Valor de parámetros (log)")
savefig("par.png")

plot(r[:, 17], label="x<eps", lw=2, ylims=[0,50])
yaxis!(0:50)
ylabel!("eps")

#plot!(r[:,17], label="c<eps", lw=2)