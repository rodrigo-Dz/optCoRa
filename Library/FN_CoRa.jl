
# --- Main functions for CoRa analysis
# --- Rodrigo Aguilar
# --- September 2025


module fn
# --- Required libraries ---
	using DelimitedFiles
	using Distributions
	using DifferentialEquations
	using ModelingToolkit
	using Sobol
	using NLsolve
	using ForwardDiff
	using LinearAlgebra
    using Statistics
    using ParameterizedFunctions
    using Plots


# --- Solve for equilibrium (fast) ---
#     p ------ Set of parameters
#     u0 ----- Initial conditions
#     system - set of ODEs
#     Returns: 
#        [0] vector with equilibrium points
#        [1] error code 
#              0 = solver converged successfully, 
#              1 = Could not find positive steady state
#              2 = possible oscillations (positive real eigenvalues),
#              3 = Error during equilibrium search
    function find_equilibrium(p, u0, system; max_iters=1000, ftol=1e-8, xtol=1e-8)
        status = 0
        par = collect(values(p))

        # --- 1. Define condition f(u) = 0 ---
        function equilibrium_condition(F, u)
            du = similar(u)
            system(du, u, par, 0.0)
            F .= du
            return F
        end
        

        # --- 2. Solve for steady state ---
        try 
            result = nlsolve(equilibrium_condition, u0, 
                            method=:newton,
                            autodiff=:forward,
                            ftol=ftol,
                            xtol=xtol,
                            iterations=max_iters,
                            show_trace=false)
            u_ss = result.zero
        

            # If negative values, try new initial conditions
            i = 1
            while any(x -> x < 0, u_ss)
                u_ = length(u0)
                u_ = fill(10.0^i, u_)

                result = nlsolve(equilibrium_condition, u_, 
                    method=:newton,
                    autodiff=:forward,
                    ftol=ftol,
                    xtol=xtol,
                    iterations=max_iters,
                    show_trace=false)

                u_ss = result.zero

                i += 1

                if i > 8
                    @warn "Could not find positive steady state; returning NaN"
                    status = 1
                    u_ss = fill(NaN, length(u0))
                    break
                end
            end


            if !result.f_converged && !result.x_converged
                @warn "Equilibrium solver failed"
                J = ForwardDiff.jacobian(u -> begin
                    du = similar(u)
                    system(du, u, par, 0.0)
                    du
                    end, u_ss)
                l = eigvals(J)

                if any(real.(l) .> 0) && any(imag.(l) .!= 0)
                    @info "Unstable equilibrium (positive real eigenvalues) — possible oscillations"
                    status = 2
                    return u_ss, status
                else
                    @warn "Equilibrium solver failed without oscillations"
                    status = 3
                    u_ss = fill(NaN, length(u0))
                    return u_ss, status
                end
            else
                status = 0
                return u_ss, status
            end
        catch e
            @warn "Error during equilibrium search: $e; returning NaN"
            u_ss = fill(NaN, length(u0))
            status = 3
            return u_ss, status
        end
        
    end

# --- Calculate Jacobian matrix at equilibrium ---
#     u_eq --- Equilibrium point
#     p ------ Set of parameters
#     system - set of ODEs
#     Returns: 
#        Jacobian matrix

    function compute_jacobian(u_eq, p, system)
        p_values = collect(values(p))
        function equilibrium_condition(u)
            du = similar(u)
            system(du, u, p_values, 0.0)
            return du
        end
        ForwardDiff.jacobian(equilibrium_condition, u_eq)
    end

    function search_oscilations(p, u0, system)
        status = 0
        par = collect(values(p))
        J = compute_jacobian(u0, par, system)
        l = eigvals(J)
        println("eigenvalues: ", l)
        if any(real.(l) .>= 0)
            @info "Unstable equilibrium (positive real eigenvalues) — possible oscillations"
            status = 2
        else
            @warn "Negative CoRa value"
            status = 4 #negCoRa
        end
        return status
    end

# --- Solve feedback system to evaluate CoRa ---

    function solve_FB(p, u0, mm)
        # If initial conditions are NaN, set to zero
        if isnan(u0[1])
            u0 = fill(0.0, length(u0))
        end 
        SS, status  = fn.find_equilibrium(p, u0, mm.FB) 
        if status != 0
            FB = NaN
        else
            FB = mm.out_FB(SS) 
        end

        return SS, FB, status
    end


# --- Evaluate CoRa ---
#     p ----- System parameters
#     u0 ---- Initial conditions
#     mm ---- Models 
#     pert -- Perturbation details
#     Returns: 
#        CoRa value, 
#        steady state of the system, 
#        steady state of the controlled variable,
#        error code 
#              0 = solver converged successfully, 
#              1 = Could not find positive steady state
#              2 = possible oscillations (positive real eigenvalues),
#              3 = Error during equilibrium search
#              4 = negative CoRa value
#              5 = CoRa value greater than 1
#              6 = NF and FB steady states differ too much

    function evalCoRa(p, u0, mm, pert)
        copy_p = copy(p)
        
    # Solve Feedback system
        SS, FB, status = solve_FB(p, u0, mm)

    # Calculate CoRa only if no errors
        if status == 0  
        # Solve NF system    
            mm.localNF(p,SS)      # Adjust parameters for no feedback system
            SS_nFB, _ = fn.find_equilibrium(p, SS, mm.nFB)
            nFB = mm.out_FB(SS_nFB)
        # If relative difference between FB and nFB is less than 0.0001, proceed to perturbation   
            if (abs(FB - nFB))/ FB < 0.01  
            # perturbation to FB
                p = copy(copy_p)
                p[pert.p] = p[pert.p]*pert.d
                SS_FBp, _ = fn.find_equilibrium(p, SS, mm.FB)
                FB_p = mm.out_FB(SS_FBp)

            # perturbation to nFB    
                p = copy(copy_p)
                p[pert.p] = p[pert.p]*pert.d
                mm.localNF(p,SS)
                SS_nFBp, _ = fn.find_equilibrium(p, SS, mm.nFB)
                nFB_p = mm.out_FB(SS_nFBp)

            # Calculate CoRa
                if FB > 0 && nFB > 0 && FB_p > 0 && nFB_p > 0
                    CoRa = log10(FB_p/FB) / log10(nFB_p/nFB)
                    #CoRa = (log(FB_p) - log(FB)) / (log(nFB_p) - log(nFB))
                    if CoRa < 0
                        status = fn.search_oscilations(p, SS, mm.FB) 
                        CoRa = NaN 
                        FB = NaN
                    elseif CoRa > 1
                        @warn "CoRa value greater than 1"
                        CoRa = NaN 
                        FB = NaN
                        status = 5
                    end
                else
                    @warn "One of the steady states is non-positive"
                    CoRa = NaN 
                    FB = NaN
                    status = 1    
                    return CoRa, SS[:, end], FB, status
                end
        # If not, report error
            else
                @warn "NF and FB steady states differ too much"
                CoRa = NaN 
                FB = NaN
                status = 6    
            end 
        else 
            CoRa = NaN
            FB = NaN
        end
        return CoRa, SS[:, end], FB, status
    end



# --- Generate CoRa curve ---
#     p ----- System parameters
#     u0 ---- Initial conditions
#     mm ---- Models
#     pert -- Perturbation details
#     Returns:
#        CoRa curve,
#        steady states,
#        errors

    function CoRacurve(p, u0, mm, pert)
        r = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
        curve = zeros(length(r))
        SSs = zeros(length(r))
        errors = zeros(length(r))
    # Iterate over perturbation range
        for j in 1:length(r)
            p[pert.c] = r[j]
            CoRa = evalCoRa(p, u0, mm, pert)
            curve[j] = CoRa[1]
            u0 =  CoRa[2]
            SSs[j] = CoRa[3]
            errors[j] = CoRa[4]
        end  
        return curve, SSs, errors
    end



# --- Calculate metrics of CoRa curve ---
#     curve -- CoRa curve
#     SSs ---- steady states of controlled variable
#     pert --- Perturbation details
#     Returns:
#        [1] robustness (fraction of points with CoRa<=eps),
#        [2] minRange (minimum rho with CoRa<=eps),
#        [3] maxRange (maximum rho with CoRa<=eps),
#        [4] min(CoRa),
#        [5] optimalRho (rho with min(CoRa)),
#        [6] neg (number of points with negative SS),
#        [7] os (number of points with oscillations),
#        [8] e (number of points with other errors),
#        [9] ss (mean steady state of controlled variable for points with CoRa<=eps)
#        [10] negative_coRa (number of points with negative CoRa)
#        [11] greater1_coRa (number of points with CoRa greater than 1)
#        [12] fb and nfb not equal

    function metrics(curve, SSs, errors, pert)
        neg_sol =     count(x -> x == 1, errors)        
        oscilations = count(x -> x == 2, errors)         #oscillations: number of points with oscillations
        other_error = count(x -> x == 3, errors)         #other: number of points with other errors
        neg_coRa =    count(x -> x == 4, errors)         #negative CoRa: number of points with negative CoRa
        greater1 =    count(x -> x == 5, errors)         #greater than 1: number of points with CoRa > 1
        not_equal =   count(x -> x == 6, errors)         #greater than 1: number of points with FB and nFB steady states differ too much

        r = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

		i = curve .<= pert.eps;                # indexes less than eps
		j = findall(i);                        # values of rho less than eps
		
        if isempty(j)
            x = copy(curve);
            x[x .=== NaN] .= Inf;                  # replace NaN with Inf for min calculation

			return [0, NaN, NaN, minimum(x), r[argmin(x)], neg_sol, oscilations, other_error, neg_coRa, greater1, not_equal,NaN]
        else
            x = copy(curve);
            x[x .=== NaN] .= Inf;                  # replace NaN with Inf for min calculation

            filtered_SSs = SSs[i]                  # steady states for values with CoRa <= eps
            ss = mean(filtered_SSs)                # ss for the values below eps
            return [length(curve[i])/length(curve), r[j[1]], r[j[end]], minimum(x), r[argmin(x)],  neg_sol, oscilations, other_error, neg_coRa, greater1, not_equal, ss]
		end    
    end



# --- Explore ---
#     iARG -- Main arguments
#     mm ---- Models
#     p ----- Core parameters
#     pert -- Perturbation details
#     expl -- Exploration details
#     u0 ---- Initial conditions
#     Returns:
#        Output file with exploration results

function explore(iARG, mm, p, pert, expl, u0)
    ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
    n = length(expl.pOp)
    
    # Open output file
    open(string("./Output/OUT_ExplCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
        if (expl.prtD==1)
            writedlm(io, [vcat([string(param) for param in expl.pOp],"robustness","minRange", "maxRange","min_CoRa", "optimalRho", "negative_sol", "oscilations", "other_errors", "negative_CoRa", "Greater_CoRa",  "not_same", "steady_state", ran)], '\t')
        else
            writedlm(io, [vcat([string(param) for param in expl.pOp],"robustness","minRange", "maxRange","min_CoRa", "optimalRho", "negative_sol", "oscilations", "other_errors", "negative_CoRa", "Greater_CoRa",  "not_same", "steady_state")], '\t')
        end
    
        if length(expl.pOp) != 2
            # Para más de 2 parámetros: usar Sobol
            sobol = SobolSeq(n)
            sobol_p = [
                round.(10.0 .^ (expl.pMin .+ (expl.pMax .- expl.pMin) .* next!(sobol)), digits=4)
                for _ in 1:expl.n_points
            ]
        else                
            # Para 2 parámetros crear grid regular
            n_per_dim = ceil(Int, sqrt(expl.n_points))
            
            g1 = range(expl.pMin[1], expl.pMax[1], length=n_per_dim)
            g2 = range(expl.pMin[2], expl.pMax[2], length=n_per_dim)
    
            sobol_p = [
                round.(10 .^ [x, y], digits = 4)
                for x in g1 for y in g2
            ]
        end
        
        n_f = length(sobol_p)
        
        # Iterate over set of parameters
        p_orig = copy(p)
        for i in 1:n_f
            println("Done ", i - 1, " out of ", n_f)
            
            # Update parameter values
            for j in 1:length(expl.pOp)
                p[expl.pOp[j]] = sobol_p[i][j]
            end
            
            # Calculate CoRa curve
            curve = CoRacurve(p, u0, mm, pert)
            m = metrics(curve[1], curve[2], curve[3], pert)    # Calculate metrics of curve
            p = copy(p_orig)
    
            # If printing full CoRa curve specified:
            if (expl.prtD==1)	
                writedlm(io, [vcat(sobol_p[i], m, curve[1])],'\t')
            else
                writedlm(io, [vcat(sobol_p[i], m)],'\t')
            end
        end
    end
end


# --- Optimize ---
#     iARG -- Main arguments
#     mm ---- Models
#     p ----- Core parameters
#     pert -- Perturbation details
#     opt --- Optimization details
#     u0 ---- Initial conditions
#     Returns:
#        Output file with optimization results

    function optimize(iARG, mm, p, pert, opt, u0)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)

        open(string("./Output/OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
            if (opt.prtD==1)
                    writedlm(io, [vcat("Iteration", [string(param) for param in opt.pOp],"robustness","minRange", "maxRange","min_CoRa", "optimalRho", "negative_sol", "oscilations", "other_errors", "negative_CoRa", "Greater_CoRa",  "not_same", "steady_state", ran)], '\t')
            else
                    writedlm(io, [vcat("Iteration", [string(param) for param in opt.pOp],"robustness","minRange", "maxRange","min_CoRa", "optimalRho", "negative_sol", "oscilations", "other_errors", "negative_CoRa", "Greater_CoRa",  "not_same", "steady_state")], '\t')
            end

            if opt.rand == 1
                for i in 1:length(opt.pOp)
			    	p[opt.pOp[i]] = round.(10 .^ (rand(Uniform(opt.pMin[i], opt.pMax[i]))), digits = 4);
                end
		    end

            p_copy = copy(p)

            curve0 =  CoRacurve(p, u0, mm, pert)
            m0 = metrics(curve0[1], curve0[2], curve0[3], pert)    # Calculate metrics of curve

        # Initial metrics

			op0 = log10(m0[3]/m0[2])    # Property to optimize, range with CoRa<=eps 
            prop0 = m0[1]
			mi0 = m0[4]                 # Minimum CoRa
            r0 = zeros(length(opt.pOp)) 
            s0 = m0[9]                  # Steady state of controlled variable for points with CoRa<=eps
                

            if mi0 == NaN
                println("Negative solutions, try to start with anther parametes")
                return NaN
            end
            if m0[6] != 0 
                println("Negative solutions, try to start with anther parametes")
                return NaN
            end
            if m0[7] != 0 # 
                println("Oscilations, try to start with anther parameters")
                return NaN
            end
            if m0[8] != 0 # there are other errors
                println("Unknown errors, try to start with anther parameters")
                return NaN
            end


            if (opt.prtD==1)	
                writedlm(io, [vcat(0, [p[i] for i in opt.pOp],  m0, curve0[1])],'\t')
            else			
                writedlm(io, [vcat(0, [p[i] for i in opt.pOp], m0)],'\t')
            end

        # Iterate optimization
            p = copy(p_copy)

            ss_locked = false

            for i in 1:opt.iter
                println(i)
            # Generate random changes in parameters according to a normal distribution
                rI = rand(MvNormal(zeros(length(opt.pOp)), zeros(length(opt.pOp)) .+ opt.cov)) 
			# Update parameter values	
                for pI in 1:length(opt.pOp)  
					r0[pI] = p[opt.pOp[pI]]                # Save previous value
					p[opt.pOp[pI]] *= (opt.M .^ rI[pI])    # Update value
                    p[opt.pOp[pI]] = round(p[opt.pOp[pI]]; sigdigits=4)
                    #println("p[",opt.pOp[pI],"] = ", p[opt.pOp[pI]])
				# Exclude values outside regime of exploration:
					if p[opt.pOp[pI]] < (10.0 ^ opt.pMin[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMin[pI])
					elseif p[opt.pOp[pI]] > (10.0 ^ opt.pMax[pI])
						p[opt.pOp[pI]] = (10.0 ^ opt.pMax[pI])
					end
				end

            # Calculate new CoRa curve
                curve1 = CoRacurve(p, u0, mm, pert)
                m1 = metrics(curve1[1], curve1[2], curve1[3], pert)    # Calculate metrics of curve
                op1 = log10(m1[3]/m1[2])  
                prop1 = m1[1]
                mi1 = m1[4]
                s1 = m1[9]



            # C1 = minimize min(CoRa) when op0 and op1 = NaN
            # NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)
                xiC = (mi0 ^ 2) / (2 * 0.083)
                xiP = (mi1 ^ 2) / (2 * 0.083)
				c1 = isnan(op0+op1) && (rand() < exp((xiC - xiP))) 


            # C8 = stabilize steady state
                if opt.mean == NaN
                    c8 = true
                else
                    if ss_locked == false
                        c8 = true
                    else
                        s0l = log10(s0)
                        s1l = log10(s1)
    
                        xiC = exp(-((s0l - opt.mean)^2) / (2 * opt.sd^2))
                        xiP = exp(-((s1l - opt.mean)^2) / (2 * opt.sd^2))
                        c8 = rand() < xiP / xiC
                    end
                end

#=  Trying to implement a more sophisticated version of C8 

# Required in opt:
                #   opt.ss_log_min, opt.ss_log_max  
                #   opt.ss_ramp, opt.ss_lambda, opt.ss_log_sd

                if opt.ss_log_min == NaN && opt.ss_log_max == NaN
                    c8 = true
                else
                    lmin = opt.ss_log_min
                    lmax = opt.ss_log_max

                    s0l = log10(s0)
                    s1l = log10(s1)

                    # lock only if step is accepted)
                    hit_range = (!ss_locked) && in_range(s1l)

                    # If not locked, no constraint
                    if !ss_locked
                        c8 = true
                    else
                        d0 = dist_to_range(s0l)
                        d1 = dist_to_range(s1l)

                        ramp = haskey(opt, :ss_ramp) ? max(opt.ss_ramp, 1) : 1
                        lmax = haskey(opt, :ss_lambda) ? opt.ss_lambda : 1.0
                        l = lmax * min(1.0, i / ramp)

                        E0 = (d0 / ss_log_sd)^2
                        E1 = (d1 / ss_log_sd)^2

                        c8 = (E1 <= E0) || (rand() < exp(-l * (E1 - E0)))
                    end
                end
=#

            # C2 = maximize range using the range
			# NOTE: As op0,op1=[0,rrO], but still correct exponential with the expected variance of ~U(0,1)
			#       !! ~U(0,1)*(rrO^2) variance resulted in very noisy runs...
				rrO = pert.r[2] - pert.r[1]
            	xiC = (rrO - op0) / (2 * 0.083)
            	xiP = (rrO - op1) / (2 * 0.083)

				c2 = rand() < exp((xiC - xiP)) && c8

            # C2 = maximize range using the proportion; we used variance of 0.01  wich gives a good balance between exploration and explotation
            #   xiC = ((prop0 - 1)^2) / (2 * 0.01)
            #   xiP = ((prop1 - 1)^2) / (2 * 0.01)
			#	c2 = rand() < exp((xiC - xiP))

            # C3 = no negative solutions
                c3 = m1[6] == 0
            # C4 = no oscilations
                c4 = m1[7] == 0
            # C5 = no other errors
                c5 = m1[8] == 0 
            # C6 = no negative CoRa
                c6 = m1[10] == 0
            # C7 = no CoRa greater than 1
               # c7 = m1[11] == 0 

            # If conditions met to accept new parameter values
				if((c1 || c2 ) && c3 && c4 && c5 && c6 && c8) 
					op0 = op1
                    prop0 = prop1
					mi0 = mi1
                    print_curve = curve1[1]
                    s0 = s1
                    m0 = m1

                    in_range(x) = (x >= (opt.mean - opt.sd)) && (x <= (opt.mean - opt.sd))

                    if in_range(log10(s0))
                        ss_locked = true
                    end
                    println("ok")
				else
				# If not, revert to previous parameter values
					for pI in 1:length(opt.pOp)
						p[opt.pOp[pI]] = r0[pI]
					end
                    println("no ok")
                    print_curve = curve0[1]
				end

            # Save results of each iteration
                if (opt.prtD==1)
                    writedlm(io, [vcat(i,[p[x] for x in opt.pOp],m0, print_curve)],'\t')
                else			
                    writedlm(io, [vcat(i,[p[x] for x in opt.pOp],m0)],'\t')
                end
            end
        end
    end



# --- Dynamics ---
#     mm ---- Models
#     u0 ---- Initial conditions
#     p ----- Core parameters
#     pert -- Perturbation details
#     dyn --- Dynamics details
#     iARG -- Main arguments
#     Returns:
#        Dynamics plot

    function dynamics(mm, u0, p, pert, dyn, iARG)
        # Extraer parámetros 
        p_values = collect(values(p))
        keys_list = collect(keys(p))
        i_pert = findfirst(isequal(pert.p), keys_list)
        println("Perturbing parameter: ", keys_list[i_pert], " at time ", dyn.pert_time, " with factor ", dyn.pert_size)
        
        # --- Search for SS ---
        SS, _ = fn.find_equilibrium(p, u0, mm.FB)
        i = 0
        while any(<(0), SS) && i < 5
            println("Negative SS found, retrying with new initial condition (10^$i)")
            u_try = fill(10.0^i, length(u0))
            try
                SS, _ = fn.find_equilibrium(p, u_try, mm.FB)
            catch
                SS .= NaN
            end
            i += 1
        end
    
        if any(isnan, SS) || any(<(0), SS)
            println("Could not find positive steady state, dynamics not simulated")
            return nothing
        end
    
        println(" Steady state found: ", SS)
    
        # --- Definir perturbación ---
        affect!(integrator) = integrator.p[i_pert] *= dyn.pert_size
        cb = PresetTimeCallback(dyn.pert_time, affect!)
    
        # --- Simulación con feedback ---
        prob_FB = ODEProblem(mm.FB, SS, dyn.tspan, p_values)
        sol_FB = solve(prob_FB, Tsit5(), callback=cb)
    
        # --- Simulación sin feedback ---
        mm.localNF(p, SS)           # Ajusta parámetros sin feedback
        p_values_nFB = collect(values(p))
        prob_nFB = ODEProblem(mm.nFB, SS, dyn.tspan, p_values_nFB)
        sol_nFB = solve(prob_nFB, Tsit5(), callback=cb)

        plt = plot(sol_FB.t, sol_FB[dyn.plot, :],
         xlabel="Time", ylabel="Concentration",
         label="Feedback", legend=:best)

        plot!(sol_nFB.t, sol_nFB[dyn.plot, :], linestyle=:dash, label="No Feedback")

        savefig("./Output/OUT_dynamics_$(iARG.mm)_$(iARG.ex)_$(iARG.pp)_$(iARG.ax).png")
        return ("Yea")
    end


    
# --- CoRa Curve ---
#     p ----- Core parameters
#     u0 ---- Initial conditions
#     mm ---- Models
#     pert -- Perturbation details
#     iARG -- Main arguments
#     Returns:
#        CoRa curve plot

    function curve(p, u0, mm, pert, iARG)
        ran = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
        c = CoRacurve(p, u0, mm, pert)
        println("SS: ", c[2])
        println("Status: ", c[3])
        plot(ran, c[1];
            xscale = :log10,
            xlims = (minimum(ran), maximum(ran)),
            ylims = (0, 2),
            xlabel = "P["*string(pert.c)*"]",
            ylabel = "CoRa")
        savefig(string("./Output/OUT_curve_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".png"))
    end





    function sensitivity2(p, u0, mm, pert, iARG, sens)
        dir = "./Output/OUT_sensitivity_$(iARG.mm)_$(iARG.ex)_$(iARG.pp)_$(iARG.ax)/"
        if !isdir(dir)
            mkdir(dir)
        end

        for i in 1:length(sens.pSens)
            ranP = 10 .^ range(sens.r[1], sens.r[2], length=sens.coras)
            ranC = 10 .^ range(pert.r[1], pert.r[2], length=pert.coras)
            curve_matrix = zeros(length(ranP), length(ranC))
            ss_matrix = zeros(length(ranP), length(ranC))
            for j in 1:length(ranP)
                p[sens.pSens] = ranP[j]
                curve = zeros(length(ranC))
                c = CoRacurve(p, u0, mm, pert)
                curve_matrix[j, :] = c[1]
                ss_matrix[j, :] = log10.(c[2])
            end
            heatmap(ranC, ranP, curve_matrix;
                xscale = :log10,
                yscale = :log10,
                xlabel = "P["*string(pert.c)*"]",
                ylabel = "P["*string(sens.pSens)*"]",
                colorbar_title = "CoRa",
                c = reverse(palette(:viridis)),
                clim = (0,1))
            savefig(string(dir,sens.pSens,"_CoRa.png"))
            heatmap(ranC, ranP, ss_matrix;
                xscale = :log10,
                yscale = :log10,
                xlabel = "P["*string(pert.c)*"]",
                ylabel = "P["*string(sens.pSens)*"]",
                colorbar_title = "SS",
                c = reverse(palette(:viridis)))
            savefig(string(dir,sens.pSens,"_SS.png"))
        end
    end




    using Plots

    function sensitivity(p, u0, mm, pert, iARG, sens)
        dir = "./Output/OUT_sensitivity_$(iARG.mm)_$(iARG.ex)_$(iARG.pp)_$(iARG.ax)/"
        isdir(dir) || mkdir(dir)

        ranP = 10 .^ range(sens.r[1], sens.r[2], length = sens.coras)
        ranC = 10 .^ range(pert.r[1], pert.r[2], length = pert.coras)

        for psym in sens.pSens   
            # guarda valor original para restaurar
            p0 = p[psym]

            curve_matrix = zeros(length(ranP), length(ranC))
            ss_matrix    = zeros(length(ranP), length(ranC))

            for (j, valP) in enumerate(ranP)
                p[psym] = valP

                c = CoRacurve(p, u0, mm, pert)  
                curve_matrix[j, :] .= c[1]
                ss_matrix[j, :] .= log10.(c[2])
            end

            # restaura el parámetro
            p[psym] = p0

            # --- Heatmap CoRa ---
            heatmap(
                ranC, ranP, curve_matrix;
                xscale = :log10,
                yscale = :log10,
                xlabel = "P[" * string(pert.c) * "]",
                ylabel = "P[" * string(psym) * "]",
                colorbar_title = "CoRa",
                c = reverse(palette(:viridis)),
                clim = (0, 1)
            )
            savefig(string(dir, string(psym), "_CoRa.png"))

            # --- Heatmap SS ---
            heatmap(
                ranC, ranP, ss_matrix;
                xscale = :log10,
                yscale = :log10,
                xlabel = "P[" * string(pert.c) * "]",
                ylabel = "P[" * string(psym) * "]",
                colorbar_title = "SS",
                c = reverse(palette(:viridis))
            )
            savefig(string(dir, string(psym), "_SS.png"))
        end

        return nothing
    end

end