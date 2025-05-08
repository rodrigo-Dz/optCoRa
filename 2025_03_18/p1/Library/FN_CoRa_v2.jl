# Steady state & CoRa calculation functions
#		Julia v.1.8

module fn
	# Required libraries
	using DelimitedFiles
	using Distributions
	using DifferentialEquations
	using ModelingToolkit
	using Sobol
	using NLsolve
	using ForwardDiff
	using LinearAlgebra

	###This function is rather simple, but it's for readability in the main code purposes, and prevention of typos
	function Restart(a, b)
		a = zeros(length(a)).+NaN
		b = zeros(length(b)).+NaN
		return a, b
	end

	###Can't have a steady state & check function without the steady state and the checking, now can we?
	function SS(syst, p, x0, rtol, strict, tag, k)
		pV = [p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps]
		tS = 0
		dXrm = 1
		while(dXrm > rtol)
			if(mod(tS, 1e7) == 0 && tS > 0)
				println("I'm on the ", tS / 1e6, "th attempt of this condition (", tag, ")")
			end
			x1 = x0
			ss = try
				fn.solve(fn.ODEProblem(syst,x0,1e6,pV); reltol=rtol)
			catch
			end
			if ss.retcode != :Success
				try
					ss = fn.solve(fn.ODEProblem(syst,x0,1e6,pV),alg_hints=[:stiff]; reltol=rtol)
				catch err
					println("WARNING: Error in ODE simulation: <<",err,">>. ss --> NaN")
					x0 = zeros(length(syst.syms)).+NaN
					break
				end
			end
			dXrm = maximum(abs.(big.(ss(1e6))-big.(ss(1e6-0.01)))./big.(ss(1e6)))
			x0 = ss(1e6)
			if(strict == false && x1 == x0)
				println("WARNING: System's solution is not strict enough according to rtol. Returning Steady State obtained")
				return x0
			end
			tS += 1e6
			if(tS>2e8)	# e7
				println("WARNING: Maximum iteration reached (simulated time 2e7). Returning NaN for k = ", k,". Max relative Delta: ",dXrm)
				x0 = zeros(length(syst.syms)).+NaN
				break
			end
		end
		return x0
	end

	###A small function to check for NaNs, again made for the purposes of readability within the code
	function NaNCheck(a, b)
		if(any(isnan.(a)))
			a, b = fn.Restart(a, b)
			return a, b
		else
			return "Valid"
		end
	end

	###Now we create the Check function!
	function Check(ssR, soR, rtol, mm)
		if(any.(isnan.(mm.outFB(ssR))) || any.(isnan.(mm.outFB(soR))) || (abs(mm.outFB(ssR) - mm.outNF(soR)) > 1e-4))
			rtol *= 1e-3
			if(rtol < 1e-24)
				println("ERROR: Check NF system (reltol=",rtol,").")
				println(vcat(pert.p,i,[p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps],mm.outFB(ssR),mm.outNF(soR)))
				if(abs(mm.outFB(ssR) - mm.outNF(soR))/mm.outFB(ssR) > 0.01)
					ssR, soR = Restart(ssR, soR)
					println("Error too large. SS results excluded!")
				end
			end
			return ssR, soR, rtol, "Insufficient"
		else
			return ssR, soR, rtol, "Sufficient"
		end
	end;

	###Here it is, the SS&Check function. It will, of course, be built upon an SS() and a Check() funciton previously defined within this document
	function SSandCheck(p, x0, rtol, mm, strict, k)
		###The rtol value basically states that this will be attempted 5 times: This is because of the comparison made in "if(abs(mm.outFB(ssR) - mm.outNF(soR)) > 1e-4)". If that is successful
		###Then the value of rtol is multiplied by 1e-3, until it is no longer greater or equal than 1e-24, as stated in our while() condition
		flag = "Insufficient"
		ssR, soR = fn.Restart(x0, x0)
		println(x0)
		while(rtol >= 1e-24 && flag == "Insufficient")
			# Reference steady state:
			ssR = fn.SS(mm.odeFB, p, x0, rtol, strict, "ssR", k)
			if(fn.NaNCheck(ssR, soR) != "Valid")
				println("Condition excluded! ssR --> NaN")
				break
			end
			# Locally analogous system reference steady state:
			mm.localNF(p,ssR)
			###Of note here, instead of using the initial condition x0, we use ssR.
			soR = fn.SS(mm.odeNF, p, ssR, rtol, strict, "soR", k)
			if(fn.NaNCheck(soR, ssR) != "Valid")
				println("Condition excluded! soR --> NaN")
				break
			end
			ssR, soR, rtol, flag = Check(ssR, soR, rtol, mm)
		end
		return ssR, soR, rtol
	end

	function Perturbation(ssR, soR, p, rtol, mm, pert, strict, k)
		p[pert.p] *= pert.d
		ssD = fn.SS(mm.odeFB, p, ssR, rtol, strict, "ssD", k)
		soD = fn.SS(mm.odeNF, p, soR, rtol, strict, "soD", k)
		p[pert.p] /= pert.d
		if(fn.NaNCheck(ssD, soD) != "Valid")
			println("Condition excluded! ssD --> NaN")
		elseif(fn.NaNCheck(soD, ssD) != "Valid")
			println("Condition excluded! SoD --> NaN")
		end
		return ssD, soD
	end

	function CoRa(ssR, ssD, soR, soD)
		if abs(log10(soD/soR)) < 1e-6
			return NaN
		end
		return log10(ssD/ssR)/log10(soD/soR)
	end

	function CoRac(p, pert, mm, x0,n)
		s = (abs(pert.r[1]) + abs(pert.r[2]))/n;

		r = 10 .^ collect(pert.r[1]:s:pert.r[2])
		CoRaValues = ones(length(r)) .+ Inf
		p[pert.p] = pert.c
		for i in 1:lastindex(r)
			p[pert.p] *= r[i]
			rtol = 1e-12
			ssR, soR = ones(length(mm.odeFB.syms)), ones(length(mm.odeNF.syms))
			try
				ssR, soR, rtol = SSandCheck(p, x0, rtol, mm)
				ssD, soD = fn.Perturbation(ssR, soR, p, rtol, mm, pert)
				CoRaValues[i] = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
			catch err
				println("WARNING: Error in ODE simulation: <<",err,">>. CoRa --> NaN")
				CoRaValues[i] = NaN
			end
			# Perturbation:
			p[pert.p] /= r[i]
		end
		return CoRaValues
	end

	function Dyn(syst, p, x0, tspan, rtol)
		pV = [p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps]
		xD = try
			fn.solve(fn.ODEProblem(syst,x0,tspan,pV); reltol=rtol)
		catch
			try
				fn.solve(fn.ODEProblem(syst,x0,tspan,pV),alg_hints=[:stiff]; reltol=rtol)
			catch err
				println("WARNING: Error in ODE simulation: <<",err,">>. ss --> NaN")
				zeros(length(syst.syms)).+NaN
			end
		end
		if xD.retcode == :Unstable
			try
				xD = fn.solve(fn.ODEProblem(syst,x0,tspan,pV),alg_hints=[:stiff]; reltol=rtol)
			catch err
				println("WARNING: Error in ODE simulation: <<",err,">>. ss --> NaN")
				zeros(length(syst.syms)).+NaN
			end
		end
		return(xD)
	end


	# MRW random initial parameters
	# INPUT: mrw - Handle for the optimization (MRW) details
	#        p   - Dictionary function with the ODE parameters & values
	# OUPUT:     - Updated p
	function mrwR(mrw,p)
		for i in 1:length(mrw.pOp)
			#p[mrw.pOp[i]] = 10.0 .^ (rand(Uniform(mrw.pMin[i], mrw.pMax[i])));
			p[mrw.pOp[i]] = 1000.0 .* (rand(Uniform(0, 1)));
		end
	end



	function mrwPar(mrw,p)
		r0 = zeros(length(mrw.pOp));			# Vector of parameters to optimize in reference system
		rI = rand(MvNormal(zeros(length(mrw.pOp)), zeros(length(mrw.pOp)) .+ mrw.cov)); # Random values to update parameters
		for pI in 1:length(mrw.pOp)
			r0[pI] = p[mrw.pOp[pI]];			 # Save previous value
			p[mrw.pOp[pI]] *= (mrw.M .^ rI[pI]); # Update value
			p[:mUs] = NaN
		#p[mrw.pOp[pI]] *= (M .^ rI[pI]); # Update value
		# Exclude values outside regime of exploration:
			if p[mrw.pOp[pI]] < (10.0 ^ mrw.pMin[pI])
				p[mrw.pOp[pI]] = (10.0 ^ mrw.pMin[pI])
			elseif p[mrw.pOp[pI]] > (10.0 ^ mrw.pMax[pI])
				p[mrw.pOp[pI]] = (10.0 ^ mrw.pMax[pI])
			end
		end
	end


	function CoRaCurve(p, pert, mm, n, gap_size, gap_tol, x0, rtol, strict)
		s = (abs(pert.r[1]) + abs(pert.r[2]))/n;
        r = 10 .^ collect(pert.r[1]:s:pert.r[2]);
        r = round.(r, sigdigits = 10);

		CoRa = zeros(length(r)).+ NaN;
		###So, then, this will repeat the for loop in however many steps for the parameters
		k = 1;
		try
			for k in 1:lastindex(r)
				if ((k > gap_size) && ((sum(isnan.(CoRa[(k-gap_size):k])))/(gap_size)) >= gap_tol)
					println("Set", i, " has failed in at least", gap_size * gap_tol, " out of ", gap_size , "continuous points, so it will no longer be computed")
					CoRa = zeros(length(r)).+ NaN
					break
				end
				###We reset our perturbation for the next run of the loop
				#p[pert.c] = p0[pert.c];
				#println(p[pert.c]);
				###The parameter to change is multiplied by the corresponding value in our steps collection, and the error tolerance is set to 1e-12 (arbitrarily)
				p[pert.c] = r[k];
				println(p[pert.c]);
				###Up next, we must find the steady states of both ssR and soR, while also checking that the process itself didn't fail. Let's make that a single, "fn.SSandCheck()" function
				ssR, soR, rtol_local = fn.SSandCheck(p, x0, rtol, mm, strict, k)
				println(ssR)
				if any(isnan.(ssR)) || any(isnan.(soR))
					continue
				end
				###Now, we have to do the Perturbation itself!
				ssD, soD = fn.Perturbation(ssR, soR, p, rtol_local, mm, pert, strict, k)
				if any(isnan.(ssD)) || any(isnan.(soD))
					continue
				end
				###And we must return CoRa now, too
				CoRa[k] = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
				###With everything done, it's time to output them into the file!
				###Check this one out
				x0 = ssR
			end
		catch
		end
		return(CoRa)
	end

	# CoRa "metrics"
	# INPUT: p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        mm    - Handle for the considered motif
	#        x0    - Vector of initial state of the ODE system
	# OUPUT: CoRas 		 - Vector of CoRa values for the range of parameters
	#   	 |CoRa<=eps| - Proportion of the CoRas vector less or equal to eps
	#		 min(CoRa)	 - Minimum CoRa value in CoRas
	function CoRam(CoRaC,eps)
		w = 0
		if any(isnan, CoRaC)
			w = 999999
		end
		return [sum(CoRaC .<= eps)/length(CoRaC), minimum(filter(!isnan,CoRaC)), sum(x^2 for x in CoRaC if x > eps), w ]
		#return [CoRas,(sum(CoRas .<= eps) + count(isnan,CoRas))/length(CoRas), minimum(filter(!isnan,CoRas))]
	end;



	function find_oscilations(mm, p)
		flag = 0
		odeFB = mm.odeFB
		u0 = [1.0, 1.0, 1.0, 1.0] 
		params = [p[:g], p[:mY], p[:gY], p[:mU], p[:gU], p[:mW], p[:gW], p[:e0], p[:eP], p[:eM], p[:mUs]]

		function equilibrium_point(odeFB, params)
			# Resolver odeFB = 0 para encontrar el equilibrio (estado estacionario)
			# Reemplaza esta parte por la llamada adecuada dependiendo de cómo esté definida la función en odeFB
			# Si f es la función que define el sistema ODE, usa esa.
			function system_eq(x)
				du = similar(x)
				odeFB.f(du, x, params, 0.0)  # Asegúrate de usar la función correcta aquí
				return du
			end
			result = nlsolve(system_eq, u0, autodiff=:forward)
			# Resolver el sistema
			return result.zero
		end
		# Encontrar el punto de equilibrio
		
		# Calcular la matriz Jacobiana en el punto de equilibrio
		def_jacobian(u) = ForwardDiff.jacobian(x -> begin
			du = similar(x)
			odeFB.f(du, x, params, 0.0)
			du
		end, u)
		
		# Bifurcation analysis: Variar mY y mW en escala logarítmica
		mY_range = 10 .^ range(-3, 3, length=50)  # Rango logarítmico para mW: 0.001 a 1000

		for mY in mY_range
			println(mY)
			params[2] = mY  # Actualizar mY
							
# Calcular la matriz Jacobiana en el estado inicial
			equilibrium = equilibrium_point(odeFB, params)

			jac_matrix_at_equilibrium = def_jacobian(equilibrium)
			eigenvalues, eigenvectors = eigen(jac_matrix_at_equilibrium) 
			oscilates = any(imag(lam) != 0 && real(lam) > 0 for lam in eigenvalues)
    		if oscilates == true
        		flag = 1  # Eigenvalues complejos con parte real positiva
    		end
		end
		return flag
	end
end




