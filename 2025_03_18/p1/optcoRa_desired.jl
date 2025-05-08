###CoRa_v2.5

###Alright, once and for all
###In the spirit of proper documentation, I'll be as explicit as I can with these comments
###Taken from v2, no need to modify this just yet
## Load functions & parameters:
using DelimitedFiles;
using Distributions;
using CSV;
using Sobol;
using NLsolve
	
#This includes the model to study itself, it must be prepared before running CoRa, with the explicit constraint that the number of d/dt are the same in the FB and NF equations
mm = include(string("Library/Md_",iARG.mm,".jl"));
###This next line is changed to use the "updated" functions .jl
fn = include(string("Library/FN_CoRa_v2.jl"));
## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file, pp : Label for perturbation type, an : Chose analysis type);
include(string("InputFiles/ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl")); # Perturbation details
include(string("InputFiles/ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters
include(string("InputFiles/ARGS_",iARG.mm,"_OptCoRa_",iARG.ex,".jl"))
key_names = (:g, :mY, :gY, :mU, :gU, :mW, :gW, :e0, :eP, :eM, :mUs); #The names of the parameters of the systems
strict = true;  # Should a steady state be found, and then refound in the next iteration, and this difference not be in accordance to the rtol given, said SS will be obtained regardless
gap_size = 20; # This is the size of a region to be examined for NaN gaps
gap_tol = 0.5; # How much of the examined region can be NaNs before the program kills it
rtol = 1e-10;
pO = copy(p);
x0 = zeros(length(mm.odeFB.syms));
rtol = rtol;


optimal = 10
renge = [optimal+(optimal*1.5), optimal-(optimal*1.5) ]

Threads.@threads for ruN in 1:mrw.runs
    eps = 0.1;
    n = 30.0;
    open(string("Results_optimization/OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,"_",ruN,".txt"), "w") do io
        if (mrw.prtD==1)
            writedlm(io, [vcat("Run","Iteration",[string(param) for param in mrw.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRa)",10 .^ collect(pert.r[1]:((abs(pert.r[1]) + abs(pert.r[2]))/n):pert.r[2]))], '\t')
        else
            writedlm(io, [vcat("Run","Iteration",[string(param) for param in mrw.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRa)")], '\t')
        end
        println("RUN #",ruN)
        p = copy(pO);
        # If random initial parameters specified:
        if (mrw.ranP==1)
            os = 1
            while os != 1
                fn.mrwR(mrw,p);
                os = fn.find_oscilations(mm,p)
            end
        end
 
        os = fn.find_oscilations(mm,p)
        if os == true
            println("choose other initial parameters")
            exit()
        end

        # Initialize system:
        CoRa0 = fn.CoRaCurve(p, pert, mm, n, gap_size, gap_tol, x0, rtol, strict)
        op0, mi0, cost0, w0 = fn.CoRam(CoRa0,eps)
 
                # Print initial state:
        if (mrw.prtD==1)	# If printing full CoRa curve specified:
            writedlm(io, [vcat(ruN,0,[p[i] for i in mrw.pOp],op0,mi0,CoRa0)],'\t')
        else			# Else:
            writedlm(io, [vcat(ruN,0,[p[i] for i in mrw.pOp],op0,mi0)],'\t')
        end
        # Optimization iterations
        println("I: minCoRa = ",mi0,"\t |CoRa<eps| = ",op0)
  

        i = 1
        eps = 0.1
        initial_eps = eps
        #temperature = 1.0  # Initial temperature for simulated annealing
        #cooling_rate = 0.995  # Cooling rate for temperature
        while (i != mrw.iter)  #size(pars, 1)   # Número de filas     
            p0 = copy(p)
            
            choose = 0
            while choose == 0
                fn.mrwR(mrw,p);
                os = fn.find_oscilations(mm,p)
                if os == true
                    p = p0   #revert
                else 
                    choose = 1
                end
            end

            # curve
            CoRa1 =  fn.CoRaCurve(p, pert, mm, n, gap_size, gap_tol, x0, rtol, strict)
            println("el eror esta aqui")
            # metrics
            # op1: points below eps / n
            # mi1: min point
            # cost1: castigo para los valores feos
            op1, mi1, cost1, w1 = fn.CoRam(CoRa1,eps)           
            println("o aqui")

            c1 = (mi1>=0);
		## If CoRa>eps for all conditions, evaluate the min(CoRa) for both sets:
		### NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)

		## If CoRa>=eps for some conditions, evaluate the |CoRa<=eps| for both sets:
		### NOTE: As op0,op1=[0,1], correct exponential with the expected variance of ~U(0,1)
			xiC = ((1 - op0)^2) / (2 * 0.083) + cost0 + w0;
			xiP = ((1 - op1)^2) / (2 * 0.083) + cost1 + w1 ; 
		    c3 = rand() < exp(xiC - xiP);
            
		# Evaluate if accept new parameter values or not:
		## Only accept in the regime of interest, i.e. CoRa>=0:
		## If CoRa>eps for all conditions, evaluate the min(CoRa) for both sets:
		### NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)

		## If CoRa>=eps for some conditions, evaluate the |CoRa<=eps| for both sets:
		### NOTE: As op0,op1=[0,1], correct exponential with the expected variance of ~U(0,1)

            if(c1 && (c3))
                # If yes, update "reference" system
                CoRa0 = CoRa1;
                op0 = op1;
                mi0 = mi1;
                cost0 = cost1
                w0 = w1
            else
                    # If not, revert to previous parameter values
                p = p0
            end

            writedlm(io, [vcat(ruN, i, [p[i] for i in mrw.pOp],op0,mi0,CoRa0)],'\t')
            flush(io)
            println(string("Line ", i, " done!"))
            i = i + 1
            #temperature *= cooling_rate
            #if op0 > 0.5  # Si se acerca al óptimo, reducir eps para mejorar ajuste fino
            #    eps = max(eps * 0.95, initial_eps * 0.01)  # Reducir eps pero con un límite inferior
            #    mrw.cov .= [x * 0.5 for x in mrw.cov]  # Reducir el tamaño del paso
            #end  # Reduce temperature
            #if i == mrw.iter
	        #    println("iteraciones alcanzadas")
            #    break
            #end     
        end
    end
end
