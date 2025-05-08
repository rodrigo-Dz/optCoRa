###CoRa_v2.5

###Alright, once and for all
###In the spirit of proper documentation, I'll be as explicit as I can with these comments
###Taken from v2, no need to modify this just yet
## Load functions & parameters:
using DelimitedFiles;
using Distributions;
using CSV;
using Sobol;
using Random, Base.Threads
using Pkg;
using CSV;
using DelimitedFiles;
#Pkg.activate(".\\CoRa");		# Activate local environment (requiere '.toml' files)
iARG = (mm = "FADv1",  # Label for motif file
ex = "Fig1",      # Label for parameters file
pp = :mY,         # Label for perturbation type
ax = :mY);    # Label for condition/environment
#pars = CSV.File("InputFiles/ARGS_FADv1_Mass_Par_1250Set2.csv"); # Core parameters
key_names = (:g, :mY, :gY, :mU, :gU, :mW, :gW, :e0, :eP, :eM, :mUs); #The names of the parameters of the systems
strict = true;  # Should a steady state be found, and then refound in the next iteration, and this difference not be in accordance to the rtol given, said SS will be obtained regardless
gap_size = 100.0; # This is the size of a region to be examined for NaN gaps
gap_tol = 0.5; # How much of the examined region can be NaNs before the program kills it
rtol = 1e-10;

@time begin;
#This includes the model to study itself, it must be prepared before running CoRa, with the explicit constraint that the number of d/dt are the same in the FB and NF equations
mm = include(string("Library/Md_",iARG.mm,".jl"));
println("a")
###This next line is changed to use the "updated" functions .jl
fn = include(string("Library/FN_CoRa_v2.jl"));
println("b")
## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file, pp : Label for perturbation type, an : Chose analysis type);
include(string("InputFiles/ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl")); # Perturbation details
include(string("InputFiles/ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters
println("c")
pO = copy(p);

x0 = zeros(length(mm.odeFB.syms));
rtol = rtol;
open(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do outfile1
    r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
    r = round.(r, sigdigits = 10);
    writedlm(outfile1, [vcat(string("Row"), r)],'\t');
    println("d")
    a = 1.0  # límite inferior para la escala logarítmica, e.g., 10^1
    b = 1000.0  # límite superior para la escala logarítmica, e.g., 10^3
    n_points = 1024 # cantidad de puntos (recomendable usar una potencia de 2)
    # Crea el generador de secuencias Sobol en 3 dimensiones
    sobol_p = SobolSeq(3)
    println("q pas1")
    # Genera los puntos y aplica la transformación logarítmica
    #sobol_points_log = [10.0 .^ (a .+ (b - a) .* next!(sobol_p)) for _ in 1:n_points]
    

    sobol_points = [1000.0 .* next!(sobol_p) for _ in 1:n_points]
    # Convierte los puntos a una matriz para facilitar el manejo
    #pars = reduce(hcat, sobol_points)'
    # Convierte los puntos a una matriz para facilitar el graficado
    pars = Array(reduce(hcat,  sobol_points)')
    println("q pas2")

    # Guarda los parámetros en un archivo CSV
    output_file = "sobol_parameters.csv"
    println("q pas3")

    open(output_file, "w") do f
        writedlm(f, pars, ',')    # Escribe los datos de `pars` en el archivo
    end
    println("q pas4")

    n = 1024
    n_chunks = 2
    chunk_size = Int(ceil(n / n_chunks))
    println("q pas5")

# Crear los 4 chunks
    chunks = [pars[(i-1)*chunk_size+1:min(i*chunk_size, n), :] for i in 1:n_chunks]
    println("q pas6")


    println(pars)
    println(chunks)
    println("sobol yei")

    Threads.@threads for piece in 1:n_chunks
    open(string("chunk_ATF_",piece,".txt"), "w") do io
    println("RUN #",piece)
    chunk = chunks[piece]
    println("para chunk", piece, "el chunk es ", chunk)
    for i in 1:size(chunk, 1)   # Número de filas
        r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
        r = round.(r, sigdigits = 10);
        p[:mU] = chunk[i,1]
        p[:mW] = chunk[i,2]
        p[:eP] = chunk[i,3]
        p[:mUs] = NaN
        println(p)
        p0 = copy(p)
        ###The tag was replaced from "io" to "outfile1" with the intention of creating a "report card" output file down the line, which would necessitate the existence of multiple output files
        ###This generates the headers for our data output, we're preparing the file beforehand
        ###This creates our "steps" for the data production, as stated in the .*_Pert_.*.jl file
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
            end
        catch
        end

        writedlm( io, [vcat(i, CoRa)],'\t');
        flush(io)
        println(string("Line ", i, " done! sim ", piece))
    end
    end
end
end
end
