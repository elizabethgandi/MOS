# ==============================================================================
# ==============================================================================

println("-"^80)
println("""\nElizabeth et Eliot \n""")
println("""\nGravityMachine : etude des generateurs \n""")

global verbose    = true
global graphic    = true
global exact      = true
global experiment = false

print("  verbose.....: "); verbose    ? println("yes") : println("no") 
print("  graphics....: "); graphic    ? println("yes") : println("no") 
print("  exact.......: "); exact      ? println("yes") : println("no") 
print("  experiment..: "); experiment ? println("yes") : println("no") 

println("\n\n-) Activation des packages necessaires\n")
using JuMP, GLPK, Gurobi
import MultiObjectiveAlgorithms as MOA
using Printf, PyPlot
println("  Fait\n")

include("MOSfiles.jl")
include("MOSjump.jl")
include("MOSsolve.jl")
include("MOStools.jl")


# ==============================================================================
# point d'entree principal

function main(fname::String)


    # --------------------------------------------------------------------------
    # Charge une instance numerique   

    println("\n0) instance et dimensions \n")
    println("  instance = $fname") 
   
    m2SPA = load2SPA(fname)
    println("  nbvar    = ",num_variables(m2SPA))
    println("  nbctr    = ",num_constraints(m2SPA, AffExpr, MOI.EqualTo{Float64}),"\n\n")  


    # --------------------------------------------------------------------------
    # Calcul de Y_N et Y_SN
 
    if exact
        println("1) calcule Y_N avec methode ϵ-constraint et x∈{0,1}")
        YN, r = solve2SPA(m2SPA, :Gurobi, :EpsilonConstraint, :Bin)
        sizeYN = size(YN,2)

        println("2) calcule Y_SN avec methode dichotomique et x∈{0,1}")
        YSN, r = solve2SPA(m2SPA, :Gurobi, :Dichotomy, :Bin)
        sizeYSN = size(YSN,2)
    end


    # --------------------------------------------------------------------------
    # Calcul de L pour Y_N (cad les generateurs selon 2 methodes de calcul) 

    println("3) calcule LB(Y_N) avec methode ϵ-constraint et 0≤x≤1")
    nbProbe = 16
    LBE, e = solve2SPA(m2SPA, :Gurobi, :EpsilonConstraint, :Con, nbPoints=nbProbe)

    println("4) calcule LB(Y_N) avec methode dichotomique et 0≤x≤1")
    nbProbe = 16
    LBD, vect_dicho = solve2SPA(m2SPA, :Gurobi, :Dichotomy, :Con, nbPoints=nbProbe)  

    # --------------------------------------------------------------------------
    # Avancées du 4/10:
    #@show LBD
    #@show vect_dicho
    #vect_fract::Vector{Int64} = findall(x->x>1,vect_dicho)
    #vect_1::Vector{Int64} = findall(x->x==1,vect_dicho)
    #vect_0::Vector{Int64} = zeros(nbVar)

    #@show vect_fract
    #@show vect_1
    #@show vect_0

    # fixer les variables à 0 et 1 et lancer une recherche exacte


    # --------------------------------------------------------------------------
    # Sortie graphique

    if graphic

        fig1 = figure("Objective Space Y",figsize=(6.5,5))
        title("single ϵ-constraint vs dichotomy | max " * string(nbProbe) * " samples")
        xlabel(L"z^1(x)")
        ylabel(L"z^2(x)")

        if exact
            # YN : all non-dominated points
            scatter(YN[1,:], YN[2,:], c="lime", s=50, label=L"$Y_N$")
            # YSN : all non-dominated supported points
            plot(YSN[1,:], YSN[2,:], c="green", mec="lime", marker="o", linestyle="dotted", label=L"$Y_{SN}$", markersize=7) 
        end

        # LBE : generateurs obtenus avec une ϵ-constraint simple avec un pas predefini (-> echantillonnage)
        scatter(LBE[1,:], LBE[2,:], c="red", marker="x", s=80, label=L"$LB$ eps") 
        # LBD : generteurs obtenus avec une dichotomie
        plot(LBD[1,:], LBD[2,:], c="blue", marker="o", linestyle="dotted", label=L"$LB$ dic", markersize=5) 

        legend() 

    end

    return nothing

end


# ==============================================================================
# Run all algorithms over all instances

function numericalExperiment(target)

    global graphic = false
    global exact   = false # some instances are very long to be solved with glpk
    
    fnames = getfname(target)
    for instance = 1:length(fnames)
        start = time()
        main( string(target,"/",fnames[instance]) )    
        elapsedTime = time()-start
        println("Elapsed time for $(fnames[instance]) : $elapsedTime (s)")
    end

    return nothing
end


# ==============================================================================

target = "../SPA/instances" # path for a standard config on macOS

if experiment
    numericalExperiment(target)
else
    #@time main(target*"/bio"*"sppaa02.txt")
    #@time main(target*"/bio"*"sppnw03.txt")
    #@time main(target*"/bio"*"sppnw04.txt")
    @time main(target*"/bio"*"sppnw10.txt")
    #@time main(target*"/bio"*"sppnw20.txt")
    #@time main(target*"/bio"*"sppnw25.txt")
    #@time main(target*"/bio"*"didactic3.txt")
    #@time main(target*"/bio"*"didactic5.txt")
    #@time main(target*"/bio"*"sppnw29.txt")
    #@time main(target*"/bio"*"sppnw19.txt")
    #@time main(target*"/bio"*"sppnw40.txt")
end

nothing
