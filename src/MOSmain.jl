# ==============================================================================
# ==============================================================================

println("-"^80)
println("""\nElizabeth et Elliot \n""")
println("""\nGravityMachine : etude des generateurs \n""")

global verbose           = true
global graphic           = true
global exact             = true
global experiment        = false
global dichotomique      = true
global deux_resolutions  = true
global penalite_ponderee = true

print("  verbose............: "); verbose           ? println("yes") : println("no") 
print("  graphics...........: "); graphic           ? println("yes") : println("no") 
print("  exact..............: "); exact             ? println("yes") : println("no") 
print("  experiment.........: "); experiment        ? println("yes") : println("no")
print("  dichotomique.......: "); dichotomique      ? println("yes") : println("no") 
print("  deux_resolutions...: "); deux_resolutions  ? println("yes") : println("no") 
print("  penalite_ponderee..: "); penalite_ponderee ? println("yes") : println("no") 

println("\n\n-) Activation des packages necessaires\n")
using JuMP, GLPK, Gurobi
import MultiObjectiveAlgorithms as MOA
using Printf, PyPlot
println("  Fait\n")

include("MOSdatastruct.jl")
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

    if (dichotomique == false) 
        println("IMPOSSIBILITE: On a besoin de la méthode dichotomique")

        return nothing
    end

    # --------------------------------------------------------------------------
    # Calcul de Y_N et Y_SN
 
    if exact
        println("1) calcule Y_N avec methode ϵ-constraint et x∈{0,1}")
        YN, r, mdYn = solve2SPA(m2SPA, :Gurobi, :EpsilonConstraint, :Bin)
        sizeYN = size(YN,2)

        println("2) calcule Y_SN avec methode dichotomique et x∈{0,1}")
        YSN, r, mdYsn = solve2SPA(m2SPA, :Gurobi, :Dichotomy, :Bin)
        sizeYSN = size(YSN,2)
    end


    # --------------------------------------------------------------------------
    # Calcul de L pour Y_N (cad les generateurs selon 2 methodes de calcul) 

    println("3) calcule LB(Y_N) avec methode ϵ-constraint et 0≤x≤1")
    nbProbe = 16
    LBE, fvar_ϵ, fsol_ϵ = solve2SPA(m2SPA, :Gurobi, :EpsilonConstraint, :Con, nbPoints=nbProbe)

    println("4) calcule LB(Y_N) avec methode dichotomique et 0≤x≤1")
    nbProbe = 16
    LBD, fvar_dico, fsol_dico = solve2SPA(m2SPA, :Gurobi, :Dichotomy, :Con, nbPoints=nbProbe)  

    # --------------------------------------------------------------------------

    # PISTE 2: FIXATION DES 0 À 0 ET RESOLUTION EXACTE SUR LES VALEURS À FRAC ET 1
    if deux_resolutions
        println("\nPISTE 2: deux_resolutions ----------------------------------------------------\n")
    end


    # PISTE 3: PÉNALISATION DES 0 ET 1 
    if penalite_ponderee
        println("\nPISTE 3: Pénalité pondérée ---------------------------------------------------\n")

        generator::tChainList{Float64} = tChainList(Float64)
        _, nbSol_LBE = size(LBE)
        # println("ϵ ->  $(fsol_ϵ)\n$(LBE)")
        for i=1:nbSol_LBE
            add!(generator, tSolution{Float64}(fsol_ϵ[:, i], LBE[:, i]))
        end
        _, nbSol_LBD = size(LBD)
        # println("ϵ ->  $(fsol_dico)\n$(LBD)")
        for i=1:nbSol_LBD
            add!(generator, tSolution{Float64}(fsol_dico[:, i], LBD[:, i]))
        end
        C, A = parse2SPA(fname)

        println("\n========================================< Newly computed solution (1 obj penality method) >========================================")
        println("ϵ-constraint and Dichotomy for a total of $(generator.length) generators:")
        new_sol_one = set2SPA_penality(generator, C, A, ONE)
        println("   total new feasible solution (may be dominated, or equal) = $(new_sol_one.length)")
        _, y_ones = to_array(new_sol_one)
        println("Solutions: $(y_ones[1, :])\n           $(y_ones[2, :]) → Note that there is $(count_equiv(new_sol_one)) equivalent solutions.")


        println("\n========================================< Newly computed solution (1 obj penality method) >========================================")
        println("ϵ-constraint and Dichotomy for a total of $(generator.length) generators:")
        new_sol_wsum = set2SPA_penality(generator, C, A, WSUM)
        println("   total new feasible solution (may be dominated, or equal) = $(new_sol_wsum.length)")
        _, y_wsum = to_array(new_sol_wsum)
        println("Solutions: $(y_wsum[1, :])\n           $(y_wsum[2, :]) → Note that there is $(count_equiv(new_sol_wsum)) equivalent solutions.")


        println("\n========================================< Newly computed solution (1 obj penality method) >========================================")
        println("ϵ-constraint and Dichotomy for a total of $(generator.length) generators:")
        new_sol_wspa = set2SPA_penality(generator, C, A, SPA)
        println("   total new feasible solution (may be dominated, or equal) = $(new_sol_wspa.length)")
        _, y_wspa = to_array(new_sol_wspa)
        println("Solutions: $(y_wspa[1, :])\n           $(y_wspa[2, :])\n → Note that there is $(count_equiv(new_sol_wspa)) equivalent solutions.")
    end


    # println("\n========================================< Newly computed solution (1 obj penality method) >========================================")
    # println("ϵ-constraint and Dichotomy for a total of $(generator.length) generators:")
    # new_sol_fixedones = set2SPA_penality_fixedzeros(generator, C, A, SPA)
    # println("   total new feasible solution (may be dominated, or equal) = $(new_sol_fixedones.length)")

    # println("sol = $(length(mdLBD)), vect_dicho = $vect_dicho")

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

        if penalite_ponderee
            scatter(y_ones[1, :], y_ones[2, :], c="orange"  , marker="1", s=120, label="new sol ones")
            scatter(y_wsum[1, :], y_wsum[2, :], c="cyan"    , marker="2", s=120, label="new sol wsum")
            scatter(y_wspa[1, :], y_wspa[2, :], c="purple"  , marker="3", s=120, label="new sol wspa")
        end

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
    #@time main(target*"/bio"*"sppnw03.txt")
    #@time main(target*"/bio"*"sppnw04.txt")
    #@time main(target*"/bio"*"sppnw10.txt")
    #@time main(target*"/bio"*"sppnw20.txt")
    #@time main(target*"/bio"*"sppnw25.txt")
    #@time main(target*"/bio"*"didactic3.txt")
    #@time main(target*"/bio"*"didactic5.txt")
    #@time main(target*"/bio"*"sppnw29.txt")
    #@time main(target*"/bio"*"sppnw19.txt")
    @time main(target*"/bio"*"sppnw40.txt")
end

nothing
