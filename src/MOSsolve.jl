
# ==============================================================================
# compute S_N for a 2SPA (JuMP/MOA model) for methods and options selected

function solve2SPA(
    m2SPA::Model,          # a 2SPA model
    solverMIP::Symbol,     # MIP solver to use (:GLPK or :Gurobi)
    methodMOA::Symbol,     # MOA method to use (:EpsilonConstraint or :Dichotomy or :Lexicographic)
    varType::Symbol;       # nature of variables (:Bin or :Con)
    nbPoints::Int64=0      # number of points to compute (optional parameter with default value)
    )


    # ---- precondition
    if methodMOA==:EpsilonConstraint && varType==:Con && nbPoints==0
        @assert false "warning: no reasonable to probe all YN with ϵ-constraint and 0≤x≤1"        
    end        


    # ---- Get the dimensions of the 2SPA instance to solve
    nbvar = num_variables(m2SPA)
    nbctr = num_constraints(m2SPA, AffExpr, MOI.EqualTo{Float64})
    start = time()


    # ---- Relax the integrality constraints on variables if varType == :Con
    if varType == :Con
        undo_relax = relax_integrality(m2SPA)
    end


    # ---- Setting the solver
    if solverMIP == :GLPK
        set_optimizer(m2SPA, () -> MOA.Optimizer(GLPK.Optimizer))
    elseif solverMIP == :Gurobi
        set_optimizer(m2SPA, () -> MOA.Optimizer(Gurobi.Optimizer))
    else
        @assert false "error: unavailable MIP solver requested"
    end   

    set_silent(m2SPA)

    if methodMOA == :EpsilonConstraint
        set_optimizer_attribute(m2SPA, MOA.Algorithm(), MOA.EpsilonConstraint())
    elseif methodMOA == :Dichotomy
        set_optimizer_attribute(m2SPA, MOA.Algorithm(), MOA.Dichotomy())
    elseif methodMOA == :Lexicographic
        set_optimizer_attribute(m2SPA, MOA.Algorithm(), MOA.Lexicographic())
    else
        @assert false "error: unavailable MOA method requested"
    end

    if nbPoints > 0
        set_optimizer_attribute(m2SPA, MOA.SolutionLimit(), nbPoints)
    end


    # ---- Run the solver
    optimize!(m2SPA)

    # ---- Querying the results
    cardSN::Int64 = result_count(m2SPA)
    verbose ? println("  cardSN = $cardSN") : nothing

    SN = Array{Number}(undef,2,cardSN)
    fsol::Matrix{Number}        = zeros(Number,nbvar,cardSN)
    SX::Vector{Vector{Float64}} = []
    sumNbFrac::Float64          = 0.0

    for i in 1:cardSN

        if varType == :Bin
            SN[1,i] = round(Int64,value(m2SPA[:obj1]; result = i))
            SN[2,i] = round(Int64,value(m2SPA[:obj2]; result = i))
            verbose ? @printf("  %3d: z=[%6d,%6d] | ", i, SN[1,i], SN[2,i]) : nothing
        else
            SN[1,i] = value(m2SPA[:obj1]; result = i)
            SN[2,i] = value(m2SPA[:obj2]; result = i)
            push!(SX, value.(m2SPA[:x]; result = i))
            verbose ? @printf("  %3d: z=[%9.2f,%9.2f] | ", i, SN[1,i], SN[2,i]) : nothing
        end

        nbUns,nbFrac = examineVectorVariables( value.(m2SPA[:x]; result = i) )
        sumNbFrac += nbFrac
        verbose ? print("nbOne=$nbUns  nbFrac=$nbFrac  \n") : nothing

    end

    verbose ? println("  nbVar = $nbvar  moyNbFrac = ",sumNbFrac/cardSN, " ⇒ ", round(100*sumNbFrac/cardSN/nbvar, digits=2),"%") : nothing

    nf::Int64              = 0
    n1::Int64              = 0
    nf_vect::Vector{Int64} = Vector{Int64}()
    x::Vector{Float64}     = Vector{Float64}()

    @printf("   i:   #0   #1   #f \n")
    for y in 1:cardSN
        n0 = 0; n1 = 0; nf = 0
        x =  value.(m2SPA[:x]; result = y)

        for j in 1:length(x)
            if     x[j] <= 0.1
                n0 +=1
            elseif x[j] >= 0.9
                n1 +=1
            else
                nf +=1
            end
        end
        @printf(" %3d: %4d %4d %4d \n", y, n0, n1, nf)
        
        if (nf != 0)
            push!(nf_vect,nf)
            fsol[:, y] =  x
        end
    end

    # ---- Restore the integrality constraints on variables if varType == :Con    
    if varType == :Con
        undo_relax()
    end

    elapsedTime::Float64 = time()-start
    println("  Elapsed time: $(round(elapsedTime,digits=3))s \n\n ")

    return SN, nf_vect, fsol, cardSN, SX
end


function fonction_deux_resolutions(cardSN::Int64, m2SPA::Model, C::Array{Int,2}, A::Array{Int,2}, SX::Vector{Vector{Float64}})
    
    nbvar::Int64 = num_variables(m2SPA)
    nbctr::Int64 = num_constraints(m2SPA, AffExpr, MOI.EqualTo{Float64})

    tot_obj_value1::Vector{Float64} = Vector{Float64}()
    tot_obj_value2::Vector{Float64} = Vector{Float64}() 
    
    λ::Float64        = rand()
    obj_value1::Int64 = 0
    obj_value2::Int64 = 0

    for i in 1:cardSN  

        obj_value1 = 0
        obj_value2 = 0

        println("\ni : $i")
        
        _, nbFrac, idf = examineVectorVariables2(SX[i]) 
        
        if nbFrac == 0
            println("Réalisable")

            obj_value1 = sum(C[1, j] * SX[i][j] for j in 1:nbvar)  
            obj_value2 = sum(C[2, j] * SX[i][j] for j in 1:nbvar)

            println("obj Val1: $obj_value1 || obj Val2: $obj_value2 ")

            push!(tot_obj_value1, obj_value1)
            push!(tot_obj_value2, obj_value2)

        else
            println("Non Réalisable")

            m2SPA_2 = Model()
            @variable(m2SPA_2, 0.0 <= xPrim[1:nbvar] <= 1.0)
            @expression(m2SPA_2, obj1λ, sum(λ * C[1, j] * xPrim[j] + (1.0 - λ) * (C[2, j]) * xPrim[j] for j in 1:nbvar))
            @objective(m2SPA_2, Min, obj1λ)
            @constraint(m2SPA_2, [h = 1:nbctr], sum(xPrim[j] * A[h, j] for j in 1:nbvar) == 1)
            
            for k in idf
                set_binary(xPrim[k])
            end
            
            set_optimizer(m2SPA_2, Gurobi.Optimizer)
            set_silent(m2SPA_2)
            optimize!(m2SPA_2)
            
            if termination_status(m2SPA_2) == OPTIMAL
                _ , nbFrac, _ = examineVectorVariables2(value.(xPrim))
                
                if nbFrac == 0
                    println("SS pb réalisable")

                    z1, z2 = evaluerFct(C, nbvar, value.(xPrim))

                    println("z1: $z1 || z2: $z2")

                    push!(tot_obj_value1, z1)
                    push!(tot_obj_value2, z2)
                    
                else
                    println("SS pb non réalisable")
                end
            end
        end
       println("----"^30)
    end

    return  tot_obj_value1, tot_obj_value2
end


function evaluerFct(C::Array{Int,2}, nbvar::Int64, x::Vector{Float64})

    z1::Float64 = 0.0
    z2::Float64 = 0.0

    for i in 1:nbvar
        z1 += C[1, i]*x[i]
        z2 += C[2, i]*x[i]
    end

    return z1, z2
end