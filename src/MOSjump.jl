
# ==============================================================================
# Create (set) a 2SPA model (bi-objective partionning problem) with JuMP/MOA 
# from a instance file,  with all variables :Bin or :Con

function set2SPA(C::Array{Int,2}, A::Array{Int,2}, varType::Symbol)

    nbctr, nbvar = size(A)
    m2SPA = Model()
    if varType == :Bin
        @variable(m2SPA, x[1:nbvar], Bin)
    elseif varType == :Con
        @variable(m2SPA, 0.0 <= x[1:nbvar] <= 1.0 )
    else
        @assert false "error: unknow type to set for the variables"        
    end
    @expression(m2SPA, obj1, sum((C[1,i])*x[i] for i in 1:nbvar))
    @expression(m2SPA, obj2, sum((C[2,i])*x[i] for i in 1:nbvar))
    @objective(m2SPA, Min, [obj1, obj2])
    @constraint(m2SPA, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1) 

    return m2SPA
end

# fixed 0s
function set2SPA_2(sols, C::Array{Int,2}, A::Array{Int,2})
    res = []
    obj = []

    n::Int64 = length(sols)

    println("n = $n")

    for s = 1:n
        not_used::Vector{Int64} = findall(x -> x == 0, sols[s])
        used::Vector{Int64} = findall(x -> x == 1, sols[s])
        frac::Vector{Int64} = findall(x -> 0 < x < 1, sols[s])

        frac = [frac; used]

        if length(frac) > 0

            nbctr, nbvar = size(A)
            m2SPA_2 = Model()
            
            @variable(m2SPA_2, x[1:nbvar], Bin)
            
            @expression(m2SPA_2, obj1, sum((C[1,i])*x[i] for i in 1:nbvar))
            @expression(m2SPA_2, obj2, sum((C[2,i])*x[i] for i in 1:nbvar))
            @objective(m2SPA_2, Min, [obj1, obj2])
            @constraint(m2SPA_2, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1) 

            for i in not_used
                @constraint(m2SPA_2, x[i] == 0)
            end

            set_silent(m2SPA_2)
            set_optimizer(m2SPA_2, () -> MOA.Optimizer(Gurobi.Optimizer))
            set_optimizer_attribute(m2SPA_2, MOA.Algorithm(), MOA.Dichotomy())
            set_optimizer_attribute(m2SPA_2, MOA.SolutionLimit(), 10)

            optimize!(m2SPA_2)

            nb_sol = result_count(m2SPA_2)
            
            println("STATUS = $(termination_status(m2SPA_2)), sol count = $(nb_sol)")
            
            for i=1:nb_sol
                tmp = zeros(Int64, nbvar)
                for (k, e) in enumerate(value.(m2SPA_2[:x]; result = i))
                    if e >= 0.7
                        tmp[k] = 1
                    end
                end
                # for e in used
                #     tmp[e] = 1
                # end
                push!(res, tmp)
                push!(obj, (sum([tmp[i] * C[1, i] for i=1:nbvar]), sum([tmp[i] * C[2, i] for i=1:nbvar])))
            end
        end
    end

    return res, obj
end


# For each generator floating point compute an integer solution.
# change the objective to minimize the number of variable with different value of "x"
function set2SPA_penality(generator::tChainList{Float64}, C::Array{Int,2}, A::Array{Int,2}, obj::tGravityWay = ONE, env::Gurobi.Env = Gurobi.Env())
    res     ::tChainList{Int64} = tChainList(Int64, true) # list of integer solution (keep all solution in lexicographic order)
    opti    ::tChainList{Int64} = tChainList(Int64, false) # list of integer solution (keep all solution in lexicographic order)

    current::Union{tNode{Float64}, Nothing} = generator.head

    max_C1 = maximum(C[1, :])
    max_C2 = maximum(C[2, :])
    while current !== nothing

        not_used    ::Vector{Int64} = findall(x -> x == 0, current.value.x)
        used        ::Vector{Int64} = findall(x -> x == 1, current.value.x)
        frac        ::Vector{Int64} = findall(x -> 0 < x < 1, current.value.x)
        penalized   ::Vector{Int64} = [not_used; used]

        if length(frac) > 0

            nbctr, nbvar = size(A)
            m2SPA_2 = Model(() -> Gurobi.Optimizer(env))
            
            @variable(m2SPA_2, x[1:nbvar], Bin)
            @variable(m2SPA_2, y[penalized], Bin)

            if obj == ONE
                @objective(m2SPA_2, Min, sum([y[i] for i in penalized]))
            elseif obj == CONE
                # TODO
            elseif obj == SPA
                @objective(m2SPA_2, Min, sum([y[i] * (((C[1, i] + C[2, i])/(max_C1 + max_C2)) + sum(A[:, i])) for i in penalized]))
            elseif obj == WSUM
                @objective(m2SPA_2, Min, sum([y[i] * ((C[1, i] + C[2, i])/(max_C1 + max_C2)) for i in penalized]))
            end

            # SPA assignment constraints
            @constraint(m2SPA_2, [i=1:nbctr], (sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1) 

            # Penality constraints
            for i in not_used
                @constraint(m2SPA_2, x[i] ≤ y[i])
            end

            for i in used
                @constraint(m2SPA_2, x[i] ≥ 1 - y[i])
            end

            set_silent(m2SPA_2)

            optimize!(m2SPA_2)

            opti_x = value.(m2SPA_2[:x])
            opti_y = [sum([opti_x[i] * C[1, i] for i=1:nbvar]), sum([opti_x[i] * C[2, i] for i=1:nbvar])]
            add!(opti, tSolution{Int64}(opti_x, opti_y))

            nb_sol = result_count(m2SPA_2)
            
            println("    → resolution status = $(termination_status(m2SPA_2)), feasible sol count = $(nb_sol)")
            
            for i=1:nb_sol
                tmp = zeros(Int64, nbvar)
                for (k, e) in enumerate(value.(m2SPA_2[:x]; result = i))
                    if e >= 0.7
                        tmp[k] = 1
                    end
                end
                add!(res, tSolution{Int64}(tmp, [sum([tmp[i] * C[1, i] for i=1:nbvar]), sum([tmp[i] * C[2, i] for i=1:nbvar])]))
            end
        end
        current = current.next
    end

    return res, opti
end

function set2SPA_penality_fixedzeros(generator::tChainList{Float64}, C::Array{Int,2}, A::Array{Int,2}, obj::tGravityWay = ONE, env::Gurobi.Env = Gurobi.Env())
    res::tChainList{Int64} = tChainList(Int64, true) # list of integer solution (keep all solution in lexicographic order)

    current::Union{tNode{Float64}, Nothing} = generator.head

    max_C1 = maximum(C[1, :])
    max_C2 = maximum(C[2, :])
    while current !== nothing

        not_used    ::Vector{Int64} = findall(x -> x == 0, current.value.x)
        used        ::Vector{Int64} = findall(x -> x == 1, current.value.x)
        frac        ::Vector{Int64} = findall(x -> 0 < x < 1, current.value.x)
        remaining_var ::Vector{Int64} = [used; frac]

        if length(frac) > 0

            nbctr, nbvar = size(A)
            m2SPA_2 = Model(() -> Gurobi.Optimizer(env))
            
            @variable(m2SPA_2, x[remaining_var], Bin)
            @variable(m2SPA_2, y[used], Bin)

            if obj == ONE
                @objective(m2SPA_2, Min, sum([y[i] for i in used]))
            elseif obj == CONE
                # TODO
            elseif obj == SPA
                @objective(m2SPA_2, Min, sum([y[i] * (((C[1, i] + C[2, i])/(max_C1 + max_C2)) + sum(A[:, i])) for i in used]))
            elseif obj == WSUM
                @objective(m2SPA_2, Min, sum([y[i] * ((C[1, i] + C[2, i])/(max_C1 + max_C2)) for i in used]))
            end

            # SPA assignment constraints
            @constraint(m2SPA_2, [i=1:nbctr], (sum((x[j]*A[i,j]) for j in remaining_var)) == 1) 

            for i in used
                @constraint(m2SPA_2, x[i] ≥ 1 - y[i])
            end

            set_silent(m2SPA_2)

            optimize!(m2SPA_2)

            nb_sol = result_count(m2SPA_2)
            
            println("    → resolution status = $(termination_status(m2SPA_2)), feasible sol count = $(nb_sol)")
            
            (nb_sol ≥ 1) || (current = current.next; continue)

            for i=1:nb_sol
                tmp = zeros(Int64, nbvar)
                for (k, e) in enumerate(value.(m2SPA_2[:x]; result = i))
                    if e >= 0.7
                        tmp[remaining_var[k]] = 1
                    end
                end
                add!(res, tSolution{Int64}(tmp, [sum([tmp[i] * C[1, i] for i=1:nbvar]), sum([tmp[i] * C[2, i] for i=1:nbvar])]))
            end
        end
        current = current.next
    end

    return res
end

function set2SPA_penality_fixedones(generator::tChainList{Float64}, C::Array{Int,2}, A::Array{Int,2}, obj::tGravityWay = ONE, env::Gurobi.Env = Gurobi.Env())
    res::tChainList{Int64} = tChainList(Int64, true) # list of integer solution (keep all solution in lexicographic order)

    current::Union{tNode{Float64}, Nothing} = generator.head

    max_C1 = maximum(C[1, :])
    max_C2 = maximum(C[2, :])
    while current !== nothing

        not_used    ::Vector{Int64} = findall(x -> x == 0, current.value.x)
        used        ::Vector{Int64} = findall(x -> x == 1, current.value.x)
        frac        ::Vector{Int64} = findall(x -> 0 < x < 1, current.value.x)
        remaining_var ::Vector{Int64} = [not_used; frac]

        if length(frac) > 0

            nbctr, nbvar = size(A)
            m2SPA_2 = Model(() -> Gurobi.Optimizer(env))
            
            @variable(m2SPA_2, x[remaining_var], Bin)
            @variable(m2SPA_2, y[not_used], Bin)

            if obj == ONE
                @objective(m2SPA_2, Min, sum([y[i] for i in not_used]))
            elseif obj == CONE
                # TODO
            elseif obj == SPA
                @objective(m2SPA_2, Min, sum([y[i] * (((C[1, i] + C[2, i])/(max_C1 + max_C2)) + sum(A[:, i])) for i in not_used]))
            elseif obj == WSUM
                @objective(m2SPA_2, Min, sum([y[i] * ((C[1, i] + C[2, i])/(max_C1 + max_C2)) for i in not_used]))
            end

            # SPA assignment constraints
            for i in 1:nbctr
                if sum([A[i,j] for j in used]) ≥ 1 # constraint already active → set all remaining var to zeros
                    for j in not_used
                        @constraint(m2SPA_2, sum((x[j] * A[i,j])) == 0)
                    end
                else # insactive constraint define as in the original problem
                    @constraint(m2SPA_2, (sum((x[j] * A[i,j]) for j in remaining_var)) == 1)
                end
            end 

            for i in not_used
                @constraint(m2SPA_2, x[i] ≤ y[i])
            end

            set_silent(m2SPA_2)

            optimize!(m2SPA_2)

            nb_sol = result_count(m2SPA_2)
            
            println("    → resolution status = $(termination_status(m2SPA_2)), feasible sol count = $(nb_sol)")
            
            (nb_sol ≥ 1) || (current = current.next; continue)

            for i=1:nb_sol
                tmp = zeros(Int64, nbvar)
                for (k, e) in enumerate(value.(m2SPA_2[:x]; result = i))
                    if e >= 0.7
                        tmp[remaining_var[k]] = 1
                    end
                end
                for e in used
                    tmp[e] = 1
                end
                add!(res, tSolution{Int64}(tmp, [sum([tmp[i] * C[1, i] for i=1:nbvar]), sum([tmp[i] * C[2, i] for i=1:nbvar])]))
            end
        end
        current = current.next
    end

    return res
end


# ==============================================================================

# ==============================================================================
# Create (load) a 2SPA model (bi-objective partionning problem) with JuMP/MOA 
# from a instance file,  with all variables set to :Bin

function load2SPA(fname::String)

    C, A = parse2SPA(fname)
    m2SPA = set2SPA(C, A, :Bin)

    return m2SPA
end
