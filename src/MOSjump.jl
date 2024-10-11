
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

# penality for zros that changes to ones and one that changes to zeros
function set2SPA_3(sols, C::Array{Int,2}, A::Array{Int,2}, env::Gurobi.Env = Gurobi.Env())
    res = []
    obj = []

    _, n::Int64 = size(sols)

    # println("n = $n")

    for s = 1:n
        not_used::Vector{Int64} = findall(x -> x == 0, sols[:, s])
        used::Vector{Int64} = findall(x -> x == 1, sols[:, s])
        frac::Vector{Int64} = findall(x -> 0 < x < 1, sols[:, s])

        penalized::Vector{Int64} = [not_used; used]

        frac = [frac; used]

        if length(frac) > 0

            nbctr, nbvar = size(A)
            m2SPA_2 = Model(() -> Gurobi.Optimizer(env))
            
            @variable(m2SPA_2, x[1:nbvar], Bin)
            @variable(m2SPA_2, y[1:nbvar], Bin)
            
            @objective(m2SPA_2, Min, sum(y))

            @constraint(m2SPA_2, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1) 

            for i in not_used
                @constraint(m2SPA_2, x[i] ≤ y[i])
            end

            for i in used
                @constraint(m2SPA_2, x[i] ≥ 1 - y[i])
            end

            set_silent(m2SPA_2)

            optimize!(m2SPA_2)

            nb_sol = result_count(m2SPA_2)
            
            println("    → resolution status = $(termination_status(m2SPA_2)), feasible sol count = $(nb_sol)")
            
            for i=1:nb_sol
                tmp = zeros(Int64, nbvar)
                for (k, e) in enumerate(value.(m2SPA_2[:x]; result = i))
                    if e >= 0.7
                        tmp[k] = 1
                    end
                end
                push!(res, tmp)
                push!(obj, (sum([tmp[i] * C[1, i] for i=1:nbvar]), sum([tmp[i] * C[2, i] for i=1:nbvar])))
            end
        end
    end

    return res, obj
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
