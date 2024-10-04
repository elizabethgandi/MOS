
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

# ==============================================================================
#=
function set2SPA_2(A::Array{Int64}, C::Array{Int,2}, varType::Symbol)#, xTilde::Array{Int,1})
    
    res::Vector{tSolution{Float64}} = (tSolution{Int64})[]

    n::Int64 = length(L)
    for s = 1:n
        sol::tSolution{Float64} = L[s]
        used::Vector{Int64} = findall(x -> x == 1, sol.x)
        frac::Vector{Int64} = findall(x -> 0 < x < 1, sol.x)

        nbctr = size(A,1)
        nbvar = size(A,2)
        #idxTilde0, idxTilde1 = split01(xTilde)

        proj = Model(Gurobi.Optimizer)
        @variable(proj, 0.0 <= x[frac] <= 1.0 )

        #@expression(proj, obj1, Min, sum(x[i] for i in idxTilde0) )
        #@expression(proj, obj2, Min, sum(x[i] for i in idxTilde1) )

        @expression(m2SPA, obj1, sum((C[1,i])*x[i] for i in 1:nbvar))
        @expression(m2SPA, obj2, sum((C[2,i])*x[i] for i in 1:nbvar))

        for i=1:nbctr
            if !(sum((x[j]A[i,j]) for j in used)) == 1
                @constraint(proj,(sum((x[j]A[i,j]) for j in frac)) == 1)
            end
        end

        optimize!(proj)
        return objective_value(proj), value.(x)
    end
end

=#

function parse_m2SPA_val(model)
    res::Vector{Vector{Float64}} = (Vector{Float64})[]

    for i=1:result_count(model)
        push(res, value.(model[:x]; result = i))
    end

    return res
end

# ==============================================================================
# Create (load) a 2SPA model (bi-objective partionning problem) with JuMP/MOA 
# from a instance file,  with all variables set to :Bin

function load2SPA(fname::String)

    C, A = parse2SPA(fname)
    m2SPA = set2SPA(C, A, :Bin)

    return m2SPA
end
