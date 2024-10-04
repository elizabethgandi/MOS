
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

function set2SPA_2(m2SPA, C::Array{Int,2}, A::Array{Int,2})
    res = []

    n::Int64 = result_count(m2SPA)
    relax = []
    for i=1:n
        push!(relax, value.(m2SPA[:x]; result = i))
    end

    for s = 1:n
        used::Vector{Int64} = findall(x -> x == 1, relax[s])
        frac::Vector{Int64} = findall(x -> 0 < x < 1, relax[s])

        if length(frac) > 0

            nbctr, nbvar = size(A)
            m2SPA_2 = Model()
            
            @variable(m2SPA_2, x[used], Bin)
            
            @expression(m2SPA_2, obj1, sum((C[1,i])*x[i] for i in frac))
            @expression(m2SPA_2, obj2, sum((C[2,i])*x[i] for i in frac))
            @objective(m2SPA_2, Min, [obj1, obj2])

            for i=1:nbctr
                if !(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1
                    @constraint(m2SPA_2, (sum((x[j]*A[i,j]) for j in frac)) == 1)
                end
            end

            optimize!(m2SPA_2)
            
            for i=1:result_count(m2SPA_2)
                tmp = zeros(Int64, nbvar)
                for (k, e) in enumerate(value.(m2SPA_2[:x]; result = i))
                    if e >= 0.7
                        tmp[k] = 1
                    end
                end
                for e in used
                    tmp[e] = 1
                end
                push!(res, tmp)
            end
        end
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
