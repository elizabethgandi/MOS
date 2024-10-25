# ==============================================================================
# Compute the performance on 2 objectives of a solution x

function examineVectorVariables(x)

    nbvar = length(x)
    nbUns = 0
    nbFrac = 0
    for j in 1:nbvar
        v = value(x[j])
        if round(v,digits=6) == 1.0
            nbUns += 1
        end
        if round(v,digits=6) != 0.0 && round(v,digits=6) != 1.0
            nbFrac += 1
        end
    end
    return nbUns,nbFrac
end

function examineVectorVariables2(x)

    nbvar = length(x)
    nbUns = 0
    nbFrac = 0

    totVariablesFrac::Vector{Float64} = []
    idf = []

    for j in 1:nbvar
        v = value(x[j])
        if round(v,digits=6) == 1.0
            nbUns += 1
        end
        if round(v,digits=6) != 0.0 && round(v,digits=6) != 1.0
            nbFrac += 1
            push!(totVariablesFrac, v)
            push!(idf, j)
        end
    end

    @show totVariablesFrac

    return nbUns,nbFrac, idf
end
