# ==============================================================================
# Compute the performance on 2 objectives of a solution x

function examineVectorVariables(x::Vector{Float64})

    nbvar::Int64  = length(x)
    nbUns::Int64  = 0
    nbFrac::Int64 = 0
    v::Float64    = 0.0

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

function examineVectorVariables2(x::Vector{Float64})

    nbvar::Int64       = length(x)
    nbUns::Int64       = 0
    nbFrac::Int64      = 0
    v::Float64         = 0.0
    idf::Vector{Int64} = Vector{Int64}()
    
    totVariablesFrac::Vector{Float64} = Vector{Float64}()

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

    return nbUns,nbFrac, idf
end
