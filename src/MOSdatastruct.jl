# ========================================< tSolution >========================================

# type corresponding to a solution
mutable struct tSolution{T}
    x :: Vector{T}                # vector variables x (1..n)
    y :: Vector{T}                # vector outcomes  y (1..p)
end

# ====================< Operator >====================

@inline function Base.:<(s1::tSolution{T}, s2::tSolution{T}) where T
    n::Int64 = length(s1.y)
    return ones(Int64, n)' * (s1.y .< s2.y) == n # 
end

function Base.:≤(s1::tSolution{T}, s2::tSolution{T}) where T
    n::Int64 = length(s1.y)
    return ones(Int64, n)' * (s1.y .≤ s2.y) == n
end

Base.:(==)(s1::tSolution{T}, s2::tSolution{T}) where T = return (s1.y == s2.y) && (s1.x == s2.x)

@inline isless_lex(s1::tSolution{T}, s2::tSolution{T}) where T = return isless(s1.y, s2.y)

@inline equiv(s1::tSolution{T}, s2::tSolution{T}) where T = return (s1.y == s2.y) && (s1.x != s2.x)

# ====================< Display >====================

Base.show(io::IO, s::tSolution{T}) where T = print("{Solution $(T): x = 1 in $(findall(x -> x == 1, s.x)) and y = $(s.y)}")



# ========================================< tChainlist >========================================
# Note 1: the dominated solution removal only work for bi-objective solutions.

global verbose_tChainList = false # debug purpose 

mutable struct tNode{T}
    value::tSolution{T}
    next::Union{tNode, Nothing}
end

mutable struct tChainList{T}
    head    ::Union{tNode{T}, Nothing}
    length  ::Int64
    const del_dominated ::Bool
end

tChainList(t::Type = Float64, b::Bool = true) = tChainList{t}(nothing, 0, b)

# add a solution to the list and remove any dominated solution from the list
function add!(list::tChainList{T}, new_sol::tSolution{T}) where T
    new_node = tNode(new_sol, nothing)
    
    current::Union{tNode{T}, Nothing} = list.head

    if list.head === nothing
        # Special case: list was empty
        verbose_tChainList && println("Create head: $(new_sol.y)")

        list.head = new_node
        list.length = 1

        return list

    elseif isless_lex(new_sol, list.head.value)
        # Special case: new value is smaller than the head
        verbose_tChainList && println("Add head: $(new_sol.y)")

        new_node.next = list.head
        list.head = new_node

        list.length += 1
    else

        # Traverse the list to find the correct position to insert
        while (current.next !== nothing) && isless_lex(current.next.value, new_sol)
            if list.del_dominated && ((current.value < new_sol) || (current.value == new_sol))
                verbose_tChainList && println("Not adding: $(new_sol.y) dominated by $(current.value.y)")
                return list
            end

            current = current.next
        end

        if list.del_dominated && ((current.value < new_sol) || (current.value == new_sol))
            verbose_tChainList && println("Not adding: $(new_sol.y) dominated by $(current.value.y)")
            return list
        end

        # Insert the new node
        verbose_tChainList && println("Add middle: $(new_sol.y)")
        new_node.next = current.next
        current.next = new_node
        list.length += 1
    end

    if list.del_dominated
        current = new_node
        while current.next !== nothing && new_sol < current.next.value
            verbose_tChainList && println("Delete: $(current.next.value.y) dominated by $(new_sol.y)")
            current.next = current.next.next
            list.length -= 1
        end
    end

    return list
end

function Base.show(io::IO, list::tChainList{T}) where T
    print("tChainlist{Int64}(length = $(list.length), ")
    if list.head !== nothing
        current::Union{tNode{T}, Nothing} = list.head
        while current != nothing
            print("{Sol: $(current.value.y)}$((current.next != nothing) ? (" -> ") : (" ->| "))")
            current = current.next
        end        
    end
    print(")")
end

rand_tSolution(lx::Int64=10, ly::Int64 = 2) = return tSolution{Int64}(rand(0:1, lx), rand(0:9, ly))

function to_array(list::tChainList{T}) where T
    i       ::Int64                     = 1                                 # storing current solution index
    current ::Union{tNode{T}, Nothing}  = list.head                         # current solution
    nb_var  ::Int64                     = (current == nothing) ? 0 : length(current.value.x)
    var     ::Array{T}                  = Array{T}(undef, nb_var, list.length)        # gather variable values "x"
    obj     ::Array{T}                  = Array{T}(undef, 2, list.length)   # gather objective values "z(x)"

    println("size -> var = $(size(var)), obj = $(size(var))")

    while current != nothing
        var[:, i] = current.value.x     # store current solution variable value
        obj[:, i] = current.value.y     # store current solution objective value 

        i += 1                          # move to next index
        current = current.next          # move to next solution
    end

    return var, obj
end

function count_equiv(list::tChainList{T}) where T
    current     ::Union{tNode{T}, Nothing}  = list.head                         # current solution
    previous    ::Union{tNode{T}, Nothing}  = nothing                           # previous solution
    nb_equiv    ::Int64 = 0

    while current != nothing
        if (!(previous === nothing) && (current.value.y == previous.value.y))
            nb_equiv += 1
        end 

        previous = current
        current = current.next          # move to next solution
    end

    return nb_equiv
end

# ========================================< tGravityWay >========================================
# type of gravity used 

@enum tGravityWay begin
    ONE     = 1 # wheight of the gravity machine will be set to one
    CONE    = 2 # wheight will be influenced by given cone
    SPA     = 3 # wheight will be set according to some SPA specific criterion
    WSUM    = 4 # wheight will be set as a wheighted sum of the objective values 
    RAND    = 5 # wheaght will be set randomly 
end
