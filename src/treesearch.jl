if !(@isdefined doimports)
    doimports = false
    println("Import libs")
    include("util.jl")
    include("data_structures/fiboheap.jl")
    println("Finish")
end

using DataStructures
using JuMP
using SCIP
using Dates
using .FibonacciHeaps

mutable struct Automata
  n::Int # Nb states
  ni::Int # Size of the alphabet
  τ #::DefaultDict{Tuple{Int64,Int64}, Array{Any, 1}} # Transition matrix τ[(s, i)] = Output states of s with input i.
  τm1 #::DefaultDict{Tuple{Int64,Int64}, Array{Any, 1}} # Inverted matrix τ^-1[(s, i)] = States for which the output state with input i is s 
  μ # Merging sequences
  μ2 # Merging sequences : μ[z][(x,y)] = μ[z][(y,x)] = Minimum size of the sequence needed to merge x and y into z
  model # Linear program
  consS # Constraints of the linear program
  consZ # Constraints of the linear program
end

Automata() = Automata(0, 0, DefaultDict(() -> []), DefaultDict(() -> []), Dict(), Dict(),Model(() -> SCIP.Optimizer(display_verblevel=0)), nothing, nothing) # Default constructor

"""Build an automata from fsm (struct from util) """
function generateAutomataFSM(fsm)
  α = Automata()
  α.n = fsm.s
  α.ni = fsm.i
  for s in 1:fsm.s
    for i in 1:fsm.i
      push!(α.τ[(s, i)], fsm.succ[s][i]+1)
      push!(α.τm1[(fsm.succ[s][i]+1, i)], s)
    end
  end

  buildMergingSequences(α)
  buildbackwardMergingSequences(α)
  buildPL(α) 

  return α
end

"""Build a deterministic automata with n states and an alphabet of size ni"""
function generateDeterministic(n, ni)
  α = Automata()
  α.n = n
  α.ni = ni
  for s in 1:n
    for i in 1:ni
      s2 = rand(1:n) 
      push!(α.τ[(s, i)], s2)
      push!(α.τm1[(s2, i)], s)
    end
  end

  buildMergingSequences(α)
  buildbackwardMergingSequences(α)
  buildPL(α) 

  return α
end

"""Build the parameter α.μ2"""
function buildbackwardMergingSequences(α)
    function backdij(z)
        α.μ2[z] = Dict()
    
        tovisit = [(z, z, 0)] 
        visited = Set()
        while length(tovisit) != 0
              r1, r2, d = popfirst!(tovisit)

              if r2 < r1
                r1, r2 = r2, r1
              end
              
              if((r1, r2) in visited)
                continue
              end
              push!(visited, (r1, r2))
              α.μ2[z][(r1, r2)] = d
              α.μ2[z][(r2, r1)] = d

              for i in 1:α.ni
                  for s1 in α.τm1[(r1, i)]
                      for s2 in α.τm1[(r2, i)]
                          if (s1, s2) in visited
                              continue
                          end
                          push!(tovisit, (s1, s2, d + 1))
                      end
                  end
              end
        end
    end

    for z in 1:α.n
        backdij(z)
    end


end

"""Build the parameter α.μ"""
function buildMergingSequences(α)
  
  function dij(s1, s2)
    tovisit = [(s1, s2, 0)] 
    visited = Set()
    while length(tovisit) != 0
      r1, r2, d = popfirst!(tovisit)
      
      if r1 == r2
        return d
      end

      if r2 < r1
        r1, r2 = r2, r1
      end
      
      if((r1, r2) in visited)
        continue
      end
      push!(visited, (r1, r2))
      
      for i in 1:α.ni
          push!(tovisit, (α.τ[(r1, i)][1], α.τ[(r2, i)][1], d + 1))
      end
    end

  end

  for s1 in 1:α.n
    for s2 in s1:α.n
      d = dij(s1, s2)
      α.μ[(s1, s2)] = d 
      α.μ[(s2, s1)] = d 
    end
  end

end

"""Check if α has an SS by checking if all the merging sequences exists"""
function hasMergingSequence(α)
  return all(!isnothing(α.μ[(s1, s2)]) for s1 in 1:α.n for s2 in (s1 + 1):α.n)
end

"""Build α.model, α.conS and α.consP"""
function buildPL(α)
  k = 25
  m = α.model
  S = α.n
  @variable(m, 0 <= x[1:k] <= 1)
  @variable(m, 0 <= s[0:S-1, 0:k] <= 1)
  @variable(m, 0 <= y1_1[0:S-1, t in 0:k-1] <= 1)
  @variable(m, 0 <= y1_2[0:S-1, t in 1:k] <= 1)
  @variable(m, 0 <= y0_1[0:S-1, t in 0:k-1] <= 1)
  @variable(m, 0 <= y0_2[0:S-1, t in 1:k] <= 1)
  @variable(m, 0 <= z[0:k] <= 1)
  @variable(m, 0 <= l[0:S-1, 0:k] <= 1 )
  
  @constraint(m, [t in 0:k], sum(l[i,t] for i in 0:S-1) <= 1 - z[t])
  @constraint(m, [i in 0:S-1, t in 0:k], l[i,t] >= s[i,t] - z[t])

	@constraint(m, [i in 0:S-1, t in 1:k], y0_2[i,t] <= 1 - x[t])
	@constraint(m, [i in 0:S-1, t in 1:k], y1_2[i,t] <= x[t])

  @constraint(m, [i in 0:S-1, t in 1:k], y0_1[i,t-1] <= y0_2[α.τ[(i+1, 1)][1] - 1, t] )
  @constraint(m, [i in 0:S-1, t in 1:k], y1_1[i,t-1] <= y1_2[α.τ[(i+1, 2)][1] - 1, t] )

  @constraint(m, [i in 0:S-1, t in 1:k], sum(y0_1[j,t-1] for j in 0:S-1 if α.τ[(j+1, 1)][1] - 1 == i) >= y0_2[i,t] )
  @constraint(m, [i in 0:S-1, t in 1:k], sum(y1_1[j,t-1] for j in 0:S-1 if α.τ[(j+1, 2)][1] - 1 == i) >= y1_2[i,t] )

	@constraint(m, [t in 1:k], sum(y0_2[i,t] for i in 0:S-1) >= 1 - x[t] )
	@constraint(m, [t in 1:k], sum(y1_2[i,t] for i in 0:S-1) >= x[t] )

	for i in 0:S-1
		bool_i0 = true
		bool_i1 = true
		for j in 0:S-1
      α.τ[(j+1, 1)][1] - 1 == i ? bool_i0 = false : nothing
      α.τ[(j+1, 2)][1] - 1 == i ? bool_i1 = false : nothing
		end
		bool_i0 ? @constraint(m, [t in 1:k], y0_2[i,t] == 0) : nothing
		bool_i1 ? @constraint(m, [t in 1:k], y1_2[i,t] == 0) : nothing
	end

	@constraint(m, [i in 0:S-1, t in 1:k], s[i,t-1] == y0_1[i,t-1] + y1_1[i,t-1])
	@constraint(m, [i in 0:S-1, t in 1:k], s[i,t] == y0_2[i,t] + y1_2[i,t])
	

  @constraint(m, consS[i in 0:S-1], s[i,0] == 0)
	@constraint(m, consZ[t in 0:k], z[t] >= 0)
  
	@objective(m, Min, sum(z[t] for t in 0:k) )
  
  α.consS = consS
  α.consZ = consZ

end

"""Print the automata α"""
function printAutomata(α)
  println(α.n, " ", α.ni)
  for s in 1:α.n
    for i in 1:α.ni
      println(s, " " , i, " -> ", α.τ[(s, i)])
    end
  end
end 


"""Function used as an initialisation of the tree used to search for an SS of α; build the root node. In this case, the root node contains all the states"""
function ssInitNode(α)
  return Set(1:α.n)
end

"""Function used to compute the successors in the tree used to search for an SS of α of the node 'states', return the list containing, for each input symbol, the set containing all the successors of each state of 'states' """
function ssNeighbors(α, states)
  return [(i, Set([s for x in states for s in α.τ[(x, i)]])) for i in 1:α.ni]
end

"""Function used to check if a set of states is a leaf of the tree used to search for an SS of α : 'states' should contain only one state."""
function ssLeaf(states)
  return length(states) == 1
end

"""Function used to return a lower bound of the sequence merging the states of 'states' in α. seq is unused. This lower bound used the linear program associated with α. """
function ssBoundPL(α, states, seq)
  
  m = α.model
  k = 25
  kmin = ssBoundRadius(α, states, seq) # Bound using the radius

  # Use the kmin to help the linear program
  for t in 0:(kmin - 1)
    set_normalized_rhs(α.consZ[t], 1)
  end
  for t in kmin:k
    set_normalized_rhs(α.consZ[t], 0)
  end

  # Change the program so that we search for the sequence merging the states in 'states' instead of all the states. 
  for i in 0:(α.n-1)
    if (i+1) in states
      set_normalized_rhs(α.consS[i], 1)
    else
      set_normalized_rhs(α.consS[i], 0)
    end
  end

  optimize!(m)

	
  if (termination_status(m) == MOI.OPTIMAL)
    ov = objective_value(m)
    return ceil(Int64, ov - 0.0001) # Remove 0.0001 in case ov is an integer with weak precision (1.000001)
  else
    return nothing
  end
end


"""Function used to return a lower bound of the sequence merging the states of 'states' in α. seq is unused. This lower bound used the radius. """
function ssBoundRadius(α, states, seq)
  if(length(states) == 1)
    return 0
  end
  m = minimum(maximum((((s1, s2) in keys(α.μ2[z])) ? α.μ2[z][(s1, s2)] : Inf)  for s1 in states for s2 in states) for z in 1:α.n)
  return m
end

"""Function used to return a lower bound of the sequence merging the states of 'states' in α. seq is unused. This lower bound used the merging sequences. """
function ssBound(α, states, seq)
  if(length(states) == 1)
    return 0
  end
  m = maximum(α.μ[(s1, s2)] for s1 in states for s2 in states)
  return m
end

"""Function used to check if a set of states is such that the search in the tree used to search for an SS of α can be stopped: 'states' should contain only one state."""
function ssStop(states)
  return length(states) == 1
end

"""Basic search tree to search for an SS of α. initNode, neighbors, leaf, bound and stop are functions that can be used to respectively build the root, find the successors of a node in the tree, check if a node is a leaf, compute the lower bound of a node (unused) and check if the search can be stopped."""
function bfSearch(α, initNode, neighbors, leaf, bound, stop)
  nbExplore = 0
  nbCut = 0
  toExplore = [(initNode(α), [])]
  best = Inf32

  visited = Set([])

  while length(toExplore) != 0
    states, seq = popfirst!(toExplore) # Each node of the tree contains the set of states and the sequence used to go to that node.

    h = length(seq)

    if(h >= best || states ∈ visited)
      nbCut += 1
      continue
    end

    nbExplore += 1
    push!(visited, states)
    
    if(stop(states))
      if h < best
        best = h
      end
      break
    end

    if(leaf(states)) # We should never go here, as stop(states) and leaf(states) usually return the same thing
      if h < best
        best = h
      end
      continue
    end

    for (i, neighbor) in neighbors(α, states)
      push!(toExplore, (neighbor, vcat(seq, [i])))
    end

  end
  return nbExplore, nbCut, best
end


"""Search tree to search for an SS of α, improved with a bound. At each iteration, the explored node is the node with the best lower bound. initNode, neighbors, leaf, bound and stop are functions that can be used to respectively build the root, find the successors of a node in the tree, check if a node is a leaf, compute the lower bound of a node and check if the search can be stopped."""
function bestSearch(α, initNode, neighbors, leaf, bound, stop)
  nbExplore = 0
  nbCut = 0
  bounds = Dict()
  s = initNode(α)
  bounds[s] = bound(α, s, [])
  toExplore = BinaryMinHeap{Any}() # To store the nnodes of the tree sorted by their lower bound.
  iteration = 0
  push!(toExplore, (bounds[s], iteration, s, []))
  best = Inf32

  visited = Set([])

  while length(toExplore) != 0
      b, it, states, seq = pop!(toExplore) # Each node of the tree contains the associated lower bound, the number of the iteration (used so that two nodes with the same bound are not equal), the set of states and the sequence used to go to that node.

    h = length(seq)

    if(h >= best || states ∈ visited)
      nbCut += 1
      #println("X")
      continue
    end
    #println()

    nbExplore += 1
    push!(visited, states)

    if(stop(states))
      if h < best
        best = h
      end
      break
    end

    if(leaf(states)) # We should never go here, as stop(states) and leaf(states) usually return the same thing
      if h < best
        best = h
      end
      continue
    end

    for (i, neighbor) in neighbors(α, states)
      if(neighbor ∉ keys(bounds))
        bounds[neighbor] = bound(α, neighbor, seq) # We compute the bound of hte neighbor if and only if it was not already computed.
      end
      if(bounds[neighbor] + 1 < bounds[states]) # It is possible that a successor has a bound greater than the bound of states minus 1 (because of the formula of the bound does not depend on 'states'). In that case, we decrease the bound of the successor
        bounds[neighbor] = bounds[states] - 1
      end
      iteration += 1
      push!(toExplore, (h + 1+ bounds[neighbor], iteration, neighbor, vcat(seq, [i])))
    end

  end
  return nbExplore, nbCut, best
end


""" Similar to best Search but do not use a set visited to register all the visited sets. Instead use a dictionnary to store the keys of all the nodes. A visited node has a low key. We use that comparison to deduce if a node was visited or not. If the key is decreased, insteaf of dupplicating in the heap, we decrease the key : each node is explore only once."""
function bestSearchUniqueExplore(α, initNode, neighbors, bound)
  nbExplore = 0 # Nb explored nodes
  identifier = 0 # Used as a unique identifier of the keys in a heap
  
  toExplore = MutableBinaryMinHeap{Tuple{Int64,Int64,Set{Int64},Array{Any,1}}}() # Binary heap, to store the nodes of the tree sorted by the sum height + lower bound of the distance to the closest leaf. Each node of the tree is associated with that sum, a unique identifier (used so that two nodes with the same sum are not considered equal), the set of states and the sequence used to go to that node.
  storage = Dict() # Store information on nodes of the tree : height, lower bound of the distance to the closest leaf, a unique identifier (just in case), and the key in the heap

  s = initNode(α) # Root of the tree
  b = bound(α, s, []) # Bound of the tree
  key = push!(toExplore, (b, identifier, s, [])) # Insert the root in the heap
  storage[s] = (0, b, identifier, key) # Store the root in storage

  size = 1 # Always contains the number of nodes in the heap

  while size != 0
    # There are nodes that need to be explored

    hb, it, states, seq = pop!(toExplore)
    # hb : sum of the height of the node plus a lower bound to the distance to the closest leaf
    # it : a unique identifier, used by the heap in case there are two nodes with the same sum
    # states : set of states associated with that node
    # seq : sequences of inputs leading to that node from the root

    size -= 1 # There is one node less in the heap
    h = length(seq) # Height of the node

    nbExplore += 1 # One new node was explored
    
    if(length(states) == 1)
        # The node is the first leaf we encounter : we stop the algorithm
        return nbExplore, 0, h
    end
    
    if(length(states) == 2)
        # (Works when the bound is the maximum merging sequence, if states contains two states, the merging sequence is exactly the synchronizing sequence)
        return nbExplore, 0, hb
    end

    for (i, neighbor) in neighbors(α, states)
        hbn = get(storage, neighbor, nothing) # Return nothing if neighbor was never stored and the stored value otherwise
        if(isnothing(hbn))
            # We compute the bound of the neighbor if and only if it was not already computed.
            identifier += 1
            b = bound(α, neighbor, seq)
            key = push!(toExplore, (h + 1 + b, identifier, neighbor, vcat(seq, [i])))
            storage[neighbor] = (h + 1, b, identifier, key)
            size += 1
        else
            # hbn contains : the last seen height of the node, its bound, a unique identifier and the key in the heap
            hn, bn, idn, keyn = hbn
            if(h + 1 < hn)
                # If the height is better than before, we decrease the key in the heap
                update!(toExplore, keyn, (h + 1 + bn, idn, neighbor, vcat(seq, [i])))
                storage[neighbor] = (h + 1, bn, idn, keyn)
            end
        end

    end

  end
  return nbExplore, 0, Inf32
end


""" Similar to best Search but use a fibonacci heap instead of a binaryheap. Also do not use a set visited to register all the visited sets. Instead use a dictionnary to store the keys of all the nodes. A visited node has a low key. We use that comparison to deduce if a node was visited or not. If the key is decreased, insteaf of dupplicating in the heap, we decrease the key : each node is explore only once."""
function bestSearchFibo(α, initNode, neighbors, bound)
  nbExplore = 0 # Nb explored nodes
  identifier = 0 # Used as a unique identifier of the keys in a heap
  
  toExplore = MakeHeap() # Fibonacci heap, to store the nodes of the tree sorted by the sum height + lower bound of the distance to the closest leaf. Each node of the tree is associated with that sum, a unique identifier (used so that two nodes with the same sum are not considered equal), the set of states and the sequence used to go to that node.
  storage = Dict() # Store information on nodes of the tree : height, lower bound of the distance to the closest leaf, a unique identifier (just in case), and the key in the fibonacci heap

  s = initNode(α) # Root of the tree
  b = bound(α, s, []) # Bound of the tree
  key = Insert(toExplore, (b, identifier, s, [])) # Insert the root in the fibonacci heap
  storage[s] = (0, b, identifier, key) # Store the root in storage

  size = 1 # Always contains the number of nodes in the fibonacci heap

  while size != 0
    # There are nodes that need to be explored

    hb, it, states, seq = ExtractMin(toExplore).value 
    # hb : sum of the height of the node plus a lower bound to the distance to the closest leaf
    # it : a unique identifier, used by the heap in case there are two nodes with the same sum
    # states : set of states associated with that node
    # seq : sequences of inputs leading to that node from the root

    size -= 1 # There is one node less in the heap
    h = length(seq) # Height of the node

    nbExplore += 1 # One new node was explored
    
    if(length(states) == 1)
        # The node is the first leaf we encounter : we stop the algorithm
        return nbExplore, 0, h
    end
    
    if(length(states) == 2)
        # (Works when the bound is the maximum merging sequence, if states contains two states, the merging sequence is exactly the synchronizing sequence)
        return nbExplore, 0, hb
    end

    for (i, neighbor) in neighbors(α, states)
        hbn = get(storage, neighbor, nothing) # Return nothing if neighbor was never stored and the stored value otherwise
        if(isnothing(hbn))
            # We compute the bound of the neighbor if and only if it was not already computed.
            identifier += 1
            b = bound(α, neighbor, seq)
            key = Insert(toExplore, (h + 1 + b, identifier, neighbor, vcat(seq, [i])))
            storage[neighbor] = (h + 1, b, identifier, key)
            size += 1
        else
            # hbn contains : the last seen height of the node, its bound, a unique identifier and the key in the fibonacci heap
            hn, bn, idn, keyn = hbn
            if(h + 1 < hn)
                # If the height is better than before, we decrease the key in the heap
                DecreaseKey(toExplore, keyn, (h + 1 + bn, idn, neighbor, vcat(seq, [i])))
                storage[neighbor] = (h + 1, bn, idn, keyn)
            end
        end

    end

  end
  return nbExplore, 0, Inf32
end


""" Similar to best Search but use a priority queue instead of a binaryheap. Also do not use a set visited to register all the visited sets. Instead use a dictionnary to store the keys of all the nodes. A visited node has a low key. We use that comparison to deduce if a node was visited or not. If the key is decreased, insteaf of dupplicating in the heap, we decrease the key : each node is explore only once."""
function bestSearchPriorityQueue(α, initNode, neighbors, bound)
  nbExplore = 0 # Nb explored nodes
  identifier = 0 # Used as a unique identifier of the keys in a heap
  
  toExplore = PriorityQueue() # Priority Queue, to store the nodes of the tree sorted by the sum height + lower bound of the distance to the closest leaf. Each node of the tree is associated with that sum, a unique identifier (used so that two nodes with the same sum are not considered equal) and the sequence used to go to that node.
  storage = Dict() # Store information on nodes of the tree : height, lower bound of the distance to the closest leaf, a unique identifier (just in case)

  s = initNode(α) # Root of the tree
  b = bound(α, s, []) # Bound of the tree
  toExplore[s] = (b, identifier, []) # Insert the root in the priority queue
  storage[s] = (0, b, identifier) # Store the root in storage

  size = 1 # Always contains the number of nodes in the queue

  while size != 0
    # There are nodes that need to be explored

    states, value = dequeue_pair!(toExplore)
    hb, it, seq = value
    # hb : sum of the height of the node plus a lower bound to the distance to the closest leaf
    # it : a unique identifier, used by the heap in case there are two nodes with the same sum
    # states : set of states associated with that node
    # seq : sequences of inputs leading to that node from the root

    size -= 1 # There is one node less in the heap
    h = length(seq) # Height of the node

    nbExplore += 1 # One new node was explored
    
    if(length(states) == 1)
        # The node is the first leaf we encounter : we stop the algorithm
        return nbExplore, 0, h
    end
    
    if(length(states) == 2)
        # (Works when the bound is the maximum merging sequence, if states contains two states, the merging sequence is exactly the synchronizing sequence)
        return nbExplore, 0, hb
    end

    for (i, neighbor) in neighbors(α, states)
        hbn = get(storage, neighbor, nothing) # Return nothing if neighbor was never stored and the stored value otherwise
        if(isnothing(hbn))
            # We compute the bound of the neighbor if and only if it was not already computed.
            identifier += 1
            b = bound(α, neighbor, seq)
            toExplore[neighbor] = (h + 1 + b, identifier, vcat(seq, [i]))
            storage[neighbor] = (h + 1, b, identifier)
            size += 1
        else
            # hbn contains : the last seen height of the node, its bound and a unique identifier
            hn, bn, idn = hbn
            if(h + 1 < hn)
                # If the height is better than before, we decrease the key in the queue
                toExplore[neighbor] = (h + 1 + bn, idn, vcat(seq, [i]))
                storage[neighbor] = (h + 1, bn, idn)
            end
        end

    end

  end
  return nbExplore, 0, Inf32
end

"""Search tree to search for PriorityQueue{Any, Any}(), using a depht first search and a branch and bound algorithm. At each iteration, we explore the left-most unexplored node. initNode, neighbors, leaf, and bound are functions that can be used to respectively build the root, find the successors of a node in the tree, check if a node is a leaf, and compute the lower bound of a node."""
function dfSearch(α, initNode, neighbors, leaf, bound)

  nbExplore = 0
  nbCut = 0
  best = Inf32
  bounds = Dict()
  heights = Dict() # Lower height in which the state appeared in the tree
  values = Dict() # Best size of the sequence found to merge each set of states

  function explore(states, seq, preds)
    h = length(seq)

    if(states in keys(values) && bounds[states] == values[states])
      v = values[states]
      if(h + v < best)
        best = h + v
      end
      nbCut += 1
      return values[states]
    end

    if(states in keys(heights))
      ph = heights[states]
      if(ph < h) # L'élément a été aperçu plus haut.
        if(states in keys(values))
          return values[states]
        else
          return Inf32
        end
      else
        heights[states] = h # On met à jour la hauteur
      end
    else
      heights[states] = h # On note que cet état est apparu à la hauteur h
    end

    b = bounds[states]
    if(best <= h + b) # On vérifie la borne
      nbCut += 1
      return nothing # On renvoie nothing plutot que l'infini, car on ne peut pas déduire la taille de la séquence descedant de l'état. 
    end

    nbExplore += 1 # Si la borne est assez élevée, on explore
    
    if(leaf(states))
      if(h < best)
        # On met à jour la meilleure solution trouvée 
        best = h
      end
      values[states] = 0
      return 0
    end

    β = nothing # Taille de la séquence descendant de cet état

    mb = Inf32

    preds2 = vcat(preds, [states])
    for (i, neighbor) in neighbors(α, states)
      seq2 = vcat(seq, [i])
      if(neighbor ∉ keys(bounds))
        bounds[neighbor] = bound(α, neighbor, seq) # Calcul des bornes des descendants si nécessaire
      end
      
      v = explore(neighbor, seq2, preds2) # Calcul de la taille de la séquence issue du descendant
      mb = min(mb, bounds[neighbor])
      if(!isnothing(v) && (isnothing(β) || v + 1 < β)) 
        β = v + 1
      end
    end

    if(mb + 1 > b)
      for i in 1:length(preds2)
        bounds[length(preds2) - i + 1] = mb + i
      end
    end

    if(!isnothing(β))
      values[states] = β
    end
    return β

  end

  states = initNode(α)
  bounds[states] = 0
  explore(states, [], [])

  return nbExplore, nbCut, best
end

"""Test all the algorithms on many automatas"""
function tests()

  NB_TESTS = 30 
  ns = [70]  # Size of the tested automatas
  nis = [2]#, 3, 4]  # Size of the tested alphabets. For each value in nis and each n in ns a bunch of NB_TESTS automatas are build

  # To store which algorithm explore less number of nodes that which algorithm 
  exploreWin = [0 0 0 0 0; 
                0 0 0 0 0 ; 
                0 0 0 0 0 ; 
                0 0 0 0 0 ; 
                0 0 0 0 0]
  

  # To store, when an algorithm explores less than another algorithm, how many nodes less did it explore.
  sume2me1 = [
              [[], [], [], [], []], 
              [[], [], [], [], []], 
              [[], [], [], [], []], 
              [[], [], [], [], []], 
              [[], [], [], [], []]]

  # Number of time the algorithms did not return the same value
  errors = 0

  index = 0

  # Number of automatas ignored because they do not have an SS
  nb_ignored_α = 0
  for n in ns
    for ni in nis
      for number in 1:50
        index += 1
        println(index, "/50");
        α = generateDeterministic(n, ni)
        if(!hasMergingSequence(α))
          nb_ignored_α += 1
          continue
        end
        
        println("Explore");
        # Init all the parameters in case some algorithms are not run.
        nbExplore1, nbCut1, v1 = 0, 0, 0
        nbExplore2, nbCut2, v2 = 0, 0, 0
        nbExplore3, nbCut3, v3 = 0, 0, 0
        nbExplore4, nbCut4, v4 = 0, 0, 0
        nbExplore5, nbCut5, v5 = 0, 0, 0
        nbExplore6, nbCut6, v6 = 0, 0, 0
        nbExplore7, nbCut7, v7 = 0, 0, 0
        nbExplore8, nbCut8, v8 = 0, 0, 0

        # Do the measures
        
        #time = @elapsed((nbExplore1, nbCut1, v1) = bfSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBoundRadius, ssStop))
        #println(time, ' ',nbExplore1, ' ',v1)
        #nbExplore2, nbCut2, v2 = dfSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBoundRadius)
        #nbExplore3, nbCut3, v3 = dfSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBoundPL)
        #time4 = @elapsed((nbExplore4, nbCut4, v4) = bestSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBound, ssStop))
        #println(time4, ' ',nbExplore4, ' ', v4)
        #nbExplore5, nbCut5, v5 = bestSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBoundPL, ssStop)
        #time6 = @elapsed((nbExplore6, nbCut6, v6) = bestSearchFibo(α, ssInitNode, ssNeighbors, ssBound))
        #println(time6, ' ',nbExplore6, ' ', v6)
        #time7 = @elapsed((nbExplore7, nbCut7, v7) = bestSearchPriorityQueue(α, ssInitNode, ssNeighbors, ssBound))
        #println(time7, ' ',nbExplore7, ' ', v7)
        time8 = @elapsed((nbExplore8, nbCut8, v8) = bestSearchUniqueExplore(α, ssInitNode, ssNeighbors, ssBound))
        println(time8, ' ',nbExplore8, ' ', v8)
        time9 = @elapsed((nbExplore9, nbCut9, v9) = bestSearchUniqueExplore(α, ssInitNode, ssNeighbors, ssBoundRadius))
        println(time9, ' ',nbExplore9, ' ', v9)
        
        # Check if there is an error
        if(length(Set([v8, v9])) != 1)
          errors += 1
          
          printAutomata(α)

          println(α.τm1)
          println(α.μ2)
          for i in 1:α.n
              for j in (i+1):α.n
                println(ssBoundRadius(α, i:j, []))
            end
        end

          println("BFS   : ", v1)
          println("DFMS   : ", v2)
          println("DFPL : ", v3)
          println("BestMS : ", v4)
          println("BestPL : ", v5)
          return
          
        else
          # Count which algorithms explored less
          explores = [nbExplore1, nbExplore2, nbExplore3, nbExplore4, nbExplore5]
          for i in 1:5
            for j in 1:5
              if i == j
                continue
              end
              if explores[i] < explores[j]
                exploreWin[i, j] += 1
                push!(sume2me1[i][j], (explores[j] - explores[i]) * 100.0 / explores[j])
              end
            end
          end

        end
        #=
        # Print the results

        println(index, " / ", 150, " (", n, ", ", ni, ", ", number, ")")
        println("NB Automatas without SS : ", nb_ignored_α)
        println("Explorations :")
        display(exploreWin)
        println("Mean exploration Word - Best / Worst * 100 : ")
        diffmeans = [0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
        for i in 1:5
          for j in 1:5
            if i == j
              continue
            end
            diff = sume2me1[i][j]
            ldiff = length(diff)
            diffmeans[i, j] = ldiff != 0 ? sum(diff) / length(diff) : 0
          end
        end
        display(diffmeans)
      =#
      end
      println("Errors : ", errors)
      println()
    end
      display(exploreWin)
  end
end

#tests()


#=

α = Automata()
α.n = 5
α.ni = 2
push!(α.τ[(1, 1)], 4)
push!(α.τ[(1, 2)], 2)
push!(α.τ[(2, 1)], 3)
push!(α.τ[(2, 2)], 2)
push!(α.τ[(3, 1)], 4)
push!(α.τ[(3, 2)], 5)
push!(α.τ[(4, 1)], 2)
push!(α.τ[(4, 2)], 5)
push!(α.τ[(5, 1)], 4)
push!(α.τ[(5, 2)], 3)

push!(α.τm1[(4, 1)], 1)
push!(α.τm1[(2, 2)], 1)
push!(α.τm1[(3, 1)], 2)
push!(α.τm1[(2, 2)], 2)
push!(α.τm1[(4, 1)], 3)
push!(α.τm1[(5, 2)], 3)
push!(α.τm1[(2, 1)], 4)
push!(α.τm1[(5, 2)], 4)
push!(α.τm1[(4, 1)], 5)
push!(α.τm1[(3, 2)], 5)
buildMergingSequences(α)
buildbackwardMergingSequences(α)
buildPL(α)

nbExplore1, nbCut1, v1 = bfSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBound, ssStop)
nbExplore2, nbCut2, v2 = dfSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBound)
nbExplore3, nbCut3, v3 = dfSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBoundPL)
nbExplore4, nbCut4, v4 = dfSearch(α, ssInitNode, ssNeighbors, ssLeaf, ssBoundRadius)
        
println(nbExplore1, " ", nbCut1, " ",  v1)
println(nbExplore2, " ", nbCut2, " ", v2)
println(nbExplore3, " ", nbCut3, " ", v3)
println(nbExplore4, " ", nbCut4, " ", v4)
=#


function test2()

  fsm = read_fsm("../data/cerny_fsm/cerny_n20.fsm")
  #println(fsm)
  
  a = generateAutomataFSM(fsm)
  #printAutomata(a)

  debms = Dates.now()
  nbExplore, nbCut, res = bestSearch(a, ssInitNode, ssNeighbors, ssLeaf, ssBound, ssStop)
  finms = Dates.now()
  println("explore: ", nbExplore,  "time: ", Dates.value(finms-debms), "ms")


#=
  dir_fsm = "../data/SS_50fsm_n100/"
  list_fsm = readdir(dir_fsm)

  count_fsm = 0
  ttl_PL = 0
  ttl_rayon = 0
  ttl_ms = 0

  for i in list_fsm
    if occursin( ".fsm", i )
      #println(dir_fsm*i)
      fsm = read_fsm(dir_fsm*i)
      a = generateAutomataFSM(fsm)
      count_fsm += 1

      if(!hasMergingSequence(a))
        continue
      end
      
      #=
      debPL = Dates.now()
      nbExplore, nbCut, res = bestSearch(a, ssInitNode, ssNeighbors, ssLeaf, ssBoundPL, ssStop)
      finPL = Dates.now()
      println("fsm: ", i, " time: ", Dates.value(finPL-debPL), "ms")

      debrayon = Dates.now()
      nbExplore, nbCut, res = bestSearch(a, ssInitNode, ssNeighbors, ssLeaf, ssBoundRadius, ssStop)
      finrayon = Dates.now()
      println("fsm: ", i, " time: ", Dates.value(finrayon-debrayon), "ms")
      =#

      debms = Dates.now()
      nbExplore, nbCut, res = bestSearch(a, ssInitNode, ssNeighbors, ssLeaf, ssBound, ssStop)
      finms = Dates.now()
      println("fsm: ", i, " time: ", Dates.value(finms-debms), "ms")

      if (count_fsm != 1)
        #ttl_PL += Dates.value(finPL-debPL)
        #ttl_rayon += Dates.value(finrayon-debrayon)
        ttl_ms += Dates.value(finms-debms)
      end
    end
  end
  println("total time ms: ", ttl_ms, "ms  average: ", ttl_ms/count_fsm, "ms")
=#

  #println("total time PL: ", ttl_PL, "ms  average: ", ttl_PL/count_fsm, "ms")
  #println("total time rayon: ", ttl_rayon, "ms  average: ", ttl_rayon/count_fsm, "ms")

end

tests()

