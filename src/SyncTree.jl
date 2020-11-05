include("util.jl")


# return min length SS (DFA, PFA and NFA)
function parcours_largeur(fsm)
	
	S = fsm.s

	#key: set of states s that we need to check / value: shortest seq to get to s
	to_visit = [(Set([i for i in 0:S-1]), "")]

	#key: set of states s already visited / value: shortest seq to get to s
	visited = Dict(Set([i for i in 0:S-1]) => "")

	while length(to_visit) != 0

		# current node / subset
		current = popfirst!(to_visit)

		# c0 and c1: successors of current
		c0 = apply_input(fsm, current[1], 0)
		c1 = apply_input(fsm, current[1], 1)
		
		# seq to get to c0 and c1
		seq_c0 = current[2]*"0"
		seq_c1 = current[2]*"1"
		# if S in a subset: undefined transition (pfa)

		if !(S in c0)
			# c0 is a SS 
			length(c0) == 1 ? (return seq_c0) : nothing

			# check if set c0 already visited
			if !(c0 in keys(visited))
				push!( to_visit, (c0,seq_c0) )
				visited[c0] = seq_c0
			end
		end

		if !(S in c1)
			# c1 is a SS
			length(c1) == 1 ? (return seq_c1) : nothing

			# check if set c1 already visited
			if !(c1 in keys(visited))
				push!( to_visit, (c1,seq_c1) )
				visited[c1] = seq_c1
			end
		end

	end

	# no SS
	return ""

end

# return min size of S' and the seq to get S' (DFA, PFA and NFA)
# same as parcours_largeur, but more infors about fsm without SS
function min_Sprime(fsm)
	
	S = fsm.s

	#key: set of states s that we need to check / value: shortest seq to get to s
	to_visit = [(Set([i for i in 0:S-1]), "")]

	#key: set of states s already visited / value: shortest seq to get to s
	visited = Dict(Set([i for i in 0:S-1]) => "")

	current_res = S
	current_seq = ""

	while length(to_visit) != 0

		# current node / subset
		current = popfirst!(to_visit)

		# check if current is better
		if length(current[1]) < current_res
			current_res = length(current[1])
			current_seq = current[2]
		end

		# c0 and c1: successors of current
		c0 = apply_input(fsm, current[1], 0)
		c1 = apply_input(fsm, current[1], 1)
		
		# seq to get to c0 and c1
		seq_c0 = current[2]*"0"
		seq_c1 = current[2]*"1"
		# if S in a subset: undefined transition (pfa)

		if !(S in c0)
			# c0 is a SS 
			length(c0) == 1 ? (return 1, seq_c0) : nothing

			# check if set c0 already visited
			if !(c0 in keys(visited))
				push!( to_visit, (c0,seq_c0) )
				visited[c0] = seq_c0
			end
		end

		if !(S in c1)
			# c1 is a SS
			length(c1) == 1 ? (return 1, seq_c1) : nothing

			# check if set c1 already visited
			if !(c1 in keys(visited))
				push!( to_visit, (c1,seq_c1) )
				visited[c1] = seq_c1
			end
		end

	end

	return current_res, current_seq

end

# return min size of S' for all k and the seq to get S'
# same as min_Sprime, but additional infos about the min length for each k (sequence size)
function min_Sprime_allk(fsm)
	
	S = fsm.s

	res_k = [S]

	#key: set of states s that we need to check / value: shortest seq to get to s
	to_visit = [(Set([i for i in 0:S-1]), "")]

	#key: set of states s already visited / value: shortest seq to get to s
	visited = Dict(Set([i for i in 0:S-1]) => "")

	current_res = S
	current_seq = ""

	while length(to_visit) != 0

		# current node / subset
		current = popfirst!(to_visit)

		# check if current is better
		if length(current[1]) < current_res
			current_res = length(current[1])
			current_seq = current[2]
		end

		# c0 and c1: successors of current
		c0 = apply_input(fsm, current[1], 0)
		c1 = apply_input(fsm, current[1], 1)
		
		# seq to get to c0 and c1
		seq_c0 = current[2]*"0"
		seq_c1 = current[2]*"1"
		# if S in a subset: undefined transition (pfa)


		# check if c0 or c1 min length for current k
		if length(res_k) < length(seq_c0)
			append!(res_k, min(last(res_k),length(c0), length(c1)))
		else
			res_k[length(seq_c0)] > length(c0) ? res_k[length(seq_c0)] = length(c0) : nothing
			res_k[length(seq_c1)] > length(c1) ? res_k[length(seq_c1)] = length(c1) : nothing
		end


		if !(S in c0)
			# c0 is a SS 
			length(c0) == 1 ? (return 1, seq_c0, res_k) : nothing

			# check if set c0 already visited
			if !(c0 in keys(visited))
				push!( to_visit, (c0,seq_c0) )
				visited[c0] = seq_c0
			end
		end

		if !(S in c1)
			# c1 is a SS
			length(c1) == 1 ? (return 1, seq_c1, res_k) : nothing

			# check if set c1 already visited
			if !(c1 in keys(visited))
				push!( to_visit, (c1,seq_c1) )
				visited[c1] = seq_c1
			end
		end

	end

	return current_res, current_seq, res_k

end
############ change res_k => get length + seq


# return shortest S' synchronizing word
function parcours_largeur_Sprime(fsm, set_Sprime)
	
	S = fsm.s

	if set_Sprime == Set([i for i in 0:S-1])
		return ""
	end

	#key: set of states s that we need to check / value: shortest seq to get to s
	to_visit = [(Set([i for i in 0:S-1]), "")]

	#key: set of states s already visited / value: shortest seq to get to s
	visited = Dict(Set([i for i in 0:S-1]) => "")

	while length(to_visit) != 0

		# current node / subset
		current = popfirst!(to_visit)

		# c0 and c1: successors of current
		c0 = apply_input(fsm, current[1], 0)
		c1 = apply_input(fsm, current[1], 1)
		
		# seq to get to c0 and c1
		seq_c0 = current[2]*"0"
		seq_c1 = current[2]*"1"
		# if S in a subset: undefined transition (pfa)

		if !(S in c0)
			# c0 is a S' Synch 
			isempty(setdiff(c0,set_Sprime)) ? (return seq_c0) : nothing

			# check if set c0 already visited
			if !(c0 in keys(visited))
				push!( to_visit, (c0,seq_c0) )
				visited[c0] = seq_c0
			end
		end

		if !(S in c1)
			# c1 is a S' Synch
			isempty(setdiff(c1,set_Sprime)) ? (return seq_c1) : nothing

			# check if set c1 already visited
			if !(c1 in keys(visited))
				push!( to_visit, (c1,seq_c1) )
				visited[c1] = seq_c1
			end
		end

	end

	# no SS
	return ""

end


######################################## TEST ########################################
#=

fsm = read_fsm("../data/fsm_hss.fsm")
#fsm = read_fsm("../data/cerny_n3.fsm")
#fsm = read_fsm("../data/pfa.fsm")
#fsm = read_fsm("../data/nfa.fsm")

#apply_input(fsm, Set([0,1,2,3]), 0)
#apply_seq(fsm, Set([0,1,2,3]), "01010")

#parcours_largeur(fsm)
#min_Sprime(fsm)
#min_Sprime_allk(fsm)
parcours_largeur_Sprime(fsm, Set([3]))
#parcours_largeur_Sprime(fsm, Set([0,2]))

=#