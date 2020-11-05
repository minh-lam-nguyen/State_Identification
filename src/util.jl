######################################## FSM ########################################

struct FSM
	type # 0: dfa (complet + determinist) / 1: nfa / 2: pfa
	s # number of states
	i # number of input
	o # number of output
	p # number of transition
	succ # successors state(s) 
	outp # successors output
end

#= struc FSM infos

	!!!! warning about succ and outp : !!!!
	to get successor of s with input i : succ[s+1][i+1]
	exemple: succ[1][2] = 4 means successor of s=0 after applying input i=1 is 4

	dfa: 
	succ[s+1][i+1] -> Int

	nfa:
	succ[s+1][i+1] -> Dict( Int => Set(state) )

=#


# read .fsm file and return FSM
function read_fsm(filename)

	#load FSM from .fsm
	open(filename, "r") do file

		#readline(file) # F
		pass_comment = true

		# 0 : complet/determinist , 1 complet/non determinist
		fsm_type = -1
		while pass_comment
			temp_line = readline(file)
			if temp_line[1] == 'F'
				pass_comment = false
				fsm_type = temp_line[3]
			end
		end

		s = parse(Int64, split(readline(file), ' ')[2])
		i = parse(Int64, split(readline(file), ' ')[2])
		o = parse(Int64, split(readline(file), ' ')[2])
		readline(file) # n0
		p = parse(Int64, split(readline(file), ' ')[2])


		# DFA
		if fsm_type == '0'
			succ = [[0 for inp in 1:i] for u in 1:s]
			outp = [[0 for inp in 1:i] for u in 1:s]
			for j in 1:p
				l = split(readline(file),' ')
				u, inp, v, out = parse(Int64, l[1]), parse(Int64, l[2]), parse(Int64, l[3]), parse(Int64, l[4])
				succ[u+1][inp+1] = v
				outp[u+1][inp+1] = out
			end

			return FSM(fsm_type,s,i,o,p,succ,outp)
		end

		# NFA or PFA
		if fsm_type == '1' || fsm_type == '2'
			succ = Dict( s => Dict( inp => Set() for inp in 1:i) for s in 1:s )
			outp = Dict( s => Dict( inp => Set() for inp in 1:i) for s in 1:s )
			for j in 1:p
				l = split(readline(file),' ')
				u, inp, v, out = parse(Int64, l[1]), parse(Int64, l[2]), parse(Int64, l[3]), parse(Int64, l[4])
				push!(succ[u+1][inp+1], v)
				push!(outp[u+1][inp+1], out)
			end

			return FSM(fsm_type,s,i,o,p,succ,outp)
		end

	end

end


# write .fsm file for the FSM
function write_fsm_file(filename, fsm)

	open(filename, "w") do file
		println(file, "F ", fsm.type)
		println(file, "s ", fsm.s)
		println(file, "i ", fsm.i)
		println(file, "o ", fsm.o)
		println(file, "n0 0")
		println(file, "p ", fsm.p)			

		for x in 0:fsm.s-1
			for x0 in fsm.succ[x+1][1]
				println(file, x, " 0 ", x0, " 0" )
			end
			for x1 in fsm.succ[x+1][2]
				println(file, x, " 1 ", x1, " 0" )
			end
		end
	end

end


# generate random nfa, nb_nondet_states: nb of non determinist states
# warning: no output (SS)
function random_nfa(nb_states, n_inp, nb_nondet_states)

	succ = Dict( s => Dict(inp => Set() for inp in 1:n_inp) for s in 1:nb_states )
	#outp = Dict( s => Dict(inp => Set() for inp in 1:n_inp) for s in 1:nb_states )

	nondet_states = Set()
	for _ in 0:nb_nondet_states-1
		push!(nondet_states, rand(0:nb_states-1))
	end

	for x in 0:nb_states-1
		if x in nondet_states
			nb_transitions = rand(1:nb_states-1)
			#println("x ", x)
			#println("nb_tr ", nb_transitions)

			for _ in 0:nb_transitions-1
				push!(succ[x+1][rand(1:2)], rand(0:nb_states-1))
				#push!(outp[x+1][rand(1:2)], rand(0:1))
			end
		end

		push!(succ[x+1][1], rand(0:nb_states-1))
		push!(succ[x+1][2], rand(0:nb_states-1))
		#push!(outp[x+1][1], rand(0:1))
		#push!(outp[x+1][2], rand(0:1))
	end	

	count_p = 0
	for i in values(succ)
		for j in values(i)
			count_p += length(j)
		end
	end

	return FSM( "1", nb_states, n_inp, 2, count_p, succ, Dict() )

end


# generate cerny automaton with n states
function cerny(n)

	succ = [[0, 0] for u in 1:n]
	outp = [[0, 0] for u in 1:n]

	succ[1][1] = 1
	succ[1][2] = 1

	for i in 2:n-1
		succ[i][1] = i-1
		succ[i][2] = i
	end

	succ[n][1] = n-1
	succ[n][2] = 0

	return FSM(0, n, 2, 2, 2*n, succ, outp)

end


# return set of state S', S_set x inp -> S'
function apply_input(fsm, S_set, inp)

	res = Set()
	for i in S_set
		# complet / determinist
		if fsm.type == '0'
			push!(res, fsm.succ[i+1][inp+1])
		end

		# nfa or pfa
		if fsm.type == '1' || fsm.type == '2'
			isempty(fsm.succ[i+1][inp+1]) ? push!(res, fsm.s) : nothing
			for succ_i in fsm.succ[i+1][inp+1]
				push!(res, succ_i)
			end
		end

	end
	return res

end


# return set of states S', S_set x seq -> S'
function apply_seq(fsm, S_set, seq)

	res = Set(S_set)
	for i in seq
		#println(res)
		res = apply_input(fsm, res, parse(Int64,i))
	end
	return res
	
end

######################################## TEST ########################################
#=

	#fsm = read_fsm("../data/fsm_hss.fsm")

	#fsm = read_fsm("../data/nfa.fsm")
	#fsm.succ[1]

	#cerny(3)

	#random_nfa(4,2,1)

=#