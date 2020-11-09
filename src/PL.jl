using JuMP, Gurobi
using LightGraphs

include("util.jl")


### !!! graph G: nodes 1:S and not 0:S-1
# return graph G of fsm
# fsm : FSM
# return SimpleDiGraph
function get_graph(fsm)
	
	S = fsm.s

	# create graph with S nodes
	G = SimpleDiGraph(S)
	#G = MetaDiGraph(g,S)

	# add edges
	for i in 1:S
		add_edge!(G, i, fsm.succ[i][1]+1)
		add_edge!(G, i, fsm.succ[i][2]+1)

	end

	return G

end


# return l the longest merging sequence that end in the set final (all states if final not defined) of fsm, 
# 		-1 if cant find merging seq between 2 states
# fsm : FSM, final : Set{Int}
# return : Int
function longest_merging(fsm, final=[i for i in 1:fsm.s])
	
	S = fsm.s
	G = SimpleDiGraph( Int64( S*(S+1) / 2 ) )
	
	# Dict: keys: pairs "uv", values: nodes in G
	pairs = Dict()
	count = 1
	for i in 1:S
		for j in i:S
			pairs["$i$j"] = count
			count += 1
		end
	end
	#println(count)
	#println(pairs)

	# set edges
	for i in 1:S
		for j in i+1:S

			i0 = fsm.succ[i][1]+1
			j0 = fsm.succ[j][1]+1
			if i0 <= j0
				add_edge!( G, pairs[string(i0)*string(j0)], pairs[string(i)*string(j)])
			else
				add_edge!( G, pairs[string(j0)*string(i0)], pairs[string(i)*string(j)])
			end

			i1 = fsm.succ[i][2]+1
			j1 = fsm.succ[j][2]+1
			if i1 <= j1
				add_edge!( G, pairs[string(i1)*string(j1)], pairs[string(i)*string(j)])
			else
				add_edge!( G, pairs[string(j1)*string(i1)], pairs[string(i)*string(j)])
			end

		end
	end


	# get longest merging seq l
	l = 0

	T = [pairs[string(i)*string(i)] for i in final]
	rest = [pairs[i] for i in keys(pairs) if !(pairs[i] in T)]

	cond = true # false if there is (at least) one element in rest that is not accessible
	while length(rest) > 0 && cond 

		temp = []
		cond = false
		for e in edges(G)
			if src(e) in T
				append!(temp, dst(e))
				ind = findfirst(x->x==dst(e), rest)
				if ind != nothing
					cond = true
					deleteat!( rest, ind )
				end
			end
		end

		T = temp
		l += 1

	end

	return l 

end
############ change dict (keys)


# return l the min rayon from any states to a state in final (all states if final not defined)
# 		-1 if cant find any path between u-v, v in final
# fsm : FSM, final : Set{Int}
# return : Int
function min_rayon(fsm, G, final=[i for i in 1:fsm.s])
	
	# min rayon
	res = (fsm.s)^2

	#bfs search
	for i in final

		temp_length = -1
		checked = []
		toCheck = [i]
		
		while length(toCheck) != 0
			for _ in 1:length(toCheck)
				temp = popfirst!(toCheck)
				append!(checked,temp)
				for j in inneighbors(G,temp)
					if !(j in checked) && !(j in toCheck)
						append!(toCheck,j)
					end
				end

			end
			temp_length += 1
		end

		println(temp_length, " ", i)
		if temp_length < res && temp_length != 0
			res = temp_length
		end

	end

	return res

end


# PL : return the min length SS (<= k), kmin : borne min, final : final states
# fsm : FSM, k : Int, kmin : Int, final : Set{Int}, relax : Bool
# return : Int
function compact(fsm, k, kmin=0, final=[i for i in 1:fsm.s], relax=false)

	if kmin > k || kmin < 0
		return -1, kmin
	end

	m = Model(Gurobi.Optimizer)

	# 0 : noCut
	set_optimizer_attribute(m, "Cuts", 0)


	########## VAR ##########

	S = fsm.s

	if relax

		## x[t] with 1 <= t <= k
		@variable(m, 0 <= x[1:k] <= 1)
		#@variable(m, x[1:k], Bin)

		## s[i,t] i in S, 0 <= t <= k
		@variable(m, 0 <= s[0:S-1, 0:k] <= 1)

		## z[t] 1 <= t <= k
		@variable(m, 0 <= z[0:k] <= 1)

		# cons_z
		#@constraint(m, cons_sum[i in 0:S-1, j in i+1:S-1, t in kmin:k], s[i,t] + s[j,t] <= 1 + z[t] )
		#lambda
		@variable(m, 0 <= l[0:S-1, kmin:k] <=1 )
		@constraint(m, cons_l1[t in kmin:k], sum(l[i,t] for i in 0:S-1) <= 1 - z[t])
		@constraint(m, cons_l2[i in 0:S-1, t in kmin:k], l[i,t] >= s[i,t] - z[t])

	else

		## x[t] with 1 <= t <= k
		@variable(m, x[1:k], Bin)

		## s[i,t] i in S, 0 <= t <= k
		@variable(m, s[0:S-1, 0:k], Bin)

		## z[t] 1 <= t <= k
		@variable(m, z[0:k], Bin)

		# cons_z
		@constraint(m, cons_sum[t in kmin:k], sum(s[i,t] for i in 0:S-1) <= 1 + z[t] * (S-1) )

	end


	########## CONS ##########

	# composante fortement connexe: on converge vers un etat de C
	@constraint(m, sum(s[i-1,k] for i in final ) == 1 )


	## s[i,t] = 1 iff i appears at time t

	@constraint(m, cons_s0[i in 0:S-1], s[i,0] == 1)

	@constraint( m, cons_j1[j in 0:S-1, t in 1:k], s[fsm.succ[j+1][2],t] >= s[j,t-1] + x[t] - 1 )

	@constraint( m, cons_j0[j in 0:S-1, t in 1:k], s[fsm.succ[j+1][1],t] >= s[j,t-1] - x[t] )

	## z[t] = 1 for all 0 <= t < kmin
	@constraint(m, [t in 0:kmin-1], z[t] == 1)

	## z_k = 0
	@constraint(m, z[k]==0)


	########## INEG ##########

	## (1) 

	#@constraint(m, [j in 0:S-1, t in 1:k], s[fsm.succ[j+1][1],t] + s[fsm.succ[j+1][2],t] >= s[j,t-1])
	#@constraint(m, [j in 0:S-1, t in 1:k], s[fsm.succ[j+1][1],t] >= Int64(fsm.succ[j+1][1] == fsm.succ[j+1][2]) * s[j,t-1])


	## (2)

	#@constraint(m, [t in 1:k], sum( s[fsm.succ[j+1][1],t] for j in 0:S-1 ) >= 1 - x[t] )
	#@constraint(m, [t in 1:k], sum( s[fsm.succ[j+1][2],t] for j in 0:S-1 ) >= x[t] )
	

	##flot 

	##(A1)+:
	#@variable(m, 0 <= f1[0:S-1, 0:k] <= 1)
	#@variable(m, 0 <= g1[0:S-1, 0:k] <= 1)
	#@constraint(m, cons_f1[i in 0:S-1, t in 0:k], f1[i,t] <= s[i,t])
	#@constraint(m, cons_g1[i in 0:S-1, t in 0:k], g1[i,t] <= s[i,t])
	#@constraint(m, cons_a1[i in 0:S-1, t in 1:k], g1[i,t] == sum(f1[j, t-1] for j in 0:S-1 if fsm.succ[j+1][2] == i) )
	#@constraint(m, cons_a1_2[t in 1:k], sum(f1[i,t-1] for i in 0:S-1) == x[t])
	##(A0)+:
	#@variable(m, 0 <= f0[0:S-1, 0:k] <= 1)
	#@variable(m, 0 <= g0[0:S-1, 0:k] <= 1)
	#@constraint(m, cons_f0[i in 0:S-1, t in 0:k], f0[i,t] <= s[i,t])
	#@constraint(m, cons_g0[i in 0:S-1, t in 0:k], g0[i,t] <= s[i,t])
	#@constraint(m, cons_a0[i in 0:S-1, t in 1:k], g0[i,t] == sum(f0[j, t-1] for j in 0:S-1 if fsm.succ[j+1][1] == i) )
	#@constraint(m, cons_a0_2[t in 1:k], sum(f0[i,t-1] for i in 0:S-1) == 1 - x[t])
	#
	#@constraint(m, [i in 0:S-1, t in 1:k], f0[i,t-1]+f1[i,t-1] <= s[i,t-1] )
	#@constraint(m, [i in 0:S-1, t in 1:k], g0[i,t]+g1[i,t] <= s[i,t] )


	## (3)

	#@constraint(m, [i in 0:S-1, j in 0:S-1, l in 0:S-1, t in 1:k], s[i,t] >= s[j,t-1] * Int64(fsm.succ[j+1][2]==i) + s[l,t-1] * Int64(fsm.succ[l+1][1]==i) - 1 )

	## (4)

	#@constraint(m, [j in 0:S-1, t in 1:k], sum( s[i,t-1] for i in 0:S-1 if fsm.succ[i+1][2]==j ) >= s[j,t] + x[t] - 1 )
	#@constraint(m, [j in 0:S-1, t in 1:k], sum( s[i,t-1] for i in 0:S-1 if fsm.succ[i+1][1]==j ) >= s[j,t] - x[t] )


	########## OBJ ########## 

	@objective(m, Min, sum(z[t] for t in 0:k) )


	########## optimize ##########

	optimize!(m)
	
	#println(value.(s))
	#println(value.(z))
	#println(value.(x))


	if (termination_status(m) == MOI.OPTIMAL)
		return objective_value(m)
	else
		return -1
	end

end


# PL : return the min size of S' we can get after applying a seq of length k
# fsm : FSM, k : Int, relax : Bool
# return : Int
function min_Sprime(fsm, k, relax=false)

	m = Model(Gurobi.Optimizer)

	# 0 : noCut
	set_optimizer_attribute(m, "Cuts", 0)


	#### VAR ####

	S = fsm.s

	if relax

		## x[t] with 1 <= t <= k
		@variable(m, 0 <= x[1:k] <= 1)
		#@variable(m, x[1:k], Bin)

		## s[i,t] i in S, 0 <= t <= k
		@variable(m, 0 <= s[0:S-1, 0:k] <= 1)

	else

		## x[t] with 1 <= t <= k
		@variable(m, x[1:k], Bin)

		## s[i,t] i in S, 0 <= t <= k
		@variable(m, s[0:S-1, 0:k], Bin)

	end


	#### CONS ####

	## s[i,t] = 1 iff i appears at time t

	@constraint( m, cons_s0[i in 0:S-1], s[i,0] == 1 )

	@constraint( m, cons_j1[j in 0:S-1, t in 1:k], s[fsm.succ[j+1][2],t] >= s[j,t-1] + x[t] - 1 )

	@constraint( m, cons_j0[j in 0:S-1, t in 1:k], s[fsm.succ[j+1][1],t] >= s[j,t-1] - x[t] )


	#### OBJ #### 

	@objective(m, Min, sum(s[i,k] for i in 0:S-1) )


	#### optimize ####

	optimize!(m)
	
	#println(value.(s))
	#println(value.(x))


	if (termination_status(m) == MOI.OPTIMAL)
		return objective_value(m)
	else
		return -1
	end

end


# PL : return the min length SS (<= k), kmin : borne min, final : final states
# fsm : FSM, k : Int, kmin : Int, final : Set{Int}, relax : Bool
# return : Int
function env_conv(fsm, k, kmin=0, final=[i for i in 1:fsm.s], relax=false)

	if kmin > k || kmin < 0
		return -1, kmin
	end

	S = fsm.s


	####################################################### PL #######################################################

	m = Model(Gurobi.Optimizer)

	# no Cut
	set_optimizer_attribute(m, "Cuts", 0)

	####################################################### VAR #######################################################

	if relax

		## x[t] with 1 <= t <= k
		@variable(m, 0 <= x[1:k] <= 1)
		#@variable(m, x[1:k], Bin)

		## s[i,t] i in S, 0 <= t <= k
		@variable(m, 0 <= s[0:S-1, 0:k] <= 1)

		## y1_1[i,t] / y1_2[i,t]  
		@variable(m, 0 <= y1_1[0:S-1, t in 0:k-1] <= 1)
		@variable(m, 0 <= y1_2[0:S-1, t in 1:k] <= 1)
		## y0_1[i,t] / y0_2[i,t]  
		@variable(m, 0 <= y0_1[0:S-1, t in 0:k-1] <= 1)
		@variable(m, 0 <= y0_2[0:S-1, t in 1:k] <= 1)

		## z[t] 1 <= t <= k
		@variable(m, 0 <= z[0:k] <= 1)

		# CONS: cons_z
		#@constraint(m, cons_sum[i in 0:S-1, j in i+1:S-1, t in kmin:k], s[i,t] + s[j,t] <= 1 + z[t] )
		#lambda
		@variable(m, 0 <= l[0:S-1, kmin:k] <= 1 )
		@constraint(m, cons_l1[t in kmin:k], sum(l[i,t] for i in 0:S-1) <= 1 - z[t])
		@constraint(m, cons_l2[i in 0:S-1, t in kmin:k], l[i,t] >= s[i,t] - z[t])

	else

		## x[t] with 1 <= t <= k
		@variable(m, x[1:k], Bin)

		## s[i,t] i in S, 0 <= t <= k
		@variable(m, s[0:S-1, 0:k], Bin)

		## y1_1[i,t] / y1_2[i,t]  
		@variable(m, y1_1[0:S-1, t in 0:k-1], Bin)
		@variable(m, y1_2[0:S-1, t in 1:k], Bin)
		## y0_1[i,t] / y0_2[i,t]  
		@variable(m, y0_1[0:S-1, t in 0:k-1], Bin)
		@variable(m, y0_2[0:S-1, t in 1:k], Bin)

		## z[t] 1 <= t <= k
		@variable(m, z[0:k], Bin)

		# CONS: cons_z
		@constraint(m, cons_sum[t in kmin:k], sum(s[i,t] for i in 0:S-1) <= 1 + z[t] * (S-1) )

	end



	####################################################### CONS #######################################################

	## s[i,t] = 1 iff i appears at time t

	@constraint(m, cons_s0[i in 0:S-1], s[i,0] == 1)

	#@constraint( m, cons_j1[j in 0:S-1, t in 1:k], s[fsm.succ[j+1][2],t] >= s[j,t-1] + x[t] - 1 )
	#@constraint( m, cons_j0[j in 0:S-1, t in 1:k], s[fsm.succ[j+1][1],t] >= s[j,t-1] - x[t] )

	## P

	@constraint(m, [i in 0:S-1, t in 1:k], y0_2[i,t] <= 1 - x[t])
	@constraint(m, [i in 0:S-1, t in 1:k], y1_2[i,t] <= x[t])

	@constraint(m, [i in 0:S-1, t in 1:k], y0_1[i,t-1] <= y0_2[fsm.succ[i+1][1], t] )
	@constraint(m, [i in 0:S-1, t in 1:k], y1_1[i,t-1] <= y1_2[fsm.succ[i+1][2], t] )

	@constraint(m, [i in 0:S-1, t in 1:k], sum(y0_1[j,t-1] for j in 0:S-1 if fsm.succ[j+1][1]==i) >= y0_2[i,t] )
	@constraint(m, [i in 0:S-1, t in 1:k], sum(y1_1[j,t-1] for j in 0:S-1 if fsm.succ[j+1][2]==i) >= y1_2[i,t] )

	@constraint(m, [t in 1:k], sum(y0_2[i,t] for i in 0:S-1) >= 1 - x[t] )
	@constraint(m, [t in 1:k], sum(y1_2[i,t] for i in 0:S-1) >= x[t] )

	for i in 0:S-1
		bool_i0 = true
		bool_i1 = true
		for j in 0:S-1
			fsm.succ[j+1][1] == i ? bool_i0 = false : nothing
			fsm.succ[j+1][2] == i ? bool_i1 = false : nothing
		end
		bool_i0 ? @constraint(m, [t in 1:k], y0_2[i,t] == 0) : nothing
		bool_i1 ? @constraint(m, [t in 1:k], y1_2[i,t] == 0) : nothing
	end

	@constraint(m, [i in 0:S-1, t in 1:k], s[i,t-1] == y0_1[i,t-1] + y1_1[i,t-1])
	@constraint(m, [i in 0:S-1, t in 1:k], s[i,t] == y0_2[i,t] + y1_2[i,t])


	@constraint(m, sum(s[i-1,k] for i in final ) == 1 )

	# k_min
	@constraint(m, [t in 0:kmin-1], z[t] == 1)

	# z_t decroissant
	@constraint(m, [t in 1:k], z[t] <= z[t-1])


	####################################################### OBJ ####################################################### 

	@objective(m, Min, sum(z[t] for t in 0:k) )


	####################################################### optimize #######################################################

	optimize!(m)
	
	#println(value.(s))
	#println(value.(z))
	#println(value.(x))

	if (termination_status(m) == MOI.OPTIMAL)
		return objective_value(m), kmin
	else
		return -1, kmin
	end

end


######################################## TEST ########################################
#=

	fsm = read_fsm("../data/fsm_hss.fsm")
	#fsm = read_fsm("../data/fsm_n20_1.fsm")

	G = get_graph(fsm)

	if length(attracting_components(G)) == 1
		#strong_conn = attracting_components(G)[1]
		
		#println(min_rayon(fsm, G, strong_conn))
		merg = longest_merging(fsm, strong_conn)
		#println( merg )
		
		#compact( fsm, 6, merg, strong_conn )
		#env_conv( fsm, 6, merg, strong_conn )

	else
		println("More than 1 attracting components: no SS")

	end

	#println(min_rayon(fsm, G))
	#println(longest_merging(fsm))

	#compact( fsm, 6 )

	#min_Sprime(fsm, 4)

	#env_conv( fsm, 6 )

=#