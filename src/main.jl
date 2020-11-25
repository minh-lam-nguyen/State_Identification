using Dates

include("SyncTree.jl")


#dir_fsm = "data/SS_50fsm_n"*string(n)*"/"
dir_fsm = "../data/SS_50fsm_n30/"
list_fsm = readdir(dir_fsm)

count_fsm = 0

debut = Dates.now()

for i in list_fsm
	if occursin( ".fsm", i )

		global count_fsm += 1

		fsm = read_fsm(dir_fsm*i)
		seq = parcours_largeur(fsm)

	end
end

fin = Dates.now()

println("total time: ", fin-debut, " / average: ", Dates.value(fin-debut)/count_fsm, " milliseconds (", count_fsm, " fsm)")
