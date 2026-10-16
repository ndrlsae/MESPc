include("data_init63.jl");

include("BQP.jl");
include("linx.jl");
#include("gamma.jl")
include("funcs.jl");
include("cuts.jl");

gamma_new = CSV.read("gamma_new.csv", DataFrame);


matriz_C = matriz_A;

matrizes_adj = [m63_02, m63_04, m63_06, m63_08];
ps = [0.2,0.4,0.6,0.8]; EPS = 1e-5;

OUT = ["obj", "n_ctr", "t_f","t_t"];
OPTMIZER = Mosek.Optimizer

results = [];
open("BQPlinx.txt", "w") do io

    matriz_adj =  matrizes_adj[4]

    println(io, "matriz 63 | usando $OPTMIZER | McCormick IT")#, p); println("p:", p);

    max_c = P_max_clique(matriz_adj);
    println(io, "clique max: ", max_c);
    println("clique max: ", max_c);

    for i in 1:max_c
    println(io, i)
    println(i)
    
    gamma_iX = gamma_new.gamma_iX[i] #compute_gamma_linx(matriz_C,i)#, gamma_init=false, tol = 5e-2,RLT = false)
    gamma_iBQP = gamma_bqp63[i]
    #println(io, "gamma X $i: $gamma_i")
    #println("gamma $i: $gamma_i")
    #gamma_iX = gammas_linxX[i]
    OBJ8 = solverLINX(matriz_C, i, gamma_iX, matriz_adj, true, ["McCormick_IT"], OUT, EPS, OPTMIZER)
    println(io, "linxX $i: $OBJ8")

    OBJ9 = solverBQP(matriz_C, i, gamma_iBQP, matriz_adj, ["it","McCormick"], OUT, EPS, OPTMIZER)
    println(io, "BQP $i: $OBJ9")

    line = (s = i, gammaX = gamma_iX, 
            linxX_obj = OBJ8[1], linxX_ctr = OBJ8[2], linxX_tf = OBJ8[3], linxX_tt = OBJ8[4],
            gammaBQP = gamma_iX,
            BQPc_obj = OBJ9[1], BQPc_ctr = OBJ9[2], BQPc_tf = OBJ9[3], BQPc_tt = OBJ9[4])
    push!(results, line);
    
    end
println(io, "results = ", results)
#df= DataFrame(results)
end
data = DataFrame(results);
CSV.write("linxXBQPc_scs_MC_IT.csv", data);




            