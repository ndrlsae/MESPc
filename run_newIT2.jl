include("data_init63.jl");

include("it2_new.jl");

# ---------------------------------------------------------------------------------------
# Example
# ---------------------------------------------------------------------------------------
using MAT, Random
C = Matrix{Float64}(matread("data63.mat")["A"])        # your file stores C under "A"
 #n, p, s = 63, 0.8, 8
 #Random.seed!(4)                                           # Julia's RNG ≠ numpy's: different graph
 #U = triu(rand(n, n) .< p, 1); G = Int.(U .| U')           # adjacency of G(n,p)
 #γ = 6.0757                                                # optimal γ for BQPc at this s
matriz_C = matriz_A;

matrizes_adj = [m63_02, m63_04, m63_06, m63_08];
matriz_adj =  matrizes_adj[4]
for s in 1:14
    println("S = ", s)
    gamma_i = gammas_bqpc63[s]
    res = it2_new(C, matriz_adj, s, gamma_i)
    println("Resultado da iteração s = $s :", res)
end