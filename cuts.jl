"""
Essa função adiciona todos os cortes de conjunto independente ao modelo,
dada a matriz de adjacência

#Argumentos
- MODEL: modelo atual
- X_m: matriz X
- x_v: vetor x
- ADJ_m: matriz de adjacência

#Saída: Modelo com cortes adicionados
"""
function indep_set_cuts_ALL(MODEL, X_m, x_v, ADJ_m)
    n = length(x_v)
    if n == 200 max_cliques_comp = CUTS200;
    #elseif n == 63 max_cliques_comp = CUTS63;
    else
        ADJ_comp = matrix_comp(ADJ_m)
        max_cliques_comp = max_cliques(ADJ_comp)
        cliques_2 = []
        for clique in max_cliques_comp
            l_c = length(clique)
            if l_c > 2
                for i in 1:(l_c-1) for j in (i+1):l_c push!(cliques_2,[clique[i],clique[j]]) end end
            end
        end
        append!(max_cliques_comp, cliques_2)
        max_cliques_comp = unique(max_cliques_comp)
    end
    for clique in max_cliques_comp
        clique_INEQ = [x_v[i] for i in clique]
        @constraint(MODEL, sum(clique_INEQ) <= 1)
                    clique_x = [0 for i in 1:n]
                    for c in clique 
                        clique_x[c] = 1; end
                    for i in 1:n
                        @constraint(MODEL, X_m[i,:]'*clique_x <= x_v[i]); end
                end
    for i in 1:n for j in (i+1):n
            @constraint(MODEL, X_m[i,j] <= x_v[i])
            @constraint(MODEL, X_m[j,i] <= x_v[j])
    end end
    return MODEL
end

"""
Essa função encontra os cortes de conjunto independente de maneira exata

#Argumentos
- X_m: matriz X
- ADJ_m: matriz de adjacência
- eps: tolerância, default: 1e-3 
- ws: warm_start (iniciar o problema exato com uma solução heurística), default: false 
- relax: retira a orbigatoriedade de os conjuntos serem dados por vetores BINÁRIOS, default: false

#Saída: conjunto de desigualdades válidas para o modelo atual
"""
function binary_indep_sets(X_m,ADJ_m, eps = 1e-3; ws = false, relax = false)
    n = size(X_m,1)
    x_BQP = diag(X_m);
    J = []
    MODEL = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(MODEL);
    set_time_limit_sec(MODEL, 3600)
    if !relax
        @variable(MODEL, x[1:n], Bin)
    else
        @variable(MODEL, x[1:n] >= 0)
    end
    
    for i in 1:n
        for j in (i+1):n
            if ADJ_m[i,j] == 1
                @constraint(MODEL, x[i] + x[j] <= 1)
            end
        end
    end
    
    if ws
        x_start = greedy_heuristic(ADJ_m, x_BQP, eps)
        #println("x_start 0 = ", x_start)
        for k in x_start 
            #println("Estou no x_start[$k] de 0")
            set_start_value(x[k], 1)
        end
    end
    @objective(MODEL, Max, x_BQP'*x)
    JuMP.optimize!(MODEL);
    x_bin = value.(x);
    obj = objective_value(MODEL)
    desig = []
    if obj > 1 + eps push!(J,[0,x_bin]); push!(desig, x_bin) end
    for i in 1:n
        if ws
            #x_start = greedy_heuristic(ADJ_m, X_m[i,:], eps)
            #for k in x_start 
                #set_start_value(x[k], 1)
                #x[k].Start = 1
            for k in 1:n
                set_start_value(x[k], x_bin[k])
            end
        end
        @objective(MODEL, Max, X_m[i,:]'*x - x_BQP[i])
        JuMP.optimize!(MODEL);
        obj = objective_value(MODEL)
        x_bin = value.(x)
        if obj > eps && !(x_bin in desig) push!(J,[i,x_bin]); push!(desig, x_bin) end
    end
    return J
end


binary_warm_start(x, adj, eps = 1e-3) = binary_indep_sets(x, adj, eps, ws = true)
relax_indep_set(x, adj, eps = 1e-3) = binary_indep_sets(x, adj, eps, ws = false, relax = true)




"""
Essa função resolve o problema do conjunto independente máximo com peso utilizando uma heurística gulosa

#Argumentos
- graph_m: matriz de adjacência 
- w: vetor de pesos
- eps: tolerância, default: 1e-3 
- advanced_search: ativa a busca local com 2 trocas, default: false

#Saída: conjunto independente de maior peso encontrado de maneira heurística
"""
function greedy_heuristic(graph_m, w, eps = 1e-3; advanced_search = false)
    n = size(w)[1]
    S = Int[]
    V_g = [i for i in 1:n]
    l_v = length(V_g)
    G_,w_ = copy(graph_m),copy(w) #1
    dg_G_ = [sum(G_[i,:]) for i in 1:n]
    while l_v > 0 #2
        sigma_v = [w_[i]/(dg_G_[i] + 1) for i in 1:n] #3
        best_v = argmax(sigma_v) #4
        push!(S,best_v) #5
        ngb_v = [i for i in 1:n if G_[best_v, i] == 1]
        for v_rem in [best_v; ngb_v]
            for i in 1:n
                if G_[v_rem, i] == 1
                    dg_G_[i] -= 1
                end
            end
        end
        G_[best_v, :] .= 0
        G_[:, best_v] .= 0
        for ngb in ngb_v
            G_[ngb, :] .= 0
            G_[:, ngb] .= 0
            w_[ngb] = -10
        end
        w_[best_v] = -10 #6
        l_v -= (length(ngb_v) + 1) #6
    end
    improved = true #7
    while improved #8
        improved = false #9
        S_set = Set(S) 
        V_S = [i for i in 1:n if !(i in S_set)] #V(G)S
        #V_S = [i for i in 1:n if !(i in S)] #V(G)/S
        for j in S #10
            N_ = []
            for v in V_S 
                if any(graph_m[v, w] == 1 for w in S if w != j)
                    continue
                end
                push!(N_, v)
            end
            if length(N_) == 0 continue; end
            w_v = [w[v] for v in N_]
            if sum(w_v) <= w[j] continue; end
            #println(w_v)
            S_ = [] #12
            G_x, w_ = graph_m[N_,N_], copy(w) #13 #subgraph induced by N_bar
            l_x0 = length(N_)
            l_x = l_x0
            dg_G_x = [sum(G_x[i,:]) for i in 1:l_x0]
            while l_x > 0 #14
                sigma_x = [w_[N_[i]]/(dg_G_x[i] + 1) for i in 1:l_x0] #15
                best_x = argmax(sigma_x) #16
                best_g_x = N_[best_x] #16
                push!(S_,best_g_x) #17
                ngb_x = [i for i in 1:l_x0 if G_x[best_x, i] == 1] #18
                for v_rem in [best_x; ngb_x]
                    for i in 1:l_x0
                        if G_x[v_rem, i] == 1
                            dg_G_x[i] -= 1
                        end
                    end
                end

                G_x[best_x, :] .= 0
                G_x[:, best_x] .= 0
                for ngb in ngb_x 
                    G_x[ngb, :] .= 0 
                    G_x[:, ngb] .= 0 
                    w_[N_[ngb]] = -10 
                end
                w_[best_g_x] = -10 #18
                l_x -=  (length(ngb_x) + 1) #18
            end
            w_S_ = [w[s] for s in S_] 
            if sum(w_S_) > w[j] #19
                S = setdiff(S, [j]) #20
                append!(S,S_); #20
               improved = true; #21
                break #22
            end
        end
        if !improved && advanced_search
            for j in S
                for l in S
                    if j>=l continue end
                    N_ = [] #27
                    for v in V_S  
                        if any(graph_m[v, w] == 1 for w in S if (w != j && w != l))
                            continue 
                        end
                        push!(N_, v) #27
                    end
                    if length(N_) == 0 continue; end
                    w_v = [w[v] for v in N_]
                    if sum(w_v) <= w[j] + w[l] continue; end
                    S_ = [] #28
                    G_x, w_ = graph_m[N_,N_], copy(w) #29 #subgraph induced by N_bar
                    l_x0 = length(N_)
                    l_x = l_x0
                    dg_G_x = [sum(G_x[i,:]) for i in 1:l_x0] #31
                    while l_x > 0 #30
                        sigma_x = [w_[N_[i]]/(dg_G_x[i] + 1) for i in 1:l_x0] #31
                        best_x = argmax(sigma_x) #32
                        best_g_x = N_[best_x] #32
                        push!(S_,best_g_x) #33
                        ngb_x = [i for i in 1:l_x0 if G_x[best_x, i] == 1] #34
                        for v_rem in [best_x; ngb_x]
                            for i in 1:l_x0
                                if G_x[v_rem, i] == 1
                                    dg_G_x[i] -= 1
                                end
                            end
                        end
                        G_x[best_x, :] .= 0
                        G_x[:, best_x] .= 0
                        for ngb in ngb_x 
                            G_x[ngb, :] .= 0 
                            G_x[:, ngb] .= 0 
                            w_[N_[ngb]] = -10 
                        end
                        w_[best_g_x] = -10 #34
                        l_x -=  (length(ngb_x) + 1) #34
                    end
                    w_S_ = [w[s] for s in S_] #35
                    if sum(w_S_) > w[j] + w[l] #35
                        S = setdiff(S, [j,l]) #36
                        append!(S,S_); #36
                        improved = true; #37
                        break #38
                    end
                end
                if improved break end
            end
        end
    end
    return S
end


"""
Essa função encontra os cortes de conjunto independente de maneira heurpistica

#Argumentos
- X_m: matriz X
- ADJ_m: matriz de adjacência
- eps: tolerância, default: 1e-3 
- advanced_search: ativa a busca local com 2 trocas da função greedy_heuristic, default: false

#Saída: conjunto de desigualdades válidas para o modelo atual
"""
function greedy_indep_sets(X_m, ADJ_m, eps = 1e-3; advanced_search = false)
    n = size(X_m, 1)
    x_v = diag(X_m)
    J = Tuple{Int, Vector{Int}}[] 
    desig = Set{Vector{Int}}() 
    
    for i in 0:n
        if i == 0
            weights = x_v
        else
            weights = X_m[i, :]
        end
        
        S = greedy_heuristic(ADJ_m, weights, eps, advanced_search = advanced_search)
        
        S_v = zeros(Int, n)
        S_v[S] .= 1 
        
        if i == 0 
            z = sum(weights[s] for s in S) - 1.0
        else 
            z = sum(weights[s] - x_v[s] for s in S)
        end
        
        if z > eps && !(S_v in desig) 
            push!(J, (i, S_v)) 
            push!(desig, S_v) 
        end
    end
    
    return J
end

greedy12(x, adj, eps = 1e-3) = greedy_indep_sets(x, adj,  eps, advanced_search = true)



"""
Essa função adiciona cortes de McCormick de maneira iterativa ao modelo

#Argumentos
- MODEL: modelo atual
- W_var: valor atual da matriz X
- X_m: matriz variável X
- x_v: vetor variável x
- eps: tolerância, default: 1e-2 
- lim_k: limite de cortes a serem adicionados, default: 500

#Saída: Modelo com cortes adicionados e um parâmetro indicando se houve adição de novos cortes
"""
function mccormick_IT(MODEL, W_var, X_m, x_v, eps = 1e-2, lim_k = 500) #acho que posso diminuir, porque o IT12ITMC ficou igual com lim_k igual a 200
    X = W_var[2:end, 2:end]
    x = W_var[2:end, 1]
    
    n = length(x_v)
    #DELTA = 0
    violated_cuts = []
    for i in 1:n
        for j in (i+1):n
            v_1 = X_m[i,j] - x_v[i]
            v_2 = X_m[i,j] - x_v[j]
            v_3 = (x_v[i] + x_v[j] - 1) - X_m[i,j]
            
            if v_1 > eps push!(violated_cuts, (v_1, i, j, "UB1")) end
            if v_2 > eps push!(violated_cuts, (v_2, i, j, "UB2")) end
            if v_3 > eps push!(violated_cuts, (v_3, i, j, "LB")) end
        end
    end
    l_vc = length(violated_cuts)
    if l_vc > 0 violated = sort(violated_cuts, rev=true)
    else return MODEL, 0 end
    if l_vc < lim_k lim_k = l_vc end
    for k in 1:lim_k
        vio, i, j, cut_type = violated[k]
        if cut_type == "UB1" @constraint(MODEL, X[i,j]<=x[i])
        elseif cut_type == "UB2" @constraint(MODEL, X[i,j]<=x[j])
        elseif cut_type == "LB" @constraint(MODEL, X[i,j]>= x[i] + x[j] - 1) end
    end
    return MODEL, 1
end


function mccormick_IT2(MODEL, W_var, X_m, x_v, eps = 1e-3)
    X = W_var[2:end, 2:end]
    x = W_var[2:end, 1]
    
    n = length(x_v)
    DELTA = 0
    for i in 1:n
        #println("i = $i")
        if DELTA == 1 break end
        for j in (i+1):n
            #println("i = j")
            if DELTA == 1 break end
            if X_m[i,j] - x_v[i] > eps @constraint(MODEL, X[i,j]<=x[j]);DELTA=1; end
            if X_m[i,j] - x_v[j] > eps @constraint(MODEL, X[i,j]<=x[j]);DELTA=1; end
            if x_v[i] + x_v[j] - 1 - X_m[i,j] > eps @constraint(MODEL, X[i,j]>=x[i]+x[j]-1);DELTA=1; end
            if X_m[i,j]  < -eps @constraint(MODEL, X[i,j]>=0);DELTA=1; end
        end
    end
    return MODEL, DELTA
end


"""
Essa função adiciona todos os cortes de McCormick ao modelo

#Argumentos
- MODEL: modelo atual
- X_m: matriz X
- x_v: vetor x

#Saída: Modelo com cortes adicionados
"""
function mccormick_ALL(MODEL, X_m, x_v)
    n = length(x_v)

    for i in 1:n
        for j in (i+1):n
            @constraint(MODEL, X_m[i,j]<=x_v[i])
            @constraint(MODEL, X_m[i,j]<=x_v[j])
            @constraint(MODEL, X_m[i,j]>=x_v[i]+x_v[j]-1)
            @constraint(MODEL, X_m[i,j]>=0)
        end
    end
    return MODEL
end

"""
Essa função adiciona todos os cortes Triangulares ao modelo

#Argumentos
- MODEL: modelo atual
- X_m: matriz X
- x_v: vetor x

#Saída: Modelo com cortes adicionados
"""
function triangle(MODEL, X_m, x_v)
    n = length(x_v)
    for i in 1:n
        for j in (i + 1):n
            #@constraint(MODEL, X_m[i,j]<=x_v[i])
            #@constraint(MODEL, X_m[i,j]<=x_v[j])
            #@constraint(MODEL, X_m[i,j]>=x_v[i]+x_v[j]-1)
            #@constraint(MODEL, X_m[i,j]>=0)
            for k in (j + 1):n
                @constraint(MODEL,X_m[i,j] - X_m[i,k] - X_m[j,k] >= - x_v[k])
                @constraint(MODEL, -X_m[i,j] + X_m[i,k] - X_m[j,k] >= - x_v[j])
                @constraint(MODEL, -X_m[i,j] - X_m[i,k] + X_m[j,k] >= - x_v[i])
                @constraint(MODEL,X_m[i,j] + X_m[i,k] + X_m[j,k] >= x_v[i] + x_v[j] + x_v[k] - 1)
            end
        end
    end
    return MODEL
end



"""
Essa função adiciona todos os cortes de McCormick e Triangulares ao modelo

#Argumentos
- MODEL: modelo atual
- X_m: matriz X
- x_v: vetor x

#Saída: Modelo com cortes adicionados
"""
function mccormick_triangle(MODEL, X_m, x_v)
    n = length(x_v)
    for i in 1:n
        for j in (i + 1):n
            @constraint(MODEL, X_m[i,j]<=x_v[i])
            @constraint(MODEL, X_m[i,j]<=x_v[j])
            @constraint(MODEL, X_m[i,j]>=x_v[i]+x_v[j]-1)
            @constraint(MODEL, X_m[i,j]>=0)
            for k in (j + 1):n
                @constraint(MODEL,X_m[i,j] - X_m[i,k] - X_m[j,k] >= - x_v[k])
                @constraint(MODEL, -X_m[i,j] + X_m[i,k] - X_m[j,k] >= - x_v[j])
                @constraint(MODEL, -X_m[i,j] - X_m[i,k] + X_m[j,k] >= - x_v[i])
                @constraint(MODEL,X_m[i,j] + X_m[i,k] + X_m[j,k] >= x_v[i] + x_v[j] + x_v[k] - 1)
            end
        end
    end
    return MODEL
end


"""
Essa função encontra os cortes de conjunto independente do modelo linxc

#Argumentos
- X_V: matriz variável X
- ADJ_m: matrzi de adjacência
- eps: tolerância, default = 1e-3

#Saída: corte
"""
function linx_cuts(X_V, ADJ_m, eps = 1e-3)
    n = size(X_V,1)
    J = []
    
    MODEL = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(MODEL);
    set_time_limit_sec(MODEL, 3600)

    @variable(MODEL, x[1:n], Bin)
     for i in 1:n
        for j in (i+1):n
            if ADJ_m[i,j] == 1
                @constraint(MODEL, x[i] + x[j] <= 1)
            end
        end
    end
    @objective(MODEL, Max, X_V'*x - 1)
    JuMP.optimize!(MODEL);
    x_cuts = value.(x);
    obj = objective_value(MODEL)
    if obj > eps
        return true, x_cuts
    else
        return false, false
    end
end
