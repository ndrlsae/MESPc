"""
Essa função calcula um lower bound para o MESPc

#Argumentos
- C : matriz de covariância
- A : matriz de adjacência
- s : tamanho da submatriz a ser selecionada

#Saída: valor da função objetivo e tempo de execução
"""
function heuristic(C,A,s)
    n = size(C,1)
    init_t = time()
    MODEL = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(MODEL);
    set_time_limit_sec(MODEL,3600);

    @variable(MODEL, x[1:n], Bin)
    for i in 1:n
        for j in i+1:n
            if A[i,j] == 0
                @constraint(MODEL, x[i] + x[j] <= 1)
            end
        end
    end
    @constraint(MODEL, sum(x) == s)
    
    f_obj = sum([C[i,i]*x[i] for i in 1:n])
    @objective(MODEL, Max, f_obj)
    
    #println("Modelo criado")
    JuMP.optimize!(MODEL);
    objval0 = objective_value(MODEL)
    #Resultados
    if has_values(MODEL)
        x = value.(x) ## provavelmente é sempre igual pois todos produziram a mesma f_obj(C,x)
    else
        println("Sem solução encontrada.")
        x = [0 for i in 1:n]
    end
    x0 = [x[i] for i in 1:n]
    eps = 1e-5
    K = [i for i in 1:n if (1-eps < x[i] < 1+eps)]
    S = [i for i in 1:n if (-eps < x[i] < eps)]
    #Local Search
    function f_obj(matriz_C, X)
        n = size(X)
        #println(n)
        matriz_X = X*X'
        Ct = matriz_C- I 
        matriz_log = Ct.*matriz_X + I
        #objval = logdet(matriz_log); #- s*log(gamma);
        objval, sinal = logabsdet(matriz_log)
        return objval
    end
    #objval0 = f_obj(C,x)
    bestGain = eps
    #println("BestGain")
    while bestGain != 0
        bestGain = 0; worst_k,best_j = 0,0;
        f_i = f_obj(C, x)
        #println(f_i)
        for k in K
            for j in S
                test = [(A[i,j] == 1) for i in K if i != k]
                if minimum(test; init = false)
                    x_i = copy(x);
                    x_i[k], x_i[j] = 0, 1
                    f_j = f_obj(C, x_i)
                    delta_ij = f_j - f_i
                    if delta_ij > bestGain
                        bestGain = delta_ij; 
                        worst_k,best_j = copy(k), copy(j);
                    end
                end
            end
        end
        if bestGain > 0
            setdiff!(K,[worst_k]); setdiff!(S,[best_j]);
            push!(K,best_j); push!(S,worst_k);
            x[worst_k], x[best_j] = 0, 1
        end
        #println(bestGain)
    end
    #println("Is x0 == x_now: ", x0==x)
    objval = f_obj(C,x)
    fim_t = time()
    delta_t = fim_t - init_t
    return objval, delta_t
end





function heuristic2(C,A,s, slv, gamma, OPTM = SCS.Optimizer)
    n = size(C,1)
    init_t = time()
    MODEL = JuMP.Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_silent(MODEL);
    set_time_limit_sec(MODEL,3600);

    @variable(MODEL, x[1:n], Bin)
    for i in 1:n
        for j in i+1:n
            if A[i,j] == 0
                @constraint(MODEL, x[i] + x[j] <= 1)
            end
        end
    end
    @constraint(MODEL, sum(x) == s)

    OUT = ["x_v"]
    opt_sol = solver(slv, [C, A], s, gamma, [], OPTM, OUT)
    #opt_sol = solver(slv, [C, A], s, gamma, ["IndepSets1_IT"], OPTM, OUT)
    x_v = opt_sol[1]
    
    f_obj = sum([x_v[i]*x[i] for i in 1:n])
    @objective(MODEL, Max, f_obj)
    
    #println("Modelo criado")
    JuMP.optimize!(MODEL);
    objval0 = objective_value(MODEL)
    #Resultados
    x = value.(x);
    x0 = [x[i] for i in 1:n]
    eps = 1e-5
    K = [i for i in 1:n if (1-eps < x[i] < 1+eps)]
    S = [i for i in 1:n if (-eps < x[i] < eps)]

    #Local Search
    function f_obj(matriz_C, X)
        n = size(X)
        #println(n)
        matriz_X = X*X'
        Ct = matriz_C- I 
        matriz_log = Ct.*matriz_X + I
        #objval = logdet(matriz_log); #- s*log(gamma);
        objval, sinal = logabsdet(matriz_log)
        return objval
    end
    
    bestGain = eps
    #println("BestGain")
    while bestGain != 0
        bestGain = 0; worst_k,best_j = 0,0;
        f_i = f_obj(C, x)
        #println(f_i)
        for k in K
            for j in S
                test = [(A[i,j] == 1) for i in K if i != k]
                if minimum(test; init = false)
                    x_i = copy(x);
                    x_i[k], x_i[j] = 0, 1
                    f_j = f_obj(C, x_i)
                    delta_ij = f_j - f_i
                    if delta_ij > bestGain
                        bestGain = delta_ij; 
                        worst_k,best_j = copy(k), copy(j);
                    end
                end
            end
        end
        if bestGain > 0
            setdiff!(K,[worst_k]); setdiff!(S,[best_j]);
            push!(K,best_j); push!(S,worst_k);
            x[worst_k], x[best_j] = 0, 1
        end
        #println(bestGain)
    end
    println("Is x0 == x_now: ", x0==x)
    objval = f_obj(C,x)
    fim_t = time()
    delta_t = fim_t - init_t
    return objval,objval0, delta_t
end