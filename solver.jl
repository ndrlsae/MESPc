"""
Essa função resolve qualquer modelo e devolve as informações definidas em output

#Argumentos
- method: string nomeando o método a ser utilizado
- params: parâmetros necessários para o modelo em questão (matriz de covariância, matriz de adjacência)
- s : tamanho da submatriz selecionada
- gamma: fator de escala
- INEQ : lista de tipos desigualdades a serem adicionadas, padrão: []
- OPTM : Otimizador a ser usado, padrão: SCS.Optmizer
- output: informações a serem retornadas, padrão: ["obj","X_m", "x_v","t_f"]
- eps: tolerância, padrão: 1e-3
"""
function solver(method, params, s, gamma, INEQ = [],
                OPTM = SCS.Optimizer,output = ["obj","X_m", "x_v","t_f"], eps = 1e-3)
    init_t = time()
    title = method
    f_obj = (model) -> error("The objective function was not defined for the method: $method")
    method = uppercase(method)
    #println(method)
    if method == "BQP"
        C, = params;
        model, W = BQP(C,s, gamma, false, OPTM);
        f_obj = (model) -> objective_value(model) - s*log(gamma);
    elseif method == "BQPC"
        C, A = params;
        model, W = BQP(C,s, gamma, A, OPTM);
        f_obj = (model) -> objective_value(model) - s*log(gamma);
    elseif method == "LINX"
        #println(title)
        C, = params
        #println(length(C))
        model, x = LINX(C, s, gamma, OPTM);
        f_obj = (model) -> (1/2)*(objective_value(model) - s*log(gamma));
    elseif method == "LINXC"
        C, A = params
        model, x = LINX(C, s, gamma, OPTM)
        JuMP.optimize!(model); x_v = value.(x);
        improved = true
        while improved
            improved, cut = linx_cuts(x_v, A)
            if improved
                @constraint(model, x'*cut <= 1);
                JuMP.optimize!(model); x_v = value.(x);
            else break end
        end
        f_obj = (model) -> (1/2)*(objective_value(model) - s*log(gamma));
    elseif method == "LINXX"
        C, A = params
        model, W = LINXX(C,s, gamma, A, OPTM)
        f_obj = (model) -> (1/2)*(objective_value(model) - s*log(gamma));
    else println("This method does not exist"); return 0; end
    
    n_cuts0 = num_constraints(model, count_variable_in_set_constraints = false)
    #println("Model $title created with $n_cuts0 constraints")
    
    if method != "LINXC" && method != "LINX"
        X = W[2:end, 2:end]
        x = W[2:end, 1]
        JuMP.optimize!(model); W_v = value.(W);
        X_m = W_v[2:end,2:end];
        x_v = W_v[2:end, 1];
    end
    list_obj = []
    for ineq in INEQ
        title *= " " *ineq
        ineq = uppercase(ineq)
        if ineq == "MCCORMICK_ALL" model = mccormick_ALL(model, X, x);#All(2) McCormick
            
        elseif ineq == "MCCORMICK_TRIANGLE" model = mccormick_triangle(model, X, x);#All(2+3) McCormick
        elseif ineq == "TRIANGLE" model = triangle(model, X, x);#All(3) McCormick
            
        elseif ineq == "INDEPSETS_ALL" model =  indep_set_cuts_ALL(model, X, x, A); #All Indep Sets
        
        elseif ineq == "INDEPSETS1_IT" model, list_obj = INEQ_IT(model, binary_indep_sets, W, A); #Binary Indep Sets
        elseif ineq == "INDEPSETS1_IT_WS" model, list_obj = INEQ_IT(model, binary_warm_start, W, A); #Binary Indep Sets
        elseif ineq == "INDEPSETS1_IT_RLX" model, list_obj = INEQ_IT(model, relax_indep_set, W, A); #Binary Indep Sets
        
        elseif ineq == "INDEPSETS2_IT"  model, list_obj = INEQ_IT(model, greedy12, W, A); #Heuristic Indep Sets 1,2
        elseif ineq == "INDEPSETS11_IT"  model, list_obj  = INEQ_IT(model, greedy_indep_sets, W, A); #Heuristic Indep Sets 1
        
        elseif ineq == "MCCORMICK_IT" model = INEQ_IT_MC(model, W, X_m, x_v, 1e-2, 200, 200) #McCormick Inequalities IT
            #INEQ_IT_MC(MODEL, W_var, X_m, x_v, eps = 1e-2, lim_it = 70, lim_k = 500)
        
        elseif ineq == "MC_INDEP1_IT" || ineq == "INDEP1_MC_IT" #alternating McCormick and Binary Indep Sets 
            model = INDEP_IT_MC(model, binary_indep_sets, W, X_m, A, x_v, 1e-2, 70, 200) #lim_it lim_k
        
        elseif ineq == "MC_INDEP1_IT_WS" || ineq == "INDEP1_MC_IT_WS" #alternating McCormick and Binary Indep Sets 
            model = INDEP_IT_MC(model, binary_warm_start, W, X_m, A, x_v, 1e-2, 70, 200) #lim_it lim_k
####################################################################################################################3        
            #elseif 
        elseif ineq == "MC_INDEP2_IT" || ineq == "INDEP2_MC_IT" #alternating McCormick and Heuristic Indep Sets
            model = INDEP_IT_MC(model, greedy_indep_sets, W, X_m, A, x_v, 1e-2, 500, 200) #lim_it lim_k
############################################################################################################################
        #elseif ineq == "I2_MC_I1" || ineq == "I2_MC_I1" #alternating McCormick and Heuristic Indep Sets + it1
        #    model = INDEP12_IT_MC(model, greedy_indep_sets, W, X_m, A, x_v, 1e-2, 200, 100) #lim_it lim_k
        end
    end

    n_cuts = num_constraints(model, count_variable_in_set_constraints = false)
    
    println(title)
    println("num constraints: ", n_cuts0, " to ", n_cuts)

    inicio = time()
    #println("pre opt")
    JuMP.optimize!(model);
    fim = time()
    fim_t = time()
    delta = fim - inicio
    delta_t = fim_t - init_t
    #println("fim")
    x_v = value.(x)
    objval = f_obj(model)
    println(termination_status(model))
    println(title," ", s, ": f_obj ",objval," em ", delta_t, "s, com tempo final = ", delta, "s")

    out = []
    for info in output #["obj","X_m", "x_v","t_f"]
        if info == "obj" push!(out,objval)
        elseif info == "X_m" X_m = value.(X_m); push!(out,X_m)
        elseif info == "x_v" push!(out,x_v)
        elseif info == "t_f" push!(out,delta)
        elseif info == "t_t" push!(out, delta_t)
        elseif info == "ctr" ctr = all_constraints(model; include_variable_in_set_constraints = false); push!(out,ctr)
        elseif info == "n_ctr" push!(out, n_cuts)
        elseif info == "add_ctr" push!(out, n_cuts-n_cuts0)
        elseif info == "l_obj" push!(out, list_obj)
        #elseif info == "W_m" push!(out, W)
        elseif info == "MODEL" push!(out, model)
        end
    end
    #println("fim solver")
    return out
end



"""
Essa função adiciona cortes de maneira iterativa, dada a função de busca de cortes válidos

#Argumentos
- MODEL: modelo atual
- sets_function: função utilizada para encontrar conjuntos independentes
- W_var: matriz variável W
- A_matrix: matriz de adjacência
- eps: tolerância, padrão: 1e-3
- lim_it: limite de iterações

#Saída: Modelo com cortes adicionados 
        *(+ a lista de valores da função objetivo para estudar a convergência)
"""
function INEQ_IT(MODEL, sets_function, W_var, A_matrix, eps = 1e-3, lim_it = 20)
    n = size(A_matrix,1)
    X = W_var[2:end, 2:end]
    x = W_var[2:end, 1]
    it, J = 0, [0]
    list_obj = []
    while length(J) > 0 && it < lim_it
        JuMP.optimize!(MODEL); W_v = value.(W_var);
        X_m = W_v[2:end,2:end];
        x_v = W_v[2:end, 1];
        obj_it = objective_value(MODEL)
        push!(list_obj, obj_it)
        J = sets_function(X_m, A_matrix, eps)
        #println(length(J), " cortes adicionados na iteração $it")
        it+=1
        for l in J
            #println(l)
            id = l[1]
            x_sub = l[2]
            for i in 1:n if x_sub[i]!=0 && x_sub[i]!=1 println("HÁ X NÃO BINARIO na casa $i ele é igual a ", x_sub[i]) end end
            if id == 0
                x_tilda = [x[i] for i in 1:n if x_sub[i] > 0.99]
                @constraint(MODEL, sum(x_tilda )<= 1)
            else
                X_tilda = [X[id,i] for i in 1:n if x_sub[i]>0.99]
                @constraint(MODEL, sum(X_tilda) <= x[id]);
            end
        end
    end
    #println("Foram feitas $it iterações de IT: ", sets_function)
    return MODEL, list_obj
end


"""
Essa função adiciona cortes de McCormick de maneira iterativa

#Argumentos
- MODEL: modelo atual
- W_var: matriz variável W
- X_m: matriz X atual
- x_v: vetor x atual
- eps: tolerância, padrão: 1e-3
- lim_it: limite de iterações, default: 70
- lim_k: limite de cortes adicionados por iteração, default: 500

#Saída: Modelo com cortes adicionados 
"""
function INEQ_IT_MC(MODEL, W_var, X_m, x_v, eps = 1e-2, lim_it = 70, lim_k = 500)
    delta, it = 1,0
    while delta == 1 && it < lim_it
        delta = 0
        MODEL,delta = mccormick_IT(MODEL,W_var, X_m,x_v, eps, lim_k)
        it+=1
        if it%10 == 0 println("$it constraints: ", num_constraints(MODEL, count_variable_in_set_constraints = false)); end
        #if k > 70 delta = 0 end
        JuMP.optimize!(MODEL);
        W_v = value.(W_var);
        X_m = W_v[2:end,2:end];
        x_v = W_v[2:end, 1];
    end
    return MODEL
end

"""
Essa função adiciona cortes de McCormick e cortes de conjunto independente de maneira iterativa, 
dada a função de busca de cortes válidos

#Argumentos
- MODEL: modelo atual
- sets_function: função utilizada para encontrar conjuntos independentes
- W_var: matriz variável W
- X_m: matriz X atual
- A: matriz de adjacência
- x_v: vetor x atual
- eps: tolerância, padrão: 1e-2
- lim_it: limite de iterações, default: 500
- lim_k: limite de cortes adicionados por iteração, default: 500

#Saída: Modelo com cortes adicionados 
"""
function INDEP_IT_MC(MODEL, sets_function, W_var, X_m, A, x_v, eps = 1e-2, lim_it = 500, lim_k = 500)
    for _ in 1:lim_it
        MODEL,delta = mccormick_IT(MODEL,W_var, X_m,x_v, eps, lim_k) 
        MODEL = INEQ_IT(MODEL, sets_function, W_var, A, eps, 30) #lim_it
        JuMP.optimize!(MODEL);
        W_v = value.(W_var);
        X_m = W_v[2:end,2:end];
        x_v = W_v[2:end, 1];
    end
    return MODEL
end;


function INDEP12_IT_MC(MODEL, sets_function, W_var, X_m, A, x_v, eps = 1e-2, lim_it = 500, lim_k = 500)
    for _ in 1:lim_it
        MODEL,delta = mccormick_IT(MODEL,W_var, X_m,x_v, eps, lim_k) 
        MODEL = INEQ_IT(MODEL, sets_function, W_var, A, eps, 1) #lim_it
        JuMP.optimize!(MODEL);
        W_v = value.(W_var);
        X_m = W_v[2:end,2:end];
        x_v = W_v[2:end, 1];
    end
    MODEL = INEQ_IT(MODEL, binary_warm_start, W, A)
    return MODEL
end;
        
##LIM_K PARA LIMITE CUTS
##LIM_IT PARA LIMITE DE ITERAÇÕES





    