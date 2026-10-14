function BQP(C,s,gamma, A = false, OPTM = SCS.Optimizer)
    n = size(C,1);
    title = "BQP"
    #Model
    model = JuMP.Model(OPTM);
    set_silent(model);
    set_time_limit_sec(model, 3600)
    #set_optimizer_attribute(model, "FeasibilityTol", 1e-5)
    #set_attribute(model, "eps_abs", 1e-5) 
    #set_attribute(model, "primal_feasibility_tolerance", 1e-5)
    #Variáveis
    @variable(model, t) 
    @variable(model, W[1:n+1, 1:n+1], PSD);
    X = W[2:end, 2:end]
    x = W[2:end, 1]
    #Constraints
    @constraint(model, W[1,1]==1)
    @constraint(model, sum(x)==s)
    @constraint(model, X*ones(n) == s*x)
    for i in 1:n
        @constraint(model, x[i]==X[i,i])
    end

    if A != false
        title *= " Clustering"
        total = 0
        for i in 1:n
        #print("$i ")
            t1 = dot(A[i,:],X[i,:])
            total += t1
        end
        #println(total)
        @constraint(model, total == s*(s-1))
    end
    n_cuts0 = num_constraints(model, count_variable_in_set_constraints = false)
    
    #@constraint(model, X .>= 1e-5*I)
        
    E = gamma*C.*X + Matrix(1.0I,n,n) - diagm(x) 

    @constraint(model, [t; 1; lower_triangular_solver(E)] in MOI.LogDetConeTriangle(n))
    @objective(model, Max, t)
    #println(title)
    return model, W
end


function solverBQP(C,s,gamma, A = false, Ineq = [false], 
                        output = ["obj","X_m", "x_v","t_f"], eps = 1e-3, OPTM = SCS.Optimizer)
    n = size(C,1);
    #println("Usando $OPTM")
    #model
    init_t = time()
    model, W = BQP(C,s,gamma, A, OPTM)
    title = "BQP"
    if A != false title*= "c" end
    X = W[2:end, 2:end]
    x = W[2:end, 1]
    
    n_cuts0 = num_constraints(model, count_variable_in_set_constraints = false)
    #println("num initial of constraints: ", n_cuts0)

    if Ineq[1] == "all"
        for ineq in Ineq[2:end]
            if ineq == "McCormick"
                title *= " McCormick ALL"
                mccormick_ALL(model, X, x)
            elseif ineq == "Indep_Set"
                title *= " IndepSet ALL"
                indep_set_cuts_ALL(model, X, x, A)
            end
        end
    elseif Ineq[1] == "it"
        if A == false println("Esse método não pode ser aplicado sem o Clustering"); end;
        
        JuMP.optimize!(model);
        W_v = value.(W);
        X_m = W_v[2:end,2:end];
        x_BQP = W_v[2:end, 1];
        #x_BQP = diag(X_m);
        if length(Ineq) >= 3
            println("ERROR: escolha apenas 1 procedimento iterativo")
        elseif Ineq[2] == "McCormick"
            title *= " McCormick IT"
            delta, k = 1, 0
            while delta == 1
                delta = 0
                model,delta = mccormick_IT(model,W, X_m,x_BQP)
                k+=1
                if k%10 == 0 println("$k constraints: ", num_constraints(model, count_variable_in_set_constraints = false)); end
                if k > 70 delta=0 end
                JuMP.optimize!(model);
                W_v = value.(W);
                X_m = W_v[2:end,2:end];
                x_BQP = W_v[2:end, 1];
                #x_BQP = diag(X_m);
            end
        elseif Ineq[2] == "Indep_Set" || Ineq[2] == "Indep_Set21" || Ineq[2] == "Indep_Set2"
                title *= " " * Ineq[2] * " IT"
                if Ineq[2] == "Indep_Set" sets_function = binary_indep_sets; end
                if Ineq[2] == "Indep_Set21" sets_function = greedy_indep_sets; end
                if Ineq[2] == "Indep_Set2" sets_function = greedy_indep_sets2; end
            J = sets_function(X_m,A, eps)
            #println(length(J))
            while length(J) > 0
                #println(length(J))
                for l in J
                    #println(l)
                    id = l[1]
                    x_sub = l[2]
                    x_1 = [i for i in 1:36 if 1 - eps < x_sub[i] < 1+ eps]
                    #println(id," : ",x_1)
                    if id == 0
                        x_tilda = [x[i] for i in 1:n if x_sub[i] > 0.99]
                        #println(x_tilda)
                        @constraint(model, sum(x_tilda )<= 1)
                        #@constraint(model, x'*x_sub <= 1)
                    else
                        X_tilda = [X[id,i] for i in 1:n if x_sub[i]>0.99]
                        #println(X_tilda)
                        @constraint(model, sum(X_tilda) <= x[id]);
                        #@constraint(model, X[id,:]'*x_sub <= x[id]);
                    end
                end
                JuMP.optimize!(model);
                W_v = value.(W);
                X_m = W_v[2:end,2:end];
                x_BQP = W_v[2:end, 1];
                #x_BQP = diag(X_m);
                J = sets_function(X_m,A)
            end
        end
    end
    #println("TESTE")
    n_cuts = num_constraints(model, count_variable_in_set_constraints = false)
    
    println(title)
    println("num constraints: ", n_cuts0, " to ", n_cuts)
    
    #Resolvendo
    inicio = time()
    JuMP.optimize!(model);
    fim = time()
    fim_t = time()
    delta = fim - inicio
    delta_t = fim_t - init_t
    W = value.(W);
    X = W[2:end,2:end];
    x = diag(X);
    W = [1 x' ; x X];
    Ct = [0 zeros(1,n);zeros(n) (gamma*C-Matrix(1.0I,n,n))];
    objval = objective_value(model) - s*log(gamma);
    println(termination_status(model))
    println(title," ", s, ": f_obj ",objval," em ", delta_t, "s, com tempo final = ", delta, "s")
    out = []
    for info in output #["obj","X_m", "x_v","t_f"]
        if info == "obj" push!(out,objval)
        elseif info == "X_m" push!(out,X)
        elseif info == "x_v" push!(out,x)
        elseif info == "t_f" push!(out,delta)
        elseif info == "t_t" push!(out, delta_t)
        elseif info == "ctr" ctr = all_constraints(model; include_variable_in_set_constraints = false); push!(out,ctr)
        elseif info == "n_ctr" push!(out, n_cuts)
        elseif info == "add_ctr" push!(out, n_cuts-n_cuts0)
        elseif info == "W_m" push!(out, W)
        elseif info == "MODEL" push!(out, model)
        end
    end
    return out
end
