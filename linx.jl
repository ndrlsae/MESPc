function LINX(C, s, gamma, OPTM = SCS.Optimizer)
    n = size(C,1)
    title = "linx"
    #Model
    model = JuMP.Model(OPTM);
    set_silent(model);
    set_time_limit_sec(model,3600)
    
    #set_optimizer_attribute(model, "max_iters", 2000) ###############################
    #set_optimizer_attribute(model, "acceleration_lookback", 0) # 0 desativa a aceleração ###############################3
    #set_optimizer_attribute(model, "scale", 1.0)
    #set_optimizer_attribute(model, "adaptive_scale", 1)

    #set_optimizer_attribute(model, "eps_abs", 1e-3)
    #set_optimizer_attribute(model, "eps_rel", 1e-3)
    #Variables
    @variable(model, 0 <= x[1:n] <= 1);
    @variable(model, t);
    #Constraints
    @constraint(model, sum(x)==s) # eTx=s

    #println("gamma :$gamma e é ",typeof(gamma))
    E = gamma*C*diagm(x)*C + Matrix(1.0I,n,n) - diagm(x) 

    @constraint(model, [t; 1; JuMP.triangle_vec(E)] in MOI.LogDetConeTriangle(n))
    @objective(model, Max, t)
    #println(title)
    return model, x
end

function LINXX(C, s, gamma, A, OPTM = SCS.Optimizer)
    n = size(C,1)
    title = "linxX"
    #Model
    model = JuMP.Model(OPTM);
    set_silent(model);
    set_time_limit_sec(model,3600)
    #set_optimizer_attribute(model, "max_iters", 2000) #######################################
    #set_optimizer_attribute(model, "acceleration_lookback", 0) # 0 desativa a aceleração ###############################3
    #set_optimizer_attribute(model, "scale", 1.0)
    #set_optimizer_attribute(model, "adaptive_scale", 1)
    
    #set_optimizer_attribute(model, "eps_abs", 1e-3)
    #set_optimizer_attribute(model, "eps_rel", 1e-3)
    #Variables
    @variable(model, t)
    @variable(model, W[1:n+1, 1:n+1], PSD);
    X = W[2:end, 2:end]
    x = W[2:end, 1]
    #Constraints
    @constraint(model, W[1,1]==1)
    @constraint(model, sum(x)==s) # eTx = s
    @constraint(model, X*ones(n) == s*x) # Xe=sx
    total = 0
    for i in 1:n
        @constraint(model, x[i]==X[i,i])
        t1 = dot(A[i,:],X[i,:]) # <A,X> = s(s-1)
        total += t1
    end
    @constraint(model, total == s*(s-1))
    
    E = gamma*C*diagm(x)*C + Matrix(1.0I,n,n) - diagm(x) 
    @constraint(model, [t; 1; JuMP.triangle_vec(E)] in MOI.LogDetConeTriangle(n))
    @objective(model, Max, t)
    #println(title)
    return model, W
end


function solverLINX(C,s,gamma, A = false, X_M = false, INEQ = [false] ,
                    output = ["obj","x_v","t_f"], eps = 1e-3, OPTM = SCS.Optimizer)
    n = size(C,1);
    #println("Usando $OPTM")
    #model
    init_t = time()
    if X_M == false
        model, x = LINX(C,s,gamma, OPTM)
        title = "linx"
        n_cuts0 = num_constraints(model, count_variable_in_set_constraints = false)
        if A != false
            title *= "-c"
            JuMP.optimize!(model); x_v = value.(x);
            improved = true
            while improved
                improved, cut = linx_cuts(x_v, A)
                if improved
                    @constraint(model, x'*cut <= 1)
                    JuMP.optimize!(model); x_v = value.(x);
                else break
                end
            end
        end
    else
        model, W = LINXX(C,s,gamma, A, OPTM)
        X = W[2:end, 2:end]
        x = W[2:end, 1]
        title = "linxX-c "
        n_cuts0 = num_constraints(model, count_variable_in_set_constraints = false)
        if INEQ[1] != false
            title *= INEQ[1]
            JuMP.optimize!(model); W_v = value.(W);
            X_m = W_v[2:end,2:end];
            x_v = W_v[2:end, 1];
            if INEQ[1] == "McCormick_ALL" #|| INEQ[2] == "McCormick_ALL"
                title *= " + MC_ALL"
                sets_function = (_,_,_) -> [] 
                model = mccormick_ALL(model, X, x);
            end
            if INEQ[1] == "Indep_Set" sets_function = binary_indep_sets; end
            if INEQ[1] == "Indep_Set2" sets_function = greedy_indep_sets2; end
            if INEQ[1] == "McCormick_IT" 
                sets_function = (_,_,_) -> [] 
                delta, k = 1,0
                while delta == 1
                    delta = 0
                    model,delta = mccormick_IT(model,W, X_m,x_v)
                    k+=1
                    if k%10 == 0 println("$k constraints: ", num_constraints(model, count_variable_in_set_constraints = false)); end
                    if k > 70 delta = 0 end
                    JuMP.optimize!(model);
                    W_v = value.(W);
                    X_m = W_v[2:end,2:end];
                    x_v = W_v[2:end, 1];
                end end 
            J = sets_function(X_m,A, eps)
             while length(J) > 0
                #println(length(J))
                for l in J
                    #println(l)
                    id = l[1]
                    x_sub = l[2]
                    if id == 0
                        x_tilda = [x[i] for i in 1:n if x_sub[i] > 0.99]
                        @constraint(model, sum(x_tilda )<= 1)
                    else
                        X_tilda = [X[id,i] for i in 1:n if x_sub[i]>0.99]
                        @constraint(model, sum(X_tilda) <= x[id]);
                    end
                end
                JuMP.optimize!(model); W_v = value.(W);
                X_m = W_v[2:end,2:end];
                x_v = W_v[2:end, 1];
                J = sets_function(X_m,A)
            end
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
    objval = (1/2)*(objective_value(model) - s*log(gamma))
    println(termination_status(model))
    println(title," ", s, ": f_obj ",objval," em ", delta_t, "s, com tempo final = ", delta, "s")

    out = []
    for info in output #["obj","X_m", "x_v","t_f"]
        if info == "obj" push!(out,objval)
        #elseif info == "X_m" push!(out,X)
        elseif info == "x_v" push!(out,x_v)
        elseif info == "t_f" push!(out,delta)
        elseif info == "t_t" push!(out, delta_t)
        elseif info == "ctr" ctr = all_constraints(model; include_variable_in_set_constraints = false); push!(out,ctr)
        elseif info == "n_ctr" push!(out, n_cuts)
        elseif info == "add_ctr" push!(out, n_cuts-n_cuts0)
        #elseif info == "W_m" push!(out, W)
        elseif info == "MODEL" push!(out, model)
        end
    end
    return out
end