MAX_ITER = 5;
#optimizer = optimizer_with_attributes(SCS.Optimizer, 
#    "eps_abs" => 1e-4,       # Relaxa a tolerância absoluta (default geralmente é 1e-5 ou 1e-4)
#    "eps_rel" => 1e-4,       # Relaxa a tolerância relativa
#    "normalize" => 1,        # Força o SCS a tentar normalizar os dados internamente
#    "acceleration_lookback" => 0 # Desliga o Anderson Acceleration, que às vezes causa instabilidade
#)
using Clarabel
optimizer = Clarabel.Optimizer

#using Hypatia
#optimizer = Hypatia.Optimizer

function info_subprob_gamma_bqp(X,x,C,s,gamma)
    n = size(C,1); 
    AUX1 = C.*X;
    #println(typeof(x))
    #Compute F1(gamma1,X)
    F1 = gamma*AUX1 - diagm(x) + Matrix(I,n,n);
    #println("F1: ", typeof(F1))
    #Compute v1(gamma1,X)
    v1 = logdet(F1) - s*log(gamma);

    # Compute inv(F1(gamma1,Xstar))
    invF1 = inv(F1);
    #println("res abaixo")
    res = tr(invF1*diagm(ones(n)-x)) - n + s;
    return AUX1,invF1,res
end


function compute_gamma_bqp(C,s, A=false, gamma_init=false, tol = 5e-2,RLT = false)
    n = size(C,1); iter = 0; max_iter = MAX_ITER; 
    # ; # Tolerance for residual and duality gap  
    # Select initial gamma
    c_sort=sort(diag(C),rev=true);
    if gamma_init== false
        global gamma=1/(c_sort[s]);
    else
        global gamma=gamma_init
    end
    # solve the linx relaxation for gamma and obtain x
    all_bound = [];
    all_gamma = [];
    all_res = [];
    iter = 0; dual_gap = Inf; inf_pnorm = Inf;
    mu=1;
    taumu=0.1;
    taudir=0.9;
    
    
    iter = 0;

    OUT = ["obj","X_m","x_v","t_f"]
    while true
            if A == false
                println("Solver MESP")
                #println("length solver : ", length(solver("BQP", [matriz_C], s, gamma, [], SCS.Optimizer, OUT)))
                bound,X, x, time = solver("BQP", [matriz_C], s, gamma, [], optimizer, OUT)
                #println("EAI")
                #bound,X,x, time = solverBQP(C,s,gamma)#,RLT)
            else
                println("Solver MESP Clustering")
                bound,X, x, time = solver("BQPc", [matriz_C, A], s, gamma, [], optimizer, OUT)
                #bound,X,x, time = solverBQP(C,s,gamma,A)#,RLT)
            end
        #println("!!")
        AUX1,invF1,res = info_subprob_gamma_bqp(X,x,C,s,gamma);
        #println("SEI LA")
        push!(all_bound,bound); push!(all_res,res); push!(all_gamma,gamma);
        # Compute MAT= H_G(gamma:=exp(psi))
        Mat = -gamma*dot((ones(n)-x),diag(invF1*AUX1*invF1));
        # Compute the search direction for Newton's method (dir=-res/Mat)
        dir = -res/Mat; edir=exp(dir);
        #Ratio test
        #println("#Ratio test")
        beta = 10^6;
        if dir < 0
            aux = -gamma/dir;
            if aux < beta
                beta = aux;
            end
        end
        #println("first if")
        if beta > 1.05
            beta = 1.;
        else
            beta = taudir * beta;
        end
        #println("if2")
        #gamma1= gamma1+ beta*dir;
        gamma = gamma * exp(beta*dir)
        # Update mu
        mu=taumu*mu
        iter += 1;
        if iter >= max_iter; break;end
        if abs(res) <= tol; break; end
    end
    #println("out of while")
    idx=argmin(all_bound)
    bound = all_bound[idx]; gamma = all_gamma[idx]; res = all_res[idx];
    #println(all_res)
    return gamma
end



function info_subprob_gamma_linx(x,C,s,gamma)
    n = size(C,1); 
    AUX1 = C*diagm(x)*C;

    #Compute F1(gamma1,X)
    F1 = gamma*AUX1 - diagm(x) + Matrix(I,n,n);

    #Compute v1(gamma1,X)
    v1 = (1/2)*(logdet(F1) - s*log(gamma));

    # Compute inv(F1(gamma1,Xstar))
    invF1 = inv(F1);

    res = tr(invF1*diagm(x - ones(n))) + n - s;
    return AUX1,invF1,res
end

function compute_gamma_linx(C,s, A=false, gamma_init=false, tol = 5e-2,RLT = false)
    n = size(C,1); iter = 0; max_iter = MAX_ITER; 
    # ; # Tolerance for residual and duality gap  
    # Select initial gamma
    c_sort=sort(diag(C),rev=true);
    if gamma_init== false
        global gamma=1/(c_sort[s]^2); #### ALTERADO/DIFERENTE DO BQP
    else
        global gamma=gamma_init
    end
    # solve the linx relaxation for gamma and obtain x
    all_bound = [];
    all_gamma = [];
    all_res = [];
    iter = 0; dual_gap = Inf; inf_pnorm = Inf;
    mu=1;
    taumu=0.1;
    taudir=0.9;
    
    
    iter = 0;

    OUT = ["obj", "x_v","t_f"]
    while true
            if A == false
                println("Solver linx")
                bound,x, time = solver("linx", [matriz_C], s, gamma, [], optimizer, OUT)
                #bound,x, time = solverLINX(C,s,gamma)#,RLT)
            else
                println("Solver linxX")
                bound,x, time = solver("linxX", [matriz_C, A], s, gamma, [], optimizer, OUT)
                #bound,x, time = solverLINX(C,s,gamma,A, true)#,RLT)
            end
        AUX1,invF1,res = info_subprob_gamma_linx(x,C,s,gamma);
        push!(all_bound,bound); push!(all_res,res); push!(all_gamma,gamma);
        # Compute MAT= H_G(gamma:=exp(psi))
        Mat = -gamma*dot((ones(n)-x),diag(invF1*AUX1*invF1));
        # Compute the search direction for Newton's method (dir=-res/Mat)
        dir = res/Mat; edir=exp(dir);
        #Ratio test
        #println("#Ratio test")
        beta = 10^6;
        if dir < 0
            aux = -gamma/dir;
            if aux < beta
                beta = aux;
            end
        end
        #println("first if")
        if beta > 1.05
            beta = 1.;
        else
            beta = taudir * beta;
        end
        #println("if2")
        #gamma1= gamma1+ beta*dir;
        gamma = gamma * exp(beta*dir)
        # Update mu
        mu=taumu*mu
        iter += 1;
        if iter >= max_iter; break;end
        if abs(res) <= tol; break; end
    end
    idx=argmin(all_bound)
    bound = all_bound[idx]; gamma = all_gamma[idx]; res = all_res[idx];
    #println(all_res)
    return gamma#, bound
end



function compute_gamma_linxc(C,s, A=false, gamma_init=false, tol = 5e-2,RLT = false)
    n = size(C,1); iter = 0; max_iter = MAX_ITER; 
    # ; # Tolerance for residual and duality gap  
    # Select initial gamma
    c_sort=sort(diag(C),rev=true);
    if gamma_init== false
        global gamma=1/(c_sort[s]^2); #### ALTERADO/DIFERENTE DO BQP
    else
        global gamma=gamma_init
    end
    # solve the linx relaxation for gamma and obtain x
    all_bound = [];
    all_gamma = [];
    all_res = [];
    iter = 0; dual_gap = Inf; inf_pnorm = Inf;
    mu=1;
    taumu=0.1;
    taudir=0.9;
    
    
    iter = 0;

    OUT = ["obj", "x_v","t_f"]
    while true
            if A == false
                println("Solver linx")
                bound,x, time = solver("linx", [matriz_C], s, gamma, [], optimizer, OUT)
                #bound,x, time = solverLINX(C,s,gamma)#,RLT)
            else
                println("Solver linxc")
                bound,x, time = solver("linxc", [matriz_C, A], s, gamma, [], optimizer, OUT)
                #bound,x, time = solverLINX(C,s,gamma,A, true)#,RLT)
            end
        AUX1,invF1,res = info_subprob_gamma_linx(x,C,s,gamma);
        push!(all_bound,bound); push!(all_res,res); push!(all_gamma,gamma);
        # Compute MAT= H_G(gamma:=exp(psi))
        Mat = -gamma*dot((ones(n)-x),diag(invF1*AUX1*invF1));
        # Compute the search direction for Newton's method (dir=-res/Mat)
        dir = res/Mat; edir=exp(dir);
        #Ratio test
        #println("#Ratio test")
        beta = 10^6;
        if dir < 0
            aux = -gamma/dir;
            if aux < beta
                beta = aux;
            end
        end
        #println("first if")
        if beta > 1.05
            beta = 1.;
        else
            beta = taudir * beta;
        end
        #println("if2")
        #gamma1= gamma1+ beta*dir;
        gamma = gamma * exp(beta*dir)
        # Update mu
        mu=taumu*mu
        iter += 1;
        if iter >= max_iter; break;end
        if abs(res) <= tol; break; end
    end
    idx=argmin(all_bound)
    bound = all_bound[idx]; gamma = all_gamma[idx]; res = all_res[idx];
    #println(all_res)
    return gamma#, bound
end
