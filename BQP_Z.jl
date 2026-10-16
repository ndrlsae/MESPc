"""
Essa função constrói os modelos BQP_Z ou BQPc_Z e devolve o modelo e a matriz variável W

#Argumentos
- C : matriz de covariância
- s : tamanho da submatriz selecionada
- gamma: fator de escala
- A : matriz de adjacência
- OPTM : Otimizador a ser usado, padrão: SCS.Optmizer
"""
function BQP(C,s,gamma, A = false, OPTM = SCS.Optimizer)
    n = size(C,1);
    title = "BQP"
    #Model
    model = JuMP.Model(OPTM);
    set_silent(model);
    set_time_limit_sec(model, 3600)
    #ZZZZZZZZZZZZZZZZZZZ
    set_optimizer_attribute(model, "eps_abs", eps)
    set_optimizer_attribute(model, "eps_rel", eps)
    set_optimizer_attribute(model, "max_iters", max_iters)
    set_optimizer_attribute(model, "verbose", 0)
    #ZZZZZZZZZZZZZZZZZZ
    #Variáveis
    @variable(model, t) 
    @variable(model, W[1:n+1, 1:n+1], PSD);
    X = W[2:end, 2:end]
    x = W[2:end, 1]
    #Constraints
    @constraint(model, W[1,1]==1)
    @constraint(model, sum(x)==s)
    @constraint(model, [i = 1:n], sum(X[i, j] for j in 1:n) == s * x[i])
    @constraint(model, [i = 1:n], X[i, i] == x[i])
    if A != false
        title *= " Clustering"
        @constraint(model, sum(X[i, j] for i in 1:n, j in 1:n if A[i, j] != 0) == s * (s - 1))
    end
    n_cuts0 = num_constraints(model, count_variable_in_set_constraints = false)
    
    M(i, j) = i == j ? gamma * C[i, i] * X[i, i] + 1 - x[i] : gamma * C[i, j] * X[i, j]
    @constraint(model, vcat(t, 1.0, [M(i, j) for j in 1:n for i in 1:j]) in
                       MOI.LogDetConeTriangle(n))
    @objective(model, Max, t)

    @constraint(model, [i = 1:n, j = i+1:n], X[i, j] >= 0)               #ZZZZZZZZZZZZZZZZZZZ
    @constraint(model, [i = 1:n, j = i+1:n; A[i, j] == 0], X[i, j] == 0) #ZZZZZZZZZZZZZZZZZZZ
    
    return model, W
end

"""
Essa função calcula o bound do BQP dada a matriz W_m

#Argumentos
- C : matriz de covariância
- s : tamanho da submatriz selecionada
- gamma: fator de escala
- W_m: matriz com resultado
"""
function BQP_bound(C, s, gamma, W_m)
    X = W_m[2:end,2:end];
    x = diag(X);
    M = gamma .* C .* X .+ Diagonal(1 .- x)
    return logabsdet(Symmetric(Matrix(M)))[1] - s * log(gamma)
end



