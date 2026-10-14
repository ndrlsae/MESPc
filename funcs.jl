#gerador de Matrizes(Grafos) Aleatórias
function NP_matrix(n,p)
    #[i+j for i in 1:3, j in 1:2]
    matrix = [0 for i in 1:n, j in 1:n]
    for i in 1:n
        for j in (i+1):n
            prob = rand()
            if prob < p
                matrix[i,j] = 1 
                matrix[j,i] = 1
            end
        end
    end 
    return matrix
end

#clique maximo
function max_clique(ADJ_M)
    GRAPH = SimpleGraph(ADJ_M)
    println("Grafo criado")
    CLIQUES = maximal_cliques(GRAPH)
    println("Cliques encontradas")
    MAX_CLIQUE = maximum(length,CLIQUES)
    cliques_max = [cli for cli in CLIQUES if length(cli) == MAX_CLIQUE]
    println("Existem $(length(cliques_max)) cliques de tamanho $MAX_CLIQUE nesse grafo")
    return MAX_CLIQUE, length(cliques_max)
end

function P_max_clique(ADJ_M)
    n = size(ADJ_M,1)

    MODEL = JuMP.Model(Mosek.Optimizer);
    set_silent(MODEL);
    set_time_limit_sec(MODEL,3600);

    @variable(MODEL, x[1:n], Bin)
    for i in 1:n
        for j in i+1:n
            if ADJ_M[i,j] == 0
                @constraint(MODEL, x[i] + x[j] <= 1)
            end
        end
    end
    @objective(MODEL, Max, sum(x))
    #println("Modelo criado")
    JuMP.optimize!(MODEL);
    
    #Resultados
    x = value.(x);
    #println("X = ", x)
    objval = sum(x)
    return round(Int,objval)
end

function matrix_comp(ADJ_M)
    n = size(ADJ_M,1)
    COMP_M = [0 for i in 1:n, j in 1:n]
    for i in 1:n
        for j in (i+1):n
            if ADJ_M[i,j] == 0
                COMP_M[i,j] = 1
                COMP_M[j,i] = 1
            #elseif ADJ_M[i,j] == 1
            #    COMP_M[i,j] = 0
            end
        end
    end
    return COMP_M
end

function max_cliques(ADJ_M)
    GRAPH = SimpleGraph(ADJ_M)
    CLIQUES = maximal_cliques(GRAPH)
    MAX_CLIQUE = maximum(length,CLIQUES)
    cliques_max = [cli for cli in CLIQUES if length(cli) == MAX_CLIQUE]
    return CLIQUES
end


function lower_triangular_solver(A)
    return [A[i, j] for i = 1:size(A, 1) for j = 1:i]
end


