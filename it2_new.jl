#=
New IT2 for MESPc: cutting-plane loop (Algorithm 1 of the paper) on the BQPc relaxation,
with the improved heuristic separation of independent-set inequalities.

Changes vs. the paper's IT2 (Algorithm 2 applied once per l):
  (a) support reduction: heuristic runs only on {i : w_i > tol}; skip l if even the whole
      support cannot give a violation; skip rows l with x_l ~ 0
  (b) GWMIN2 greedy rule  w_i / w(N[i])  instead of  w_i / (deg+1)
  (c) incremental best-improvement 1-exchange using cnt[i] = #neighbours of i in S
  (d) every violated set visited by the search is a candidate cut; up to k_per most
      violated distinct sets are added per l
  (e) a pool of all sets found so far is scored for every l with one matrix product

Cuts:  l = 0  ->  x' v <= 1          l >= 1  ->  X[l,:] v <= x[l]
(v = characteristic vector of an independent set of G)

Packages: JuMP, SCS, LinearAlgebra (+ MAT, Random for the example at the bottom)
=#

using JuMP, SCS, LinearAlgebra

# ---------------------------------------------------------------------------------------
# Heuristic separation
# ---------------------------------------------------------------------------------------

"""GWMIN2 greedy on the vertices `cand` of graph `Ab` (Bool adjacency) with weights `w`."""
function greedy2(Ab::AbstractMatrix{Bool}, w::Vector{Float64}, cand::Vector{Int})
    S = Int[]
    alive = trues(length(cand))
    while any(alive)
        c = cand[alive]
        nbw = Ab[c, c] * w[c] .+ w[c]                    # w(N[i]) in the remaining graph
        k = argmax(w[c] ./ max.(nbw, 1e-15))
        i = c[k]
        push!(S, i)
        alive .&= .!Ab[i, cand]                          # remove neighbours of i
        alive[cand .== i] .= false                       # remove i
    end
    return S
end

"""Greedy + incremental best-improvement 1-exchange.
Returns all distinct independent sets visited (sorted index vectors)."""
function heuristic_sets(Ab::AbstractMatrix{Bool}, w::Vector{Float64})
    m = length(w)
    S = greedy2(Ab, w, collect(1:m))
    visited = Set{Vector{Int}}([sort(S)])
    inS = falses(m)
    inS[S] .= true
    while true
        cnt = vec(sum(Ab[:, inS], dims = 2))             # neighbours of each vertex in S
        bestgain, bestj, bestSp = 1e-12, 0, Int[]
        Sidx = findall(inS)
        for j in Sidx
            # vertices that become free when j leaves S: outside S, only neighbour in S is j
            Np = findall(.!inS .& (cnt .== 1) .& Ab[:, j])
            isempty(Np) && continue
            Sp = greedy2(Ab, w, Np)
            gain = sum(w[Sp]) - w[j]
            push!(visited, sort(vcat(filter(!=(j), Sidx), Sp)))
            if gain > bestgain
                bestgain, bestj, bestSp = gain, j, Sp
            end
        end
        bestj == 0 && return collect(visited)
        inS[bestj] = false
        inS[bestSp] .= true
    end
end

"""Separation round. Returns a vector of (l, v) with l = 0 (x'v <= 1) or l >= 1 (X[l,:]v <= x[l])."""
function separate_new(A::AbstractMatrix, x::Vector{Float64}, X::Matrix{Float64},
                      pool::Vector{Vector{Float64}}; eps = 1e-4, tol = 1e-3, k_per = 3)
    n = length(x)
    Ab = A .!= 0
    W = hcat(x, X)                                       # column 1: x̂ ; column l+1: X̂[:, l]
    rhs = vcat(1.0, x)
    P = isempty(pool) ? zeros(n, 0) : reduce(hcat, pool)
    scores = P' * W                                      # (e) pool scored for every l at once
    cuts = Tuple{Int,Vector{Float64}}[]
    for l in 0:n
        col = l + 1
        (l > 0 && x[l] < tol) && continue                # (a) row l numerically zero
        w = max.(W[:, col], 0.0)
        V = findall(w .> tol)                            # (a) support
        sum(w[V]) - rhs[col] <= eps && continue
        cands = Dict{Vector{Int},Float64}()
        for p in findall(scores[:, col] .- rhs[col] .> eps)
            cands[findall(P[:, p] .> 0.5)] = scores[p, col] - rhs[col]
        end
        for Sloc in heuristic_sets(Ab[V, V], w[V])       # (b)(c)(d)
            Sg = V[Sloc]
            viol = sum(W[Sg, col]) - rhs[col]
            viol > eps && (cands[Sg] = viol)
        end
        for (Sg, _) in first(sort(collect(cands), by = kv -> -kv[2]), k_per)
            v = zeros(n)
            v[Sg] .= 1.0
            push!(cuts, (l, v))
        end
    end
    return cuts
end

# ---------------------------------------------------------------------------------------
# BQPc relaxation
# ---------------------------------------------------------------------------------------

"""Build BQPc with scale factor γ. Returns (model, x, X)."""
function build_bqpc(C::Matrix{Float64}, A::AbstractMatrix, s::Int, γ::Float64;
                    eps = 1e-5, max_iters = 100_000)
    n = size(C, 1)
    model = Model(SCS.Optimizer)
    set_optimizer_attribute(model, "eps_abs", eps)
    set_optimizer_attribute(model, "eps_rel", eps)
    set_optimizer_attribute(model, "max_iters", max_iters)
    set_optimizer_attribute(model, "verbose", 0)

    @variable(model, Y[1:n+1, 1:n+1], PSD)               # Y = [1 x'; x X] ⪰ 0
    x = Y[2:end, 1]
    X = Y[2:end, 2:end]
    @constraint(model, Y[1, 1] == 1)
    @constraint(model, sum(x) == s)
    @constraint(model, [i = 1:n], sum(X[i, j] for j in 1:n) == s * x[i])
    @constraint(model, [i = 1:n], X[i, i] == x[i])
    @constraint(model, sum(X[i, j] for i in 1:n, j in 1:n if A[i, j] != 0) == s * (s - 1))

    # objective: ldet(γ C∘X + Diag(e - x)) via the log-det cone (u = 1)
    M(i, j) = i == j ? γ * C[i, i] * X[i, i] + 1 - x[i] : γ * C[i, j] * X[i, j]
    @variable(model, t)
    @constraint(model, vcat(t, 1.0, [M(i, j) for j in 1:n for i in 1:j]) in
                       MOI.LogDetConeTriangle(n))
    @objective(model, Max, t)
    return model, x, X
end

"""Bound recomputed from the primal point (same as the Python code)."""
function bqp_bound(C, s, γ, xv, Xv)
    M = γ .* C .* Xv .+ Diagonal(1 .- xv)
    return logabsdet(Symmetric(Matrix(M)))[1] - s * log(γ)
end

# ---------------------------------------------------------------------------------------
# New IT2 cutting-plane loop
# ---------------------------------------------------------------------------------------

function it2_new(C::Matrix{Float64}, A::AbstractMatrix, s::Int, γ::Float64;
                 eps = 1e-4, tol = 1e-3, k_per = 3, max_rounds = 200)
    n = size(C, 1)
    model, x, X = build_bqpc(C, A, s, γ)
    seen = Set{Tuple{Int,Vector{Int}}}()
    pool = Vector{Vector{Float64}}()
    poolkeys = Set{Vector{Int}}()
    t_solve = t_sep = 0.0
    rounds, ncuts, bound, final_time = 0, 0, NaN, 0.0
    while true
        t = @elapsed optimize!(model)
        t_solve += t
        final_time = t
        rounds += 1
        xv = value.(x)
        Xv = value.(X)
        Xv = (Xv + Xv') / 2
        bound = bqp_bound(C, s, γ, xv, Xv)
        rounds > max_rounds && break

        new = 0
        t_sep += @elapsed begin
            for (l, v) in separate_new(A, xv, Xv, pool; eps = eps, tol = tol, k_per = k_per)
                I = findall(v .> 0.5)
                key = (l, I)
                key in seen && continue
                push!(seen, key)
                new += 1
                if l == 0
                    @constraint(model, sum(x[i] for i in I) <= 1)
                else
                    @constraint(model, sum(X[l, i] for i in I) <= x[l])
                end
                if !(I in poolkeys)
                    push!(poolkeys, I)
                    push!(pool, v)
                end
            end
        end
        ncuts += new
        new == 0 && break
    end
    return (bound = bound, ncuts = ncuts, rounds = rounds, total_time = t_solve + t_sep,
            sep_time = t_sep, final_time = final_time,
            status = termination_status(model))
end

# ---------------------------------------------------------------------------------------
# Example
# ---------------------------------------------------------------------------------------
# using MAT, Random
# C = Matrix{Float64}(matread("data63.mat")["A"])        # your file stores C under "A"
# n, p, s = 63, 0.8, 8
# Random.seed!(4)                                           # Julia's RNG ≠ numpy's: different graph
# U = triu(rand(n, n) .< p, 1); G = Int.(U .| U')           # adjacency of G(n,p)
# γ = 6.0757                                                # optimal γ for BQPc at this s
# res = it2_new(C, G, s, γ)
# println(res)
