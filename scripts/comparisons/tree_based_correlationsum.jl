using DynamicalSystems, Neighborhood


function correlationsum_tree(X, ε; q = 2, metric = Euclidean(), w = 0)
    tree = KDTree(X, metric)
    qf = q == 2 ? identity : x -> float(x)^(q-1)
    _correlationsum_new(tree, X, ε, w, qf)
end

function _correlationsum_tree(tree, X, ε, w, qf)
    theiler = Neighborhood.Theiler(w)
    searchtype = Neighborhood.WithinRange(maximum(ε))
    N = length(X)
    Cq = similar(ε, Float64)
    @inbounds for i in 1:N
        idxs = Neighborhood.isearch(tree, X[i], searchtype, theiler(i))
        isempty(idxs) && continue
        dist = [Neighborhood.evaluate(tree.metric, X[i], X[j]) for j in idxs]
        sort!(dist)
        for k in 1:length(ε)-1
            j = searchsortedfirst(dist, ε[k])
            Cq[k] += qf(j-1)
        end
        Cq[end] += qf(length(dist))
    end
    return Cq ./ (N*qf(N-w-1))
end
