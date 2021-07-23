% Modified from Mathworks Function 2015

function h = WattsStrogatz(N,K,beta)
% H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
% nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
%
% beta = 0 is a ring lattice, and beta = 1 is a random graph.

% Connect each node to its K next and previous neighbors. This constructs
% indices for a ring lattice.
s = repelem((1:N)',1,K);
t = s + repmat(1:K,N,1);
t = mod(t-1,N)+1;

% Rewire the target node of each edge with probability beta
for source=1:N    
    switchEdge = rand(K, 1) < beta;
    
    newTargets = rand(N, 1);
    newTargets(source) = 0;
    newTargets(s(t==source)) = 0;
    newTargets(t(source, ~switchEdge)) = 0;
    
    [~, ind] = sort(newTargets, 'descend');
    t(source, switchEdge) = ind(1:nnz(switchEdge));
end

h1 = graph(s,t);

adj = adjacency(h1);
Hidden_Ind = zeros(N);
for ii = 1:N
    for jj = 1:N
        if adj(ii,jj) == 1
            Hidden_Ind(ii,jj) = 1;
        end
    end
end

h = zeros(N);
I = randperm(N);
for ii = 1:N
    for jj = 1:N
        if Hidden_Ind(I(ii), I(jj)) == 1
            h(ii,jj) = 1;
        end
    end
end
h = graph(Hidden_Ind);

end