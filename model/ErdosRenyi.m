function [h, adjMatrix] = ErdosRenyi(N,p)
% [h, adjMatrix] = ErdosRenyi(N,p) returns an Erdős–Rényi model random graph
% with N nodes, and p the probability of connection of every possible edge.

% To get average degree k, set p = k/(N-1). Exact target not guaranteed due
% to stochasticity.

adjMatrix = zeros(N,N);
% s = [];
% t = [];
for v = 1:N
    for w = v+1:N
        if rand() < p
            % s(end + 1) = 
            adjMatrix(v,w) = 1;
            adjMatrix(w,v) = 1;
        end
    end
end

h = graph(adjMatrix);
end