function [h, adjMatrix] = HolmeKim(N,m,mt,N0)
% [h, adjMatrix] = HolmeKim(N,m,mt,N0) returns a Holme-Kim model scale-free clustered graph
% with N nodes, m edges added at each time step, average mt edges added 
% via triad-formation steps during a time step, and N0 initial vertices.

% Based on Petter Holme, Beom Jun Kim (2001), Growing scale-free networks 
% with tunable clustering. Physical Review E, Volume 65, 026107

% mt = 0 yields the Barabási–Albert model scale-free graph. 
% Clustering coefficient increases linearly with mt
% To get average degree k, set m = k/2, N0 = k + 1

assert(mt <= m - 1, 'Average number of triad formation edges can be atmost one less than number of edges being added in a time step');
Pt = mt / (m - 1); 	% probability of choosing TF step
adjMatrix = zeros(N,N);

% Start with fully connected subgraph on first N0 nodes
for v = 1:N0
    for w = 1:N0
        if v ~= w
            adjMatrix(v,w) = 1;
        end
    end
end

% Add remaining nodes via PA and TF steps
for v = N0 + 1:N	% each time step a new vertex is connected to graph
	wPAPrev = 0;	% resetting for safety
	for i = 1:m	% running through PA + TF steps
		x = rand(2,1);
		degVector = sum(adjMatrix, 2);
		w = 0;

		if i > 1 && x(1) <= Pt	% do TF step
			validNeighbours = find(adjMatrix(wPAPrev,:) ~= 0 & adjMatrix(v,:) == 0 & (1:N) ~= v);    % nodes that are neighbours of node taht was attached to in recent-most PA step, while being non-neighbours of current node v and not v itself
			if ~isempty(validNeighbours) % if no valid neighbour, do PA step
				w = validNeighbours(ceil(x(2) * numel(validNeighbours)));
			end
		end
		if w == 0 % do PA step if TF was not tried or failed
			% w = find(cumsum(degVector + 1) / sum(degVector + 1) >= x(2), 1);
            nonNeighbours = find(adjMatrix(v,:) == 0 & (1:N) ~= v);   % nodes that are not connected to v, and excluding v
            degVector_NN = degVector(nonNeighbours);
            w = nonNeighbours(find(cumsum(degVector_NN) / sum(degVector_NN) >= x(2), 1));
			wPAPrev = w;
		end
		
        assert(v~=w);
		adjMatrix(v,w) = 1;
		adjMatrix(w,v) = 1;
	end
end

h = graph(adjMatrix);
end