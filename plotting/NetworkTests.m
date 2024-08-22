   
    % Holme-Kim
N = 1000; k = 12; m = k/2; mt = m-1; N0 = k+1; G = HolmeKim(N, m, mt, N0); plot_subtitle = sprintf('N = %d, m = %d, m_t = %d, N_0 = %d', N,m,mt,N0);
graph_plot(G, N, plot_subtitle);

    %Erdos-Renyi
N = 1000; k = 12; p = k/(N-1); G = ErdosRenyi(N, p); plot_subtitle = sprintf('N = %d, p = %0.2f', N,p);
graph_plot(G, N, plot_subtitle);

    %% Functions
function graph_plot(G, N, plot_subtitle)
    defaultNodeColor = [1 0 0];

    GM = graph_metrics(G);
    degreeList = GM.degreeList;
    GCC = GM.ClusteringCoefficient_Global % Global Clustering Coefficient
    CPL = GM.ShortestPath_Average % Characteristic Path Length
    AD = GM.Degree_Average % average degree
    
    figure(), clf
        NodeColors = 1 - (0 + 1 * degreeList / max(degreeList)) * (1 - defaultNodeColor);
        NodeMarkers = repmat("*", N, 1);
        plot(G, 'NodeColor', NodeColors, 'MarkerSize', 7, 'LineWidth', 1.5, 'Marker', NodeMarkers, 'EdgeColor', 'Black');
        title(['Graph Network: ' plot_subtitle]);    
        subtitle(sprintf('Color intensity indicates degree of node. \\gamma = %.2f, l = %.2f, <k> =  %.2f',GCC, CPL, AD));
            
    figure(), clf
        histogram(degreeList, 'Normalization', 'probability')
        xscale('log'); xlim([1 N]);
        yscale('log'); ylim([1/N 1]);
        title(['Graph Network: ' plot_subtitle]);
        subtitle(sprintf('\\gamma = %.2f, l = %.2f, <k> =  %.2f',GCC, CPL, AD));
        xlabel('Degree, k')
        ylabel('Fraction of Nodes, P(k)')
end
    
function GM = graph_metrics(G)
%%% Calculating Graph Metrics

GM = struct();
GM.degreeList = degree(G);
N = numnodes(G);

GM.ClusteringCoefficient_Local = zeros(N,1);   % local clustering coefficient
num_triplets = 0;
num_triplets_closed = 0;
for node = 1:N
    neighborhood = subgraph(G, neighbors(G, node));
    node_degree = GM.degreeList(node);
    
    num_triplets_local = (node_degree * (node_degree - 1) / 2);
    num_triplets_closed_local = numedges(neighborhood);
    
    GM.ClusteringCoefficient_Local(node) = num_triplets_closed_local / num_triplets_local ;
    num_triplets = num_triplets + num_triplets_local;
    num_triplets_closed = num_triplets_closed + num_triplets_closed_local;
end

GM.ClusteringCoefficient_Local_Avg = mean(GM.ClusteringCoefficient_Local);  % average local clustering coefficient
GM.ClusteringCoefficient_Global = num_triplets_closed / num_triplets; % global clustering coefficient
GM.ShortestPath_Average = sum(distances(G),'all') / (2 * N*(N-1)/2);   % (i,i) has 0 length. Dividing by two to account for double-counting via (i,j) and (j,i)
GM.Degree_Average = mean(GM.degreeList);

end
