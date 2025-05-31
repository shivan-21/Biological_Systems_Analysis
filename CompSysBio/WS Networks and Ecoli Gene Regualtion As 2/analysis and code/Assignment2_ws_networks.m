rng(0); % seed the random value for consistency in analysis

% Generate graphs using the WattsStrogatz function
g_low_p  = WattsStrogatz(100, 5, 0.3);
g_mid_p  = WattsStrogatz(100, 5, 0.7);
g_high_p = WattsStrogatz(100, 5, 1.0);

% Store the graphs and corresponding beta values in cell/array for easy looping
graphs = {g_low_p, g_mid_p, g_high_p};
betas  = [0.3, 0.7, 1.0];

% Preallocate a structure array to store results
results = struct('beta', [], 'graph', [], 'uniqueDegrees', [], 'p_deg', [], 'avgClustering', [], 'charPathLength', []);

% Loop over each graph to calculate properties and plot degree distribution
for i = 1:length(graphs)
    G = graphs{i};
    
    % --- Degree Distribution ---
    deg = degree(G);              % Get degree of each node
    edges = min(deg):max(deg);      % Unique degree values (bins)
    counts = histcounts(deg, [edges, max(edges)+1]); % Count nodes in each bin
    p_deg = counts / numnodes(G);   % Normalize to get probability
    
    % Plot the degree distribution
    figure;
    bar(edges, p_deg, 'FaceColor', [0.2 0.6 0.5]);
    xlabel('Degree');
    ylabel('Probability');
    title(sprintf('Degree Distribution for WS Network BETA = %.1f, (SEED 0)', betas(i)));
    
    % --- Clustering Coefficient ---
    % clustering coefficients for each node.
    clusts = clustering_coefficients(sparse(adjacency(G))); 
    avgClust = mean(clusts); % avg clustering coeficient 
    
    % --- Characteristic Path Length ---
    % Calculate the pairwise shortest paths and compute the average (ignoring zeros and infinities)
    D = distances(G);
    D = D(~isinf(D) & D > 0);       % Remove infinite distances and self-distances (zeros)
    charPathLength = mean(D);
    
    % --- Store results for future use ---
    results(i).beta = betas(i);
    results(i).graph = G;
    results(i).uniqueDegrees = edges;
    results(i).p_deg = p_deg;
    results(i).avgClustering = avgClust;
    results(i).charPathLength = charPathLength;
end


%% Create Random networks for comparision 
% disp(numedges(g_low_p)) | All WaatStrogatts graphs have 500 edges 
rng(0);  rand_g1 = ErdosRenyi(100,500);
rng(1); rand_g2 = ErdosRenyi(100, 500);
rng(2); rand_g3 = ErdosRenyi(100, 500); 

% repeat analogous process to calculate important grpah properties 
ER_graphs = {rand_g1, rand_g2, rand_g3};

% Preallocate a structure array to store results for ER graphs
er_results = struct('seed', [], 'graph', [], 'degreeBins', [], ...
                    'degreeProb', [], 'avgClust', [], 'charPathLen', []);

% Loop over each ER graph to calculate and store properties
for j = 1:length(ER_graphs)
    % Retrieve the current ER graph
    G_er = ER_graphs{j};
    
    % --- Degree Distribution ---
    % Compute the degrees of all nodes in the graph
    nodeDegrees = degree(G_er);
    % Define bins from the minimum to maximum degree found
    degreeBins = min(nodeDegrees):max(nodeDegrees);
    % Count the number of nodes for each degree bin
    counts = histcounts(nodeDegrees, [degreeBins, max(degreeBins)+1]);
    % Normalize counts to get the probability distribution
    degreeProb = counts / numnodes(G_er);
    
    % Plot the degree distribution with a new bar color [0.8 0.3 0.3]
    figure;
    bar(degreeBins, degreeProb, 'FaceColor', [0.8, 0.3, 0.3]);
    xlabel('Degree');
    ylabel('Probability');
    title(sprintf('Degree Distribution for Erdos-Renyi Graph (Seed = %d)', j-1));
    
    % --- Clustering Coefficient ---
    % Compute clustering coefficients for each node.
    % (Assuming you have a function 'clustering_coefficients' available.)
    clusterVals = clustering_coefficients(sparse(adjacency(G_er)));
    avgCluster = mean(clusterVals);  % Compute the average clustering coefficient
    
    % --- Characteristic Path Length ---
    % Calculate pairwise shortest path lengths
    D = distances(G_er);
    % Exclude infinite distances and self-distances (zeros)
    validD = D(~isinf(D) & D > 0);
    charPathLen = mean(validD);  % Compute the average path length
    
    % --- Store the results for future use ---
    er_results(j).seed = j-1;          % Based on the rng seeds used (0, 1, 2)
    er_results(j).graph = G_er;
    er_results(j).degreeBins = degreeBins;
    er_results(j).degreeProb = degreeProb;
    er_results(j).avgClust = avgCluster;
    er_results(j).charPathLen = charPathLen;
end

%% Comparison of Watts–Strogatz (WS) Networks with Erdős–Rényi (ER) Graphs

% First, display a text summary comparing each WS network with each ER graph
fprintf('================ Comparison Report ================\n\n');
for wsIdx = 1:length(results)
    fprintf('WS Network (beta = %.1f):\n', results(wsIdx).beta);
    fprintf('  Avg. Clustering Coefficient: %.4f\n', results(wsIdx).avgClustering);
    fprintf('  Characteristic Path Length : %.4f\n', results(wsIdx).charPathLength);
    fprintf('  Degree Distribution: [min = %d, max = %d]\n', ...
        results(wsIdx).uniqueDegrees(1), results(wsIdx).uniqueDegrees(end));
    fprintf('---------------------------------------------------\n');
    for erIdx = 1:length(er_results)
        fprintf('   ER Graph (seed = %d):\n', er_results(erIdx).seed);
        fprintf('     Avg. Clustering Coefficient: %.4f\n', er_results(erIdx).avgClust);
        fprintf('     Characteristic Path Length : %.4f\n', er_results(erIdx).charPathLen);
        fprintf('     Degree Distribution: [min = %d, max = %d]\n', ...
            er_results(erIdx).degreeBins(1), er_results(erIdx).degreeBins(end));
    end
    fprintf('\n');
end



%% Additionally, plot comparisons of the key parameters

% For plotting, extract the WS beta values and properties.
ws_betas    = [results.beta];
ws_clustering = [results.avgClustering];
ws_path     = [results.charPathLength];

% For the ER graphs, we will use all three values (they were computed with different seeds)
er_clustering = [er_results.avgClust];
er_path       = [er_results.charPathLen];

% Plot Comparison of Average Clustering Coefficient
figure;
subplot(1,2,1);
% Plot WS networks as a function of beta
bar(ws_betas, ws_clustering, 0.5, 'FaceColor', [0.2, 0.6, 0.5]);
hold on;
% Also plot each ER graph's value as horizontal lines for clarity
for erIdx = 1:length(er_results)
    yline(er_results(erIdx).avgClust, '--', sprintf('ER (seed=%d)', er_results(erIdx).seed),...
        'LineWidth',1.5, 'LabelHorizontalAlignment','left');
end
xlabel('WS Rewiring Probability (beta)');
ylabel('Avg. Clustering Coefficient');
title('Clustering Coefficient Comparison');
hold off;

% Plot Comparison of Characteristic Path Length
subplot(1,2,2);
bar(ws_betas, ws_path, 0.5, 'FaceColor', [0.3, 0.4, 0.9]);
hold on;
for erIdx = 1:length(er_results)
    yline(er_results(erIdx).charPathLen, '--', sprintf('ER (seed=%d)', er_results(erIdx).seed),...
        'LineWidth',1.5, 'LabelHorizontalAlignment','left');
end
xlabel('WS Rewiring Probability (beta)');
ylabel('Characteristic Path Length');
title('Path Length Comparison');
hold off;

%% Interpretation (to be done by the researcher)
% You can now compare:
% - How the clustering coefficient decreases as beta increases (especially when beta = 1)
% - How the characteristic path length changes as beta increases.
% - How these properties of the WS networks compare with the properties of the ER graphs.
%
% For example, as beta tends to 1, the WS networks should approach the properties
% of a random graph (i.e. similar clustering and path length to the ER graphs).

%% Interpretation (to be done by the researcher)
% You can now compare:
% - How the clustering coefficient decreases as beta increases (especially when beta = 1)
% - How the characteristic path length changes as beta increases.
% - How these properties of the WS networks compare with the properties of the ER graphs.
%
% For example, as beta tends to 1, the WS networks should approach the properties
% of a random graph (i.e. similar clustering and path length to the ER graphs).


%% Function to generate a Watts Strogatz Network 
function h = WattsStrogatz(N, K, beta)
    % H = WattsStrogatz(N,K,beta) returns a Watts-Strogatz model graph with N
    % nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
    % ADDS SYMMETRICALLY ON EITHER SIDE. therefore, for 10 nearest
    % neighbours put K =5! 
    % beta = 0 is a ring lattice, and beta = 1 is a random graph.
    
    % Provide default values if not enough inputs are given
    if nargin < 3
        beta = 0.3;  % default rewiring probability
    end
    if nargin < 2
        K = 5;       % default: mean degree will be 2*K = 10 neighbors total
    end
    if nargin < 1
        N = 100;     % default number of nodes
    end
    
    % Connect each node to its K next neighbors. This constructs indices for a ring lattice.
    s = repelem((1:N)', 1, K); % repelem((1:N)', 1, K) repeats each node index K times horizontally, resulting in an N×K edge list
    t = s + repmat(1:K, N, 1); % repmat(1:K, N, 1) creates an N×K matrix where each row is the vector [1, 2, ..., K]
    t = mod(t-1, N) + 1;% modulo operation ensures that the node indices wrap around.

    
    % Rewire the target node of each edge with probability beta
    for source = 1:N    
        switchEdge = rand(K, 1) < beta; % rand(K, 1) generates a K×1 vector of random numbers between 0 and 1
        
        newTargets = rand(N, 1); % A random number is generated for each node (creating an N×1 vector) as a score
        newTargets(source) = 0; % no self-connections 
        newTargets(s(t==source)) = 0; % no double connections 
        newTargets(t(source, ~switchEdge)) = 0; % allows only those nodes that can rewire
        
        [~, ind] = sort(newTargets, 'descend'); % sort candidate nodes in newTargets in descending order
        t(source, switchEdge) = ind(1:nnz(switchEdge)); 
    end

    % Create the graph. (s and t are matrices but will be linearized automatically.)
    h = graph(s, t);

end

%% Comparing each WattStrogats graph to all thress random networks.

%% Function to Create Erdos Renyi random Network | n nodes and m edges 
function G = ErdosRenyi(n, m)


%   The function randomly selects m edges from all possible n*(n-1)/2
%   edges between n nodes.

    % Check that m is not greater than the maximum number of possible edges
    maxEdges = n*(n-1)/2;
    if m > maxEdges
        error('The number of edges m exceeds the maximum possible (%d) for %d nodes.', maxEdges, n);
    end

    % Generate a list of all possible edges (each row is a pair [i, j] with i < j)
    allEdges = nchoosek(1:n, 2);
    
    % Randomly select m edges from all possible edges
    idx = randperm(size(allEdges, 1), m);
    selectedEdges = allEdges(idx, :);
    
    % Create the graph using the selected edges
    G = graph(selectedEdges(:,1), selectedEdges(:,2));
end
