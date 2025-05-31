

%Read the file and construct the adjacency matrix
filename = 'bio-DM-LC.edges'; % load file
data = load(filename);

%% Extract edges and weights
% Extract edges and weights
edges = data(:, 1:2); % one has nodes and 2 has edges
weights = data(:, 3);
nodes = unique(edges(:));
n = max(nodes); % Assuming nodes are labeled from 1 to max node
fprintf('Number of nodes in the graph are: %d\n', n)

% Construct weighted adjacency matrix
g = graph(edges(:, 1), edges(:,2), weights)
g_adj = adjacency(g); 
degree_from_graph = degree(g);
fprintf("The max degree from the graph is: %d", max(degree_from_graph));

A = zeros(n, n);
for i = 1:size(edges, 1)
    A(edges(i,1), edges(i,2)) = weights(i);
    A(edges(i,2), edges(i,1)) = weights(i); % Assuming an undirected graph is not required 
end
%disp(g_adj == A)


% (a) Compute nodes with the highest degree
degree_from_mat = sum(A > 0, 2); % sums along the row, 2 tells that it's row wise
[max_degree, max_nodes] = max(degree_from_mat); % returns indice as well 
fprintf('Node(s) with highest degree: %s with degree %d\n', num2str(max_nodes), max_degree);

% (b) Plot Degree Distribution
unique_degrees = unique(degree_from_mat);
distribution = histcounts(degree_from_mat, [unique_degrees; max(unique_degrees)+1]); %  histcounts(data, bin_edges) 
figure;
stem(unique_degrees, distribution, 'filled');
xlabel('Degree k');
ylabel('N(k)');
title('Degree Distribution');

%% (c) Compute Clustering Coefficient Distribution
clustering_coeffs = zeros(n, 1);
for i = 1:n
    neighbors = find(A(i, :));
    k_i = length(neighbors);
    if k_i > 1
        subgraph = A(neighbors, neighbors) > 0;
        actual_edges = sum(subgraph(:)) / 2;
        clustering_coeffs(i) = (2 * actual_edges) / (k_i * (k_i - 1));
    else
        clustering_coeffs(i) = 0;
    end
end
% clustering coeficient distribution
figure;
histogram(clustering_coeffs, 10);
xlabel('Clustering Coefficient');
ylabel('Frequency');
title('Clustering Coefficient Distribution');

% (d) Compute Characteristic Path Length
G = graph(A);
D = distances(G); % Compute shortest paths
D(D == Inf) = 0; % Ignore disconnected components
characteristic_path_length = sum(D(:)) / (n * (n - 1));
fprintf('Characteristic Path Length: %.4f\n', characteristic_path_length);

%% Alternate part one: with In-built function

%% Create a Regular Lattice 


g_lat= make_ring_lattice(10,4);
figure;
h = plot(g_lat,'Layout', 'force', 'UseGravity', true); 
title('g_$lat$');


function G = make_ring_lattice(n, k)
    halfk = floor(k / 2);
    edges = [];
    
    for u = 1:n
        for j_offset = 1:halfk
            v = mod((u - 1) + j_offset, n) + 1;
            edges = [edges; u v]; % vertically concatenates two vectors u and v 
        end
    end
    
    G = graph(edges(:,1), edges(:,2));
end

