A = [
    0 1 1 0; 
    1 0 1 1; 
    1 1  0 1;
    0 1 1 0];  % adjacency matric of an undirected graph is always symmetirc 



% load the graph from adjacency matrix
G = graph(A);

% Compute the degrees of G
degree_of_g = degree(G) ; 
matrix_of_degree = diag(degree_of_g) ; 


% plot the histogram of the degrees
figure; 
histogram(degree_of_g); 
xlabel('Degree')
ylabel('Frequency')
title('DEgree  distribution of G')


% compute the clustering coeficient
% directly uses the adj matrix as a sparse matrix! 

C = clustering_coefficients(sparse(A)); % uses BGL! 


% calculate the characteristic path length
% Dij represents the shortest path bw noes i and j
distance_matrix = all_shortest_paths(sparse(A)); % MATLAB BGL USES SPARSE MATRIX
distance_matrix(distance_matrix == inf) = NaN; 
% find mean path length 
mean_path_length = mean(distance_matrix(:)) ; 
fprintf("The distance matrix is \n {}") ;

%% Do it for an edge list 

% Example edge list (each row represents an edge)
edge_list = [
    1 2; 
    1 3; 
    2 3; 
    2 4; 
    3 4]; 

G_edge = graph(edge_list(:, 1), edge_list(:, 2)) ; 

% Compute the degrees of G 
degree_of_g_edge= degree(G_edge); 

% Plot degree dist of G-edge
figure; 
histogram(degree_of_g_edge); 
xlabel('degree of g_edge')
ylabel('Frequency')
title('Degree Distribution')

% Calculate the clustering coefiicent
adj_edge_list = adjacency(G_edge);
C_edge = clustering_coefficients(sparse(adj_edge_list));

% Charateristic path lenght 
D_edge = all_shortest_paths(sparse(adj_edge_list)); 
Dedge(D_edge == inf) = NaN; 
L_edge = mean(D_edge) ; 



 