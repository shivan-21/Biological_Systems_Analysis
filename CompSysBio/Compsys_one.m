%% create Graphs using list-maps and an adjacency matrix 

figure; 

G0 = graph(); 
G0 = addnode(G0, {'1', '2', '3', '4', '5'}) ; 
G0 = addedge(G0, [1 2 3 3], [3,5,4,5]); % adds edges like 1->3, 2->5, 3->4 etc. i.e. pairwise mapping. 
G0.Edges; 
subplot(1,2,1);
plot(G0) ; 
title('Graph G0')
disp(G0)

A = [0,0,1,0,1; 0,0,1,1,1; 1,1,0,1,0; 0,1,1,0,1; 1,1,0,1,0];
G1 = graph(A); 
n_vertices = num_vertices(A); % this is an inbuilt BGL function otherwise | size(A, 1);
n_edge = num_edges(A); % also inbuilt | sum(A(:)) / 2
subplot(1,2,2);
plot(G1);
title('Graph G1')
fprintf('The number of vertices and edges are respectively %d, %d\n', n_vertices, n_edge);

%% 

data = readtable('edge_list.csv'); % read data into a table

nodes = unique([data.Var1, data.Var2]); % unique() removes duplicates.


n = length(nodes) ; % Stores the number of unique nodes

% ismember(A, B) finds the index of each element of A in B
%{
ismember(A, B) checks if elements of A exist in B and returns:
tf: Logical array (1 if found, 0 if not).
target/source: Indices of data.Var1 and data.Var2 in nodes 
%}
[tf, target] = ismember(data.Var1, nodes); 
[tf, source] = ismember(data.Var2, nodes); 


adj = sparse(source, target, data.Var3);  % constructs a sparse adjacency matrix
% row, col, weight 
figure; 
spy(adj);  % plots the sparsity patterns of adjacency matrix
title('Sparsity Patter of Adj Matrix')


figure; 
G2 = graph(source, target, data.Var3); 
G3 = graph(adj); 
% Both are the same grpah constructed differently 
subplot(1,2,1);
plot(G2); 
title('Graph G2')
subplot(1,2,2); 
plot(G3)
title('Graph G3')


%%
g5 = load("Graph1.mat")
a = adjacency(g5.G); 
deg = degree(g5.G); 

figure; 
histogram(deg); 
xlabel('Degree'); 
ylabel('Frequency'); 
title('Degree Distribution of loaded file Graph1.mat') ; 

betweenness = centrality(g5.G, 'betweenness'); 
closeness= centrality(g5.G, 'closeness') ; 