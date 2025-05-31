n = 100; % Number of nodes
k = 4;   % Each node is connected to k nearest neighbors 

% Initialize adjacency matrix
A = zeros(n);

% Create edges for k/2 nearest neighbors on both sides
for shift = 1:k/2
    % Connect each node to its shift-th neighbor (only works for non-boundary indices)
    A = A + diag(ones(n-shift, 1), shift) + diag(ones(n-shift, 1), -shift);
end

% Add wrap-around connections (A(1,99) & A(2,100)) and (A(99,1) and A(100,2)
% This adds ones at positions:
%   - diag(ones(k/2, 1), n-k/2) places ones in A(1, n-k/2+1) ... A(k/2, n)
%   - diag(ones(k/2, 1), -(n-k/2)) places ones in A(n-k/2+1, 1) ... A(n, k/2)
A = A + diag(ones(k/2, 1), n-k/2) + diag(ones(k/2, 1), -(n-k/2));

% Manually add the missing wrap-around edge between node 1 and node n
A(1, n) = 1;
A(n, 1) = 1;

% Convert to a graph object
G = graph(A);


% Plot the graph
figure;
plot(G, 'Layout', 'circle');
title('Regular Lattice Graph (n ==100, k== 4)');


%% operations on the graph 

%{
threshold =4; 
if all(degree_of_g == threshold)
    disp('All values are below the threshold.')
else
    disp('weird behaviour')
end
%}

% Display the Degree Distribution 
node_degrees = sum(A, 2); % col sum and row sum yield the same thing. 



disp("Degree distribution:");
disp(tabulate(node_degrees)); % Should show only degree k = 4


% plot histogram 
figure; 
histogram(degree_of_g); 
xlabel('Degree') ; ylabel('frequency'); title('Degree Distribution of Lattice Graph'); 
