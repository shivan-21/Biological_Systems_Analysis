
g_baby = ErdosRenyi(10, 20); 

disp(numedges(g_baby))

function G = ErdosRenyi(n, m)
%ErdosRenyi Generates an Erdős-Rényi random network with n nodes and m edges.
%
%   G = ErdosRenyi(n, m) returns a MATLAB graph object representing an
%   undirected random network with n nodes and exactly m edges.
%
%   The function randomly selects m edges from all possible n*(n-1)/2
%   edges between n nodes.
%
%   Example:
%       G = ErdosRenyi(100, 250); % 100 nodes and 250 edges

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