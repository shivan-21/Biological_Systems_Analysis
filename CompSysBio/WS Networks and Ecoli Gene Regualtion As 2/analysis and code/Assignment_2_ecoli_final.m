%% Seting up the Graph from provided table 

% Read the tab-separated file
T = readtable('E_coli_transcription_network.txt','Delimiter','\t');

% Create a directed graph: source is transcription factor, target is operon
G = digraph(T.transcription_factor, T.operon);

% Store the regulation type as an edge attribute
G.Edges.Regulation = T.regulation_type;

% Define a color mapping for each regulation type
colorMap = containers.Map({'activator','repressor','dual'},...
    {[0, 0.8, 0], [1, 0, 0], [0, 0, 1]});

%%  Visulise the network using a force-directed layout
figure;
h = plot(G,'Layout', 'force', 'UseGravity', true); 
title('E. coli Transcription Network');

% Set color based on the regulation type of each edge
for i = 1:numedges(G)
    regType = lower(G.Edges.Regulation{i});
    if isKey(colorMap, regType)
        edgeColor = colorMap(regType);
    else
        edgeColor = [0, 0, 0];  % default to black if unrecognized
    end
    highlight(h, 'Edges', i, 'EdgeColor', edgeColor, 'LineWidth', 1.5);
end

% Add manual legend for edge colors
hold on;
hActivator = plot(nan, nan, 'Color', [0, 0.8, 0], 'LineWidth', 2); % activator (green)
hRepressor = plot(nan, nan, 'Color', [1, 0, 0], 'LineWidth', 2);   % repressor (red)
hDual     = plot(nan, nan, 'Color', [0, 0, 1], 'LineWidth', 2);    % dual (blue)
legend([hActivator, hRepressor, hDual], {'Activator', 'Repressor', 'Dual'}, 'Location', 'Best');
hold off;

%%  Calculate Degree Centrality, Closeness Centrality and Shortest Path betweenness centrality | Report top 5 

% Compute centrality measures for each node in graph G
degreeCentrality = indegree(G) + outdegree(G); % sum of indegree and outdegree
closenessCentrality = centrality(G, 'outcloseness'); % can be changed to incloseness or avg later
betweennessCentrality = centrality(G, 'betweenness'); % inbuilt func for node bw shortest path 

% Create a table with node names and their centrality values
centralityTable = table(G.Nodes.Name, degreeCentrality, closenessCentrality, betweennessCentrality, ...
    'VariableNames', {'Node', 'Degree', 'Closeness', 'Betweenness'});

% Identify transcription factors from the original data table T
transcriptionFactors = unique(T.transcription_factor);

% Filter the centrality table to include only transcription factors
isTF = ismember(centralityTable.Node, transcriptionFactors);
tfCentralityTable = centralityTable(isTF, :);

% --- Top Five by Degree Centrality ---
[~, idxDegree] = sort(tfCentralityTable.Degree, 'descend');
top5Degree = tfCentralityTable(idxDegree(1:min(5, end)), :);
disp('Top 5 Transcription Factors by Degree Centrality:');
disp(top5Degree);

% --- Top Five by Closeness Centrality ---
[~, idxCloseness] = sort(tfCentralityTable.Closeness, 'descend');
top5Closeness = tfCentralityTable(idxCloseness(1:min(5, end)), :);
disp('Top 5 Transcription Factors by Closeness Centrality:');
disp(top5Closeness);

% --- Top Five by Betweenness Centrality ---
[~, idxBetweenness] = sort(tfCentralityTable.Betweenness, 'descend');
top5Betweenness = tfCentralityTable(idxBetweenness(1:min(5, end)), :);
disp('Top 5 Transcription Factors by Betweenness Centrality:');
disp(top5Betweenness);

%--- Degree Distribution of Original Network --- %
totalDegree = indegree(G) + outdegree(G);
figure;
histogram(totalDegree);
title('Degree Distribution for Complete Network');
xlabel('Degree');
ylabel('Number of Nodes');

%% Recomputing properties of activator only network after removing top 5 activator TFs
% (Revised using supervisor's approach)

% Compute total number of activating edges in the entire network
total_network_activation = sum(strcmp(T.regulation_type, 'activator'));

% Initialize an array for activation fraction for each TF (for all nodes in G)
activationFraction = zeros(size(G.Nodes, 1), 1);

for i = 1:size(G.Nodes, 1)
    tf = G.Nodes.Name(i);
    tf_activation = sum(strcmp(T.transcription_factor, tf) & strcmp(T.regulation_type, 'activator'));
    activationFraction(i) = tf_activation / total_network_activation; 
end

% Create a table pairing each TF with its activation fraction and sort it
tfActivationTable = table(G.Nodes.Name, activationFraction, 'VariableNames', {'TF', 'ActivationFraction'});
tfActivationTable = sortrows(tfActivationTable, 'ActivationFraction', 'descend');

% Select the top five transcription factors
top5ActivatorTFs = tfActivationTable.TF(1:min(5, height(tfActivationTable)));

disp('Top 5 Transcription Factors with Highest Fraction of Activation:');
disp(tfActivationTable(1:min(5, height(tfActivationTable)), :));
disp('These activators will be removed from the network. Only activating edges will remain.');

% Remove the Top Five Activator TFs from the Original Network
G_removed = rmnode(G, top5ActivatorTFs);

% Create a New Graph Using Only 'Activator' Edges from the reduced network
edgeFilter = strcmpi(G_removed.Edges.Regulation, 'activator');
filteredSource = G_removed.Edges.EndNodes(edgeFilter, 1);
filteredTarget = G_removed.Edges.EndNodes(edgeFilter, 2);
G_actOnly = digraph(filteredSource, filteredTarget);

% Recompute Network Properties on the Activator-Only Network

% (i) Degree Distribution
totalDegree_actOnly = indegree(G_actOnly) + outdegree(G_actOnly);
figure;
histogram(totalDegree_actOnly);
title('Degree Distribution (Activator-Only Network)');
xlabel('Degree');
ylabel('Number of Nodes');

% (ii) Degree Centrality: sum of indegree and outdegree
degreeCentrality_actOnly = indegree(G_actOnly) + outdegree(G_actOnly);
degreeCentralityTable_actOnly = table(G_actOnly.Nodes.Name, degreeCentrality_actOnly, ...
    'VariableNames', {'Node', 'DegreeCentrality'});
disp('Degree Centrality in Activator-Only Network:');
disp(degreeCentralityTable_actOnly);

% (iii) Closeness Centrality: using the 'outcloseness' measure
closenessCentrality_actOnly = centrality(G_actOnly, 'outcloseness');
closenessCentralityTable_actOnly = table(G_actOnly.Nodes.Name, closenessCentrality_actOnly, ...
    'VariableNames', {'Node', 'ClosenessCentrality'});
disp('Closeness Centrality in Activator-Only Network:');
disp(closenessCentralityTable_actOnly);

%% Recomputing properties of repressor/dual network after removing top 5 repressors
% (Revised using supervisor's approach)

% Compute total number of repressing edges in the entire network
total_network_repression = sum(strcmp(T.regulation_type, 'repressor'));

% Initialize an array for repression fraction for each TF
repressionFraction = zeros(size(G.Nodes, 1), 1);

for i = 1:size(G.Nodes, 1)
    tf = G.Nodes.Name(i);
    tf_repression = sum(strcmp(T.transcription_factor, tf) & strcmp(T.regulation_type, 'repressor'));
    repressionFraction(i) = tf_repression / total_network_repression; 
end

% Create a table pairing each TF with its repression fraction and sort it
tfRepressionTable = table(G.Nodes.Name, repressionFraction, 'VariableNames', {'TF', 'RepressionFraction'});
tfRepressionTable = sortrows(tfRepressionTable, 'RepressionFraction', 'descend');

% Select the top five transcription factors
top5RepressorTFs = tfRepressionTable.TF(1:min(5, height(tfRepressionTable)));

disp('Top 5 Transcription Factors with Highest Fraction of Repression:');
disp(tfRepressionTable(1:min(5, height(tfRepressionTable)), :));
disp('These repressors will be removed from the network. Only repressing and dual edges will remain.');

% Remove the Top Five Repressor TFs from the Original Network
G_removed = rmnode(G, top5RepressorTFs);

% Create a new network with only 'repressor' and 'dual' edges
edgeFilter = strcmpi(G_removed.Edges.Regulation, 'repressor') | strcmpi(G_removed.Edges.Regulation, 'dual');
filteredSource = G_removed.Edges.EndNodes(edgeFilter, 1);
filteredTarget = G_removed.Edges.EndNodes(edgeFilter, 2);
G_rep = digraph(filteredSource, filteredTarget);

% Recompute network properties on the Repressor/Dual Network

% (i) Degree Distribution
totalDegree_rep = indegree(G_rep) + outdegree(G_rep);
figure;
histogram(totalDegree_rep);
title('Degree Distribution (Repressor/Dual Network)');
xlabel('Degree');
ylabel('Number of Nodes');

% (ii) Degree Centrality: Sum of indegree and outdegree for each node
degreeCentrality_rep = indegree(G_rep) + outdegree(G_rep);
degreeCentralityTable_rep = table(G_rep.Nodes.Name, degreeCentrality_rep, ...
    'VariableNames', {'Node', 'DegreeCentrality'});
disp('Degree Centrality in Repressor/Dual Network:');
disp(degreeCentralityTable_rep);

% (iii) Closeness Centrality: Using 'outcloseness' measure
closenessCentrality_rep = centrality(G_rep, 'outcloseness');
closenessCentralityTable_rep = table(G_rep.Nodes.Name, closenessCentrality_rep, ...
    'VariableNames', {'Node', 'ClosenessCentrality'});
disp('Closeness Centrality in Repressor/Dual Network:');
disp(closenessCentralityTable_rep);
