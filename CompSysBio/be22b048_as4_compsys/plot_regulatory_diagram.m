function plotGeneRegulatoryNetwork()
% plotGeneRegulatoryNetwork  Draws the MDM2/MYC–P53/RB–Apoptosis/Proliferation network
% 
%   Activation edges (traditionally, pointed arrows) are colored dark green.
%   Inhibition edges (traditionally, flat bars) are colored red.
%
% Nodes:       MDM2, MYC, P53, RB, Apoptosis, Proliferation
% Activations: MDM2→MDM2,  P53→Apoptosis,  MYC→P53,  MYC→Proliferation
% Inhibitions: MDM2⊣P53,   MDM2⊣Apoptosis, MYC⊣RB,  RB⊣Proliferation

    % Define nodes
    nodeNames = {'MDM2','MYC','P53','RB','Apoptosis','Proliferation'};

    % List of activation edges (green arrows)
    actSrc = { 'MDM2',      'P53',          'MYC',         'MYC'      };
    actTgt = { 'MDM2',      'Apoptosis',    'P53',         'Proliferation' };

    % List of inhibition edges (arrows arrows)
    inhSrc = { 'MDM2',      'MDM2',         'MYC',         'RB', 'RB'     };
    inhTgt = { 'P53',       'Apoptosis',    'RB',          'Proliferation', 'MYC'};

    % Build directed graph
    s = [actSrc, inhSrc];
    t = [actTgt, inhTgt];
    G = digraph(s, t, [], nodeNames);

    % Plot with a layered layout
    h = plot(G, ...
             'Layout','layered', ...
             'NodeFontSize',12, ...
             'ArrowSize', 15, ...
             'LineWidth', 1.5);

    % Color the activation edges green
    highlight(h, actSrc, actTgt, ...
              'EdgeColor', [0 0.5 0]);

    % Color the inhibition edges red
    highlight(h, inhSrc, inhTgt, ...
              'EdgeColor', [1 0 0]);

    % tweak axes
    axis off
    title('Gene Regulatory Network','FontWeight','normal')
    
end
