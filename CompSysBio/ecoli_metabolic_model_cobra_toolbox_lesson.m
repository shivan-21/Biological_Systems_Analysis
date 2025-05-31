clear all; 
% Initialise the model 
model = readCbModel('e_coli_core.xml');
% model parameters
% rxns: A cell array of reaction identifiers
% mets: A cell array of metabolite names.
%S: The stoichiometric matrix 
%lb and ub: Vectors containing the lower and upper bounds for each reaction.
%c: A vector where nonzero entries indicate the reaction(s) used as the objective (e.g., biomass reaction).


model.rxns(find(model.c));
BiomassRxn = checkObjective(model);

% find the reaction IDs
BiomassID = findRxnIDs(model, BiomassRxn);
printUptakeBound(model); 

%% Change the reachtion bounds and solve
model1 = changeRxnBounds(model, 'EX_o2_e', 0, 'l');
printUptakeBound(model1)
solution1 = optimizeCbModel(model1);


%% Change Reaction Bounds and run 
model2 = changeRxnBounds(model, {'EX_glc__D_e', 'EX_nh4_e'}, [-5, -20], 'l');
printUptakeBound(model2);
solution2 = optimizeCbModel(model2)


%% Change the objective function 
model3 = changeObjective(model, 'EX_ac_e');
bio_rxn = model.rxns(find(model.c));

model3 = changeRxnBounds(model3, bio_rxn, 0.5*solution1.f, 'l');
solution3 = optimizeCbModel(model3)


%% Delete a single gene and run the model 

model4 = readCbModel('e_coli_core.xml');

[grRatio, grKO, grWT, hasEffect] = singleGeneDeletion(model4);
EssGenes = model4.genes(grRatio < 1e-3);
% Reaction Deletion 
[grRatio_r, gr_KO_r, grWT_r, hasEffect_r] = singleRxnDeletion(model4); 
% Essential Reactions 
EssRxns = model4.rxns(grRatio_r < 1e-3); 
solution4 = optimizeCbModel(model4)

%% Flux variability analysis 
model5 = readCbModel('e_coli_core.xml') ; 
EssRxns = model5.rxns(grRatio_r < 1e-3); 

% Flux variability 
[min, max] = fluxVariability(model5)

% flux vatiabiity for various reactions 
[min, max] = fluxVariability(model, 'rxnNameList', {'PGK', 'PGM'});