function [model,condGF,targetGF,relaxGF] = runGapfillingFunctions(model,objectiveFunction,biomassReaction,osenseStr,database,curateDefMedia)
% This function runs a set of gapfillling functions on a reconstruction to
% refined as part of the DEMETER pipeline. Reactions are filled in to
% enable flux through an objective function, e.g., the biomass objective
% function.
%
% USAGE:
%
%   [model,condGF,targetGF,relaxGF] = runGapfillingFunctions(model,objectiveFunction,biomassReaction,osenseStr,database)
%
% INPUTS
% model:               COBRA model structure
% objectiveFunction    Reaction ID of the objective function for which flux
%                      will be enabled through gapfilling
% biomassReaction:     Reaction ID of the biomass objective function that
%                      will be restored afterwards (if objectiveFunction is 
%                      not biomassReaction)
% osenseStr:           Maximize ('max')/minimize ('min') linear part of the 
%                      objective
% database:            rBioNet reaction database containing min. 3 columns:
%                      Column 1: reaction abbreviation, Column 2: reaction
%                      name, Column 3: reaction formula.
% OPTIONAL INPUT
% curateDefMedia:      boolean indicating that growth on defined medium is
%                      being curated (default=false)
%
% OUTPUT
% model:               COBRA model structure
% condGF               Reactions added based on conditions (recognizing
%                      certain patterns of reactions)
% targetGF             Reactions added based on tagrted gapfilling 
%                      (specific metabolites that could not be produced)
% relaxGF              Reactions added based on relaxFBA (lowest level of
%                      confidence)
%
% .. Author:
%       - Almut Heinken, 2016-2020

tol=0.0000001;

model = changeObjective(model, objectiveFunction);
modelOld=model;

if nargin<6
    curateDefMedia=0;
end

% Perform gapfilling to enable growth
model = conditionSpecificGapFilling(model, database);

% If model is still unable to grow-test exchanges based on targeted
% metabolite analysis
FBA = optimizeCbModel(model,osenseStr);
if abs(FBA.f) < tol
model = targetedGapFilling(model,osenseStr,database);
end

FBA = optimizeCbModel(model,osenseStr);
if abs(FBA.f) < tol
    if curateDefMedia==1
        % special case: curate defined medium growth
        model = untargetedGapFilling(model,osenseStr,database,1,1,1);
    else
        % try gapfilling based on relaxFBA
        model = untargetedGapFilling(model,osenseStr,database,1,1);
        
        FBA = optimizeCbModel(model,osenseStr);
        if abs(FBA.f) < tol
            % try gapfilling without excluding sink and demand reactions
            model = untargetedGapFilling(model,osenseStr,database,0,0);
        end
    end
end

% test if gapfilled reactions are really needed
[model,condGF,targetGF,relaxGF] = verifyGapfilledReactions(model,osenseStr);

% if the changes make no difference, reverse the changes
FBA = optimizeCbModel(model,osenseStr);
if abs(FBA.f) < tol || FBA.stat==0
    model=modelOld;
    condGF = {};
    targetGF = {};
    relaxGF = {};
end

% change back to biomass reaction
model=changeObjective(model,biomassReaction);

% relax constraints-cause infeasibility problems
relaxConstraints=model.rxns(find(model.lb>0));
model=changeRxnBounds(model,relaxConstraints,0,'l');

% change back to unlimited medium
% list exchange reactions
exchanges = model.rxns(strncmp('EX_', model.rxns, 3));
% open all exchanges
model = changeRxnBounds(model, exchanges, -1000, 'l');
model = changeRxnBounds(model, exchanges, 1000, 'u');

end