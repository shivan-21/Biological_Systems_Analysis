function [model, transportersWithoutExchanges] = findTransportersWithoutExchanges(model)
% Part of the DEMETER pipeline. Finds transporters to extracellular space 
% that are blocked because they have no exchange reaction associated with 
% them.
%
% USAGE   
%       [model, transportersWithoutExchanges] = findTransportersWithoutExchanges(model)
%
% INPUT
% model                             COBRA model structure
%
% OUTPUT
% model                             COBRA model structure
% transportersWithoutExchanges      Removed transport reactions
%
% .. Authors:
% Almut Heinken and Stefania Magnusdottir, 2016-2019

% find extracellular metabolites
ExMets=model.mets(find(contains(model.mets, '[e]')));

% check if an exchange reaction exists
cnt=1;
transportersWithoutExchanges={};
for i = 1:length(ExMets)
    [rxnList, rxnFormulaList] = findRxnsFromMets(model,ExMets{i});
    if ~any(strncmp(rxnList,'EX_',3))
        for j=1:length(rxnList)
            transportersWithoutExchanges{cnt,1}=rxnList{j};
        end
    end
end

% remove transporters without exchanges
transportersWithoutExchanges = unique(transportersWithoutExchanges);
if ~isempty(transportersWithoutExchanges)
% keep extracellular reactions that are not transporters
subs=model.subSystems(find(strcmp(model.rxns,transportersWithoutExchanges)));
transportersWithoutExchanges(~strcmp(subs,'Transport,extracellular'))=[];
model = removeRxns(model, transportersWithoutExchanges);
end

end
