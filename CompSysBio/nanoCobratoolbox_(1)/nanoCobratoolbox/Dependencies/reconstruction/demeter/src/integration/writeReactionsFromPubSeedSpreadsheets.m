function writeReactionsFromPubSeedSpreadsheets(infoFilePath,inputDataFolder,spreadsheetFolder)
% Prepares input file for the comparative genomics part
% Write reaction spreadsheets from InReactions and PubSEED spreadsheets.
%
% USAGE:
%   writeReactionsFromPubSeedSpreadsheets(infoFilePath,inputDataFolder,spreadsheetFolder)
%
% INPUTS
% infoFilePath          File with information on reconstructions to refine
% inputDataFolder       Folder to save propagated data to (default: folder 
%                       in current path called "InputData")                
% spreadsheetFolder     Folder with comparative genomics data retrieved 
%                       from PubSEED in spreadsheet format if available. 
%                       For an example of the required format, see 
%                       cobratoolbox/papers/2021_demeter/exampleSpreadsheets.
%
% .. Authors:
%       - Almut Heinken, 06/2020

% get PubSEED IDs of new organisms to reconstruct
infoFile = readtable(infoFilePath, 'ReadVariableNames', false);
infoFile = table2cell(infoFile);

% find folder with annotation versions in the information file
versionCol = find(strcmp(infoFile(1,:),'Annotation version ID'));

% load reactions
reactions=readtable('InReactions.txt', 'ReadVariableNames', false,'FileType','text','delimiter','tab');
reactions = table2cell(reactions);
exchanges=reactions(find(strcmp(reactions(:,1),'Exchange')),2);

reactionDatabase = readtable('ReactionDatabase.txt', 'Delimiter', 'tab','TreatAsEmpty',['UND. -60001','UND. -2011','UND. -62011'], 'ReadVariableNames', false);
reactionDatabase=table2cell(reactionDatabase);
database.reactions=reactionDatabase;
for i=1:length(exchanges)
    exchanges{i,2}=database.reactions{ismember(database.reactions(:, 1), exchanges{i,1}), 3};
end

% get all spreadsheets
dInfo = dir(spreadsheetFolder);
fileList={dInfo.name};
fileList=fileList';
fileList(find(strcmp(fileList(:,1),'.')),:)=[];
fileList(find(strcmp(fileList(:,1),'..')),:)=[];

genomeAnnotation={};
cnt=1;

for i=1:length(fileList)
    i
    spreadsheet=readtable([spreadsheetFolder filesep fileList{i}], 'ReadVariableNames', false,'FileType','text','delimiter','tab');
    spreadsheet = table2cell(spreadsheet);
    for j=2:size(spreadsheet,1)
        % replace PubSeed with AGORA Model IDs
        % some entries in the comparative genomics spreadsheet have no
        % reconstruction -> skip
        if ~isempty(find(strcmp(infoFile(:,2),spreadsheet{j,1})))
            spreadsheet{j,1}=infoFile{find(strcmp(infoFile(:,2),spreadsheet{j,1})),1};
            for k=2:size(spreadsheet,2)
                if ~isempty(spreadsheet{j,k})
                    getRxns=reactions(find(strcmp(reactions(:,1),spreadsheet{1,k})),2);
                    for r=1:length(getRxns)
                        genomeAnnotation{cnt,1}=spreadsheet{j,1};
                        genomeAnnotation{cnt,2}=getRxns{r};
                        % include GPRs
                        peg=infoFile(find(strcmp(infoFile(:,1),spreadsheet{j,1})),versionCol);
                        genes=strsplit(spreadsheet{j,k},',');
                        if strcmp(genes,'][')
                            genomeAnnotation{cnt,3}='gap_filled';
                        else
                        if length(genes) ==1  
                            genes{1}=genes{1}(~isspace(genes{1}));
                        genomeAnnotation{cnt,3}=[peg{1} '.' genes{1} '.peg'];
                        else
                            genomeAnnotation{cnt,3}=[peg{1} '.' genes{1} '.peg'];
                            for g=2:length(genes)
                                genes{g}=genes{g}(~isspace(genes{g}));
                                genomeAnnotation{cnt,3}=[genomeAnnotation{cnt,3} ' or ' peg{1} '.' genes{g} '.peg'];
                            end
                        end
                        end
                        cnt=cnt+1;
                        % get exchanges if needed
                        formula=database.reactions{ismember(database.reactions(:, 1), getRxns{r}), 3};
                        if contains(formula,'[e]')
                            mets=strsplit(formula,' ');
                            exMets=mets(contains(mets,'[e]'));
                            for t=1:length(exMets)
                                exRxn=['EX_' strrep(exMets{t},'[e]','(e)')];
                                genomeAnnotation{cnt,1}=spreadsheet{j,1};
                                genomeAnnotation{cnt,2}=exRxn;
                                genomeAnnotation{cnt,3}='gap_filled';
                                cnt=cnt+1;
                            end
                        end
                    end
                end
            end
        end
    end
end

% remove duplicate reactions for organisms
if size(genomeAnnotation,1)>0
delArray=[];
cnt=1;
orgs=unique(genomeAnnotation(:,1));
for i=1:length(orgs)
    getRxns=genomeAnnotation(find(strcmp(genomeAnnotation(:,1),orgs{i})),2);
    getRxnInds=find(strcmp(genomeAnnotation(:,1),orgs{i}));
    [uniqueRxns, ~, J]=unique(getRxns);
    cntRxns = histc(J, 1:numel(uniqueRxns));
    duplRxns=uniqueRxns(cntRxns>1);
    if ~isempty(duplRxns)
        for j=1:length(duplRxns)
            findInds=getRxnInds(find(strcmp(getRxns(:,1),duplRxns{j})),1);
            for k=2:length(findInds)
                delArray(cnt,1)=findInds(k);
                cnt=cnt+1;
            end
        end
    end
end
genomeAnnotation(delArray,:)=[];
end

writetable(cell2table(genomeAnnotation),[inputDataFolder filesep 'genomeAnnotation'],'FileType','text','WriteVariableNames',false,'Delimiter','tab');

end
