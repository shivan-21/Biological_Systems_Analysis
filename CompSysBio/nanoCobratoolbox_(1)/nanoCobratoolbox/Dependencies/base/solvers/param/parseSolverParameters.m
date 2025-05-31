function [problemTypeParams, solverParams] = parseSolverParameters(problemType, varargin)
% Gets default cobra solver parameters for a problem of type problemType, unless 
% overridden by cobra solver parameters provided by varagin either as parameter
% struct, or as parameter/value pairs.
%
% USAGE:
%    [problemTypeParams, solverParams] = parseSolverParameters(problemType,varargin)
%
% INPUT:
%    problemType:       The type of the problem to get parameters for
%                       ('LP','MILP','QP','MIQP','NLP')
%
% OPTIONAL INPUTS:
%    varargin:          Additional parameters either as parameter struct, or as
%                       parameter/value pairs. A combination is possible, if
%                       the parameter struct is either at the beginning or the
%                       end of the optional input.
%                       All fields of the struct which are not COBRA parameters
%                       (see `getCobraSolverParamsOptionsForType`) for this
%                       problem type will be passed on to the solver in a
%                       solver specific manner.
%
% OUTPUTS:
%    problemTypeParams: The COBRA Toolbox specific parameters for this
%                       problem type given the provided parameters
%
%    solverParams:      Additional parameters provided which are not part
%                       of the COBRA parameters and are assumed to be part
%                       of direct solver input structs.

cobraSolverParameters = getCobraSolverParamsOptionsForType(problemType); % build the default Parameter Structure

% set the solver Type
eval(['global CBT_' problemType '_SOLVER;'])
eval(['defaultSolver = CBT_' problemType '_SOLVER;']);

% initialize the solver variables
solverVars = cell(numel(cobraSolverParameters),1);

% get the default variables for the correct solver.
[solverVars{:}] = getCobraSolverParams(problemType,cobraSolverParameters,struct('solver',defaultSolver));
defaultParams = [cobraSolverParameters',solverVars];

nVarargin = numel(varargin);
% parse the supplied parameters
if nVarargin > 0
    % we should have a struct at the end
    if mod(nVarargin,2) == 1
        optParamStruct = varargin{end};
        if ~isstruct(optParamStruct)
            % but it could also be at the first position, so test that as well.
            optParamStruct = varargin{1};
            varargin(1) = [];
            nVarargin = numel(varargin); %added this in case varagin{1} is the parameter structure
            if ~isstruct(optParamStruct)
                error(['Invalid Parameters supplied.\n',...
                       'Parameters have to be supplied either as parameter/Value pairs, or as struct.\n',...
                       'A combination is possible, if the last or first input argument is a struct, and all other arguments are parameter/value pairs'])
            end
        else
            varargin(end) = [];
        end
    else
        % no parameter struct. so initialize an empty one.
        optParamStruct = struct();
    end
    nVarargin = numel(varargin); %added this in case varagin{1} is the parameter structure
    % now, loop through all parameter/value pairs.
    for i = 1:2:nVarargin
        cparam = varargin{i};
        if ~ischar(cparam)
            error('Parameters have to be supplied as ''parameterName''/Value pairs');
        end
        % the param struct overrides the only use the parameter if it is
        % not a field of the parameter struct.
        if ~isfield(optParamStruct,cparam)
            try
                optParamStruct.(cparam) = varargin{i+1};
            catch
                error('All parameters have to be valid matlab field names. %s is not a valid field name',cparam);
            end
        else
            warning('Duplicate parameter %s, both supplied as a field name and a parameter/value pair!',cparam);
        end
    end
else
    % no potional parameters.
    optParamStruct = struct();
end

% set up the cobra parameters
problemTypeParams = struct();

for i = 1:numel(defaultParams(:,1))
    % if the field is part of the optional parameters (i.e. explicitly provided) use it.
    if isfield(optParamStruct,defaultParams{i,1})
        problemTypeParams.(defaultParams{i,1}) = optParamStruct.(defaultParams{i,1});
        % and remove the field from the struct for the solver specific parameters.
        optParamStruct = rmfield(optParamStruct,defaultParams{i,1});
    else
        % otherwise use the default parameter
        problemTypeParams.(defaultParams{i,1}) = defaultParams{i,2};
    end
end

% assign all remaining parameters to the solver parameter struct.
solverParams = optParamStruct;