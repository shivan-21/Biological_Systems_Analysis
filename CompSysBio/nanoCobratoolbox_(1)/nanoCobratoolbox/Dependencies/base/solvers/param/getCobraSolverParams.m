function varargout = getCobraSolverParams(solverType, paramNames, parameters)
% This function gets the specified parameters in `paramNames` from
% parameters, the global cobra paramters variable or default values set within
% this script. 
%
% It will use values with the following priority
% parameters > solver type parameters > default parameters
%
% The specified parameters will be delt to the specified output arguements.
% See examples below.
%
% USAGE:
%
%    varargout = getCobraSolverParams(solverType, paramNames, parameters)
%
% INPUTS:
%    solverType:    Type of solver used: 'LP', 'MILP', 'QP', 'MIQP'
%    paramNames:    Cell array of strings containing parameter names OR one
%                   parameter name as string
%
% OPTIONAL INPUTS:
%    parameters:    Structure with fields pertaining to parameter values that
%                   should be used in place of global or default parameters.
%                   parameters can be set to 'default' to use the default
%                   values set within this script.
%
% OUTPUTS:
%    varargout:     Variables which each value corresponding to paramNames
%                   is output to.
%
% EXAMPLE:
%    parameters.saveInput = 'LPproblem.mat';
%    parameters.printLevel = 1;
%    [printLevel, saveInput] = getCobraSolverParams('LP', {'printLevel', 'saveInput'}, parameters);
%
%    %Example using default values
%    [printLevel, saveInput] = getCobraSolverParams('LP', {'printLevel','saveInput'}, 'default');
%
% .. Authors:
%       - Richard Que (12/01/2009)
%       - Ronan (16/07/2013) default MPS parameters are no longer global variables

if nargin < 2
    error('getCobraSolverParams: No parameters specified')
end
if nargin <3
    parameters = '';
end
% Persistence will make names specific to one type of solver.
% Default Values
% For descriptions of the different settings please have a look at 
% getCobraSolverParamsOptionsForType
valDef.minNorm = 0;

%These default tolerances are based on the default values for the Gurobi LP solver
%https://www.gurobi.com/documentation/9.0/refman/parameters.html
valDef.feasTol = 1e-6; % (primal) feasibility tolerance
valDef.optTol = 1e-6;  % (dual) optimality tolerance

valDef.objTol = 1e-6; % this should be used only when comparing the values of two objective functions

valDef.minNorm = [];

valDef.printLevel = 0;
valDef.verify = 0;

valDef.primalOnly = 0;

valDef.timeLimit = 1e36; 
valDef.iterationLimit = 1000;

%valDef.logFile = ['Cobra' solverType 'Solver.log'];
valDef.logFile = []; %log file should be empty to avoid creating it by default
valDef.saveInput = [];
valDef.PbName = [solverType 'problem'];
valDef.debug = 0;
valDef.lifting = 0;
   
valDef.method = -1; 

% CPLEX parameters
valDef.DATACHECK = 1;
valDef.DEPIND = 1;
valDef.checkNaN = 0;
valDef.warning = 0;

% tolerances
valDef.intTol = 1e-12;
valDef.relMipGapTol = 1e-12;
valDef.absMipGapTol = 1e-12;
valDef.NUMERICALEMPHASIS = 1;    

if ~iscell(paramNames)
    paramNames = {paramNames};
end

switch solverType
    case 'LP'
        global CBT_LP_PARAMS
        parametersGlobal = CBT_LP_PARAMS;
    case 'QP'
        global CBT_QP_PARAMS
        parametersGlobal = CBT_QP_PARAMS;
    case 'MILP'
        global CBT_MILP_PARAMS
        parametersGlobal = CBT_MILP_PARAMS;
    case 'EP'
        global CBT_EP_PARAMS
        parametersGlobal = CBT_EP_PARAMS;
        valDef.feasTol = 1e-6; % (primal) feasibility tolerance
        valDef.optTol = 1e-6;  % (dual) optimality tolerance
    case 'MIQP'
        global CBT_MIQP_PARAMS
        parametersGlobal = CBT_MIQP_PARAMS;
    case 'NLP'
        global CBT_NLP_PARAMS
        parametersGlobal = CBT_NLP_PARAMS;
    otherwise
        display('Unrecognized solver type')
        return;
end

varargout = cell(1, numel(paramNames));
paramNames = columnVector(paramNames);
for i=1:length(paramNames)
    % set values to default
    if isfield(valDef,paramNames{i})
        varargout{i} = valDef.(paramNames{i});
    end
    if ~strcmp(parameters,'default') % skip of using default values
        % set values to global values
        if isfield(parametersGlobal,paramNames{i})
           varargout{i} = parametersGlobal.(paramNames{i});
        end
        % set values to specified values
        if isfield(parameters,paramNames{i})
           varargout{i} = parameters.(paramNames{i});
        end
    end
end
