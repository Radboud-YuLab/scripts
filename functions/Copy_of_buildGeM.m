% Load the base HUMAN-GEM
load('Human-GEM.mat');
setRavenSolver("gurobi");
ihuman = addBoundaryMets(ihuman);
essentialTasks = parseTaskList("C:\Users\t.schaefers\Documents\Human-GEM-main\Human-GEM-main\data\metabolicTasks\metabolicTasks_Essential.txt");
checkTasks(ihuman, [], true, false, false, essentialTasks);

% Test dataset from GTEX
tpms = {'C:\Users\t.schaefers\Documents\MeBo-project\TPMs\E-MTAB-4801-query-results-MDAMB231.txt',...
        'C:\Users\t.schaefers\Documents\MeBo-project\TPMs\E-MTAB-4801-query-results-MCF7.txt'};
    
expr_threshold = 0.1;    
project = 'E-MTAB-4801'

% Create two different models
for i = 1:length(tpms)
    % Create data struct
    expr_data = readtable(tpms{i});
    cell_line = expr_data.Properties.VariableNames(3:end);
    data_struct.genes = expr_data.GeneID;
    data_struct.tissues = cell_line;
    data_struct.levels = table2array(expr_data(:,3:end));
    data_struct.threshold = expr_threshold;
    % Create output file
    out_file = strcat(cell_line,"_",project,".mat");
    disp("Generating GEM: " + cell_line);
    disp("Expression data: " + tpms{i});
    disp("Genes > expression threshold: " + numel(data_struct.levels(data_struct.levels > expr_threshold)));
    % Setup model
    refModel = ihuman;  % the reference model from which the GEM will be extracted
    tissue = data_struct.tissues{1};  % must match the tissue name in data_struct.tissues
    celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
    hpaData = [];  % data structure containing protein abundance information (not used here)
    arrayData = data_struct;  % data structure with gene (RNA) abundance information
    metabolomicsData = [];  % list of metabolite names if metabolomic data is available
    removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
    taskFile = [];  % we already loaded the task file, so this input is not required
    useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
    printReport = true;  % (default) print status/completion report to screen
    taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
    params = [];  % additional optimization parameters for the INIT algorithm
    paramsFT = [];  % additional optimization parameters for the fit-tasks algorithm
    % Build reference model
    GEM_RAVEN = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);
    GEM_RAVEN.id = data_struct.tissues{1};
    GEM_RAVEN.name = data_struct.tissues{1};
    % Create cobra model as well
    GEM_COBRA = ravenCobraWrapper(GEM_RAVEN);
    GEM_COBRA.modelID = data_struct.tissues{1};
    save(out_file, 'GEM_RAVEN', 'GEM_COBRA');
end