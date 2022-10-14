function [model_raven, model_cobra] = buildGEM(tpm, task_list, expr_threshold)

% Load the base HUMAN-GEM
ihuman = importYaml('Human-GEM.yml');
ihuman = addBoundaryMets(ihuman);
setRavenSolver("gurobi");
 
% Load essential TASK list
essentialTasks = parseTaskList(task_list);
% Create data struct
expr_data = readtable(tpm);
vars = expr_data.Properties.VariableNames;
% Replace NaN's with 0
expr_data.(vars{2})(isnan(expr_data.(vars{2}))) = 0;
% Build struct
cell_line = expr_data.Properties.VariableNames(2);
data_struct.genes = expr_data.(vars{1});
data_struct.tissues = cell_line;
data_struct.levels = table2array(expr_data(:,2:end));
data_struct.threshold = expr_threshold;
% Create output file
disp("Generating GEM: " + cell_line);
disp("Expression data: " + tpm);
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
model_raven = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);
model_raven.id = data_struct.tissues{1};
model_raven.name = data_struct.tissues{1};
% Create cobra model as well
model_cobra = ravenCobraWrapper(model_raven);
model_cobra.modelID = data_struct.tissues{1};

end