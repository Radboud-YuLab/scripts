function [model] = prep_model(model, missing_rxns, task_list, media)

ihuman = importYaml('Human-GEM.yml');
setRavenSolver("gurobi");

% Prepare the reference model (ihuman)
refModel = ihuman;
refModel = addBoundaryMets(refModel);

% Check if both models can perform essential tasks
essentialTasks = parseTaskList(task_list);
checkTasks(model, [], true, false, false, essentialTasks);

% Simplify model
model = simplifyModel(model);

% Correct models for missing ATP hydrolisis/ ROS detox reactions
missing_rxn_idxs = getIndexes(refModel,missing_rxns,"rxns");
model = addrxnBack(model,refModel, missing_rxn_idxs);

% Show the biomass objective
printObjective(model);

% Filter only relevant metabolites;
media_mets = table2array(media(:,1));
media_lb = table2array(media(:,2));
media_ub = table2array(media(:,3));
model = setExchangeBounds(model,media_mets,media_lb,media_ub,true,true);

end
