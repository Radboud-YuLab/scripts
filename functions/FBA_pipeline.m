% toDo next
% 1) Run additional model for Taximofen mutant (apply correct O2 uptake constraints
% 2) Count how much ROX detox reaction are differentially changed between
% the two conditions
% Cobra stuff
initCobraToolbox;
changeCobraSolver('gurobi');
% Build models
build_models = true;
apply_seahore_constraints = true;
random_sampling = true;
nRandomSamples =1000;
tpm_threshold = 1;
% Model settings
essentialTasks = "C:\Users\t.schaefers\Documents\MetabolicModelling\resources\Human-GEM\data\metabolicTasks\metabolicTasks_Essential.txt";
%tpm_mdamb231 = 'C:\Users\t.schaefers\Documents\MetabolicModelling\data\E-MTAB-2770\E-MTAB-2770-MDAMB231.txt';
%tpm_mcf7 = 'C:\Users\t.schaefers\Documents\MetabolicModelling\data\E-MTAB-2770\E-MTAB-2770-MCF7.txt';
tpm_mcf7_wt = 'C:\Users\t.schaefers\Documents\MetabolicModelling\data\TPM_Paul\MCF7_WT_TPM.txt';
tpm_mcf7_tamr = 'C:\Users\t.schaefers\Documents\MetabolicModelling\data\TPM_Paul\MCF7_TAMr_TPM.txt';
medium_list = 'C:\Users\t.schaefers\Documents\MetabolicModelling\data\trace_nutrients_long.csv';

% Output settings
out_dir_root = "C:\Users\t.schaefers\Documents\MetabolicModelling\results\final-results-paul-latest2";
out_file = 'TPM_Paul_newest_original.mat';

% Create output directiroy if not existent
if ~exist(out_dir_root, 'dir')
   mkdir(out_dir_root);
end

% Convert Seahorse FX constraints for O2 Uptake
% pmol↔mmol 1 mmol = 1000000000 pmol.
convert = @(x) (1.0000e-12*x)/0.001*60;
% Load essential tasks & media constraints
missing_rxns = {"MAR03964","MAR06608","MAR04767","MAR12125","MAR08415","MAR08413","MAR08410","MAR08409","MAR03982","MAR03980","MAR03960"};
% Build metabolic models
if build_models == true
    % E-MTAB-2770
    %mdamb231 = buildGEM(tpm_mdamb231,essentialTasks, tpm_threshold);
    mcf7_wt = buildGEM(tpm_mcf7_wt,essentialTasks, tpm_threshold);
    mcf7_tamr = buildGEM(tpm_mcf7_tamr,essentialTasks, tpm_threshold);
    save(fullfile(out_dir_root,out_file),'mcf7_wt','mcf7_tamr');
end
% Load base models (without constraints applied)
model_data = load(fullfile(out_dir_root,out_file));
model_ref = model_data.mcf7_wt;
model_alt = model_data.mcf7_tamr;
%%%%% Create MCDF7-parental model %%%%%
baseMedium = readtable(medium_list);
if apply_seahore_constraints == true
    idx = baseMedium.MET_ID == "MAM02630e";
    MCF7_FCCP = {399.32; 420.50 ;432.01; 426.34}; %Maxiumum oxygen respiration (pmol/min)
    MCF7_ANROT = {47.54; 49.70; 48.19; 45.90}; % Minimum oxygen respiration (pmol/min)
    % Set bounds for oxygen exchange
    baseMedium(idx,"lb") = {-max(cellfun(convert,MCF7_FCCP))};
    baseMedium(idx,"ub") = {-min(cellfun(convert,MCF7_ANROT))};
end
model_mcf7 = prep_model(model_ref, missing_rxns, essentialTasks, baseMedium);
sol_mcf7 = solveLP(model_mcf7); % FBA analysis
%%%%%%% Create MCF7-tamoxifen model %%%%%%
baseMedium = readtable(medium_list);
if apply_seahore_constraints == true
    idx = baseMedium.MET_ID == "MAM02630e";
    MCF7_FCCP_TM = {260.85; 255.33; 237.61; 217.22}; %Maxiumum oxygen respiration (pmol/min)
    MCF7_ANROT_TM = {43.25; 41.90; 39.68; 37.98}; %Minimum oxygen respiration (pmol/min)
    baseMedium(idx,"lb") = {-max(cellfun(convert,MCF7_FCCP_TM))};
    baseMedium(idx,"ub") = {-min(cellfun(convert,MCF7_ANROT_TM))};
end
model_mcf7_tm = prep_model(model_alt, missing_rxns, essentialTasks, baseMedium);
sol_mcf7_tm = solveLP(model_mcf7_tm);
% Save intermediate analysis
save(fullfile(out_dir_root,'AfterConstraintsApplied.mat'));
%%%% Perform random sampling for mdamb231 %%%%%
if random_sampling
    %MCF7 sampling
    [~,good_rxns_mcf7]=randomSampling(model_mcf7,nRandomSamples,true,false,true,[],false);
    solutions_mcf7=randomSampling(model_mcf7,nRandomSamples,true,false,true,good_rxns_mcf7,false);
    %MDAMB231 random sampling
    [~,good_rxns_mcf7_tm]=randomSampling(model_mcf7_tm,nRandomSamples,true,false,true,[],false);
    solutions_mcf7_tm=randomSampling(model_mcf7_tm,nRandomSamples,true,false,true,good_rxns_mcf7_tm,false);
    % Compare stats
    res = evaluateSampling(model_mcf7, model_mcf7_tm, solutions_mcf7, solutions_mcf7_tm, 0.01,0,out_dir_root);
    res_subs = countAltSubsystems(model_ref,res, out_dir_root);
    res_subs2 = countAltSubsystems(model_alt,res, out_dir_root);
    %Create a bar chart of differential reactions (normalized by metabolic
    save(fullfile(out_dir_root,'AfterRandomSampling.mat'));
end
% Plot random sampling results
% res = res(ismember(res.rxn,{'MAR04757';'MAR03980';'MAR03939';'MAR04657'}),:);
% for i=1:length(res.rxn);
%     %f = figure('visible','off','WindowState','minimized');
%     target = res.rxn{i};
%     stat = res.ks_pVal(i);
%     rxnName = model_mcf7.rxnNames(getIndexes(model_mcf7,target,"rxns"))
%     plotSampleHist({target}, {full(solutions_mcf7),full(solutions_mcf7_tm)}, {model_mcf7, model_mcf7_tm})
%     legend({'MCF7', 'MCF7-TAMr'}, 'Location', 'Best')
%     title(append(rxnName,' (',vertcat(res.subsystem{i}),')'));
%     subtitle(append('{\it p}=',vertcat(string(stat)),',','{\it n}=',string(nRandomSamples)));
%     out_file = append(string(target),'.pdf');
%     print(fullfile(out_dir_root,'plots',out_file),'-dpdf');
% end