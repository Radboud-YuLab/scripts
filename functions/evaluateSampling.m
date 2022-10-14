function [res] = evaluateSampling(model_ref, model, x, y, alpha, eff_thres, out_dir)
%EVALUATESAMPLING Summary of this function goes here
%   Detailed explanation goes here
format long
x=full(x);
y=full(y);
%Calculate set of overlapping reactions
rxn_overlap = intersect(model_ref.rxns, model.rxns);
disp("Total reactions ref/alt: " + numel(model_ref.rxns)+'/'+numel(model.rxns));
disp("Reaction overlap: " + numel(rxn_overlap));
x1 = x(getIndexes(model_ref, rxn_overlap, "rxns"),:);
y1 = y(getIndexes(model, rxn_overlap, "rxns"),:);
%Compare sample stats
[stats, pVals] = compareTwoSamplesStat(x1,y1);
%Caluclate effect size
eff_sizes = cell(length(rxn_overlap),1);
for i=1:length(x1)
    eff_sizes{i,1} = abs(computeCohen_d(y1(i,:),x1(i,:)));
end
% select stats
idx = pVals.ks < alpha & abs(cell2mat(eff_sizes)) > eff_thres;
ks_pvals = pVals.ks(idx);
ks_stats = stats.ks(idx);
eff_filtered = eff_sizes(idx);
disp("Significantly altered reactions: " + numel(ks_pvals));
% Select significant reactions
rxn_changed = rxn_overlap(idx,:);
rxn_idx_changed = getIndexes(model_ref,rxn_changed,"rxns");
subsystem_changed =  model_ref.subSystems(rxn_idx_changed);
disp(subsystem_changed);
% Create final result table
stats_x = calcSampleStats(x1(idx,:));
stats_y = calcSampleStats(y1(idx,:));
% Write final results table
res = table(rxn_changed, subsystem_changed, stats_x.mean, stats_y.mean, stats_x.std, stats_y.std, ks_stats, ks_pvals, eff_filtered);
res.Properties.VariableNames = {'rxn','subsystem','mean_sample_x','mean_sample_y','std_sample_x','std_sample_y','ks_stats','ks_pVal','cohen_effect'};
res.Properties.VariableNames{'mean_sample_x'} = append('mean_sample_',model_ref.id);
res.Properties.VariableNames{'mean_sample_y'} = append('mean_sample_',model.id);
res.Properties.VariableNames{'std_sample_x'} = append('std_sample_',model_ref.id);
res.Properties.VariableNames{'std_sample_y'} = append('std_sample_',model.id);
res = sortrows(res,8,'ascend');
% return res and idx of changes reactions
writetable(res,fullfile(out_dir,'SamplingRxnsAltered.csv'),'Delimiter',',');  

end

