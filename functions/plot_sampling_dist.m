function [] = plot_sampling_dist(model_ref,x,y,rxns,idx)
%PLOT_SAMPLING_DIST Summary of this function goes here
%   Detailed explanation goes here

rxn_plotting = getIndexes(model_ref,rxns,"rxns");
full_x = full(x);
full_y = full(y);
% Test 
rxn_test_x = full_x(rxn_plotting(idx),:);
rxn_test_y = full_y(rxn_plotting(idx),:);
f = figure('visible','off');hold on;
plot1 = scatter(1:1000,rxn_test_x,'filled');
hold on;
plot2 = scatter(1:1000,rxn_test_y,'filled');
xlabel('Sampling iteration');
ylabel('flux (mmol/gDW/h)');
legend({'MCF7', 'MCF7-TAMr'}, 'Location', 'northeast');
out_file = append('test','.jpeg');
saveas(f,fullfile("./",out_file));
end

