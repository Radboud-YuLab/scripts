function model = addrxnBack(model,model_original,rxnsList)
if isfield(model,'rxnMiriams')
    model = ravenCobraWrapper(model);
end
for j = 1:length(rxnsList)
rxn.rev =  model_original.lb(rxnsList(j));
    if rxn.rev == -1000
        rxn.rev = 1;
    else
        rxn.rev = 0;
    end
    mets = find(model_original.S(:,rxnsList(j)));
    matrix.coef = model_original.S(mets,rxnsList(j));
    matrix.metID = model_original.mets(mets);
    matrix.metNames = model_original.metNames(mets);
    %mapping mets to model.metnames, get s_ index for new mets
    [~,metindex] = ismember(matrix.metNames,model.metNames);
    if any(metindex == 0)
        for k = 1:length(metindex)
            if metindex(k) == 0
                model = addMetabolite(model,matrix.metID{k}, ...
                    'metName',matrix.metNames{k});
            end
        end
    end
    %add reaction back
    rxnID   = model_original.rxns{rxnsList(j)};
    model = addReaction(model,...
        rxnID,...
        'reactionName', rxnID,...
        'metaboliteList',matrix.metID,...
        'stoichCoeffList',matrix.coef,...
        'reversible',rxn.rev);
    disp(['rxn: ',rxnID, ' has been added back to the model '])
    printRxnFormula(model,'rxnAbbrList',rxnID,'metNameFlag',true);
end
end