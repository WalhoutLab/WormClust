
%% Loading model table
addpath('/data/nandas/Coflux_matrix/COBRA/cobratoolbox')
initCobraToolbox;
changeCobraSolver ('gurobi', 'all', 1);
%load('final_model_1315_COBRA.mat')
parpool(10);
%% Loading the network & calculating flux
model=iCEL2model_xl("Final_model_table.xlsx")
model=convertToIrreversible(model)
iCel1315_model=model
%% Gene to reaction mapping 
Tbl= loadNetwork("Final_model_table.xlsx","Reaction List",1);
Tbl=Tbl([1:3022],[1:9])
rxn2geneMapping = Tbl([1:3022],[1,5]);
G2RMap = containers.Map;
for i = 1:length(rxn2geneMapping);
   rxn = rxn2geneMapping{i,1};
   rxnIndexfun = cellfun(@(x)isequal(x,rxn),iCel1315_model.rxns);
   rxnIndex = find(rxnIndexfun);
   genes = rxn2geneMapping{i,2};
   genes = strrep(genes, '(','');
   genes= strrep(genes, ')', '');
   genes= strrep(genes, '[', '');
   genes= strrep(genes, ']', '');
   genes= strrep(genes, '|', ',');
   genes= strrep(genes, '&', ',');
   genes = strsplit(genes, ',');
   for j = 1:length(genes);
       gene = strtrim(genes{j});
       if(isKey(G2RMap, gene))
           reactionsIndexList = G2RMap(gene);
           reactionsIndexList = [reactionsIndexList, rxnIndex];
           G2RMap(gene) = reactionsIndexList
       else
           G2RMap(gene) = rxnIndex;
       end
   end   
end

%% Genes List
genesList = G2RMap.keys;
%lb_original = M.LB;
%ub_original = M.UB;
Vmax={}

%% Computing Vmax(i)
for i=1:length(iCel1315_model.rxns);
    iCel1315_model=model;
    iCel1315_model=changeObjective(iCel1315_model,iCel1315_model.rxns(i));
    V=optimizeCbModel(iCel1315_model,'max');
    %V=changeRxnBounds(iCel1315_model,iCel1315_model.rxns(i),V.obj,'l')
    Vmax{i}=V.v(i);
    i
end
iCel1315_model=model;
M_original=model;
%% Computing Vmax(i,j)
correlation_matrix = {};
%environment = getEnvironment();
%restoreEnvironment(environment);
pctRunOnAll("initCobraToolbox(false)");
%parfor j = 1:length(iCel1315_model.rxns);
parfor j = 1:length(iCel1315_model.rxns);
    iCel1315_model=model;
    %Updating the lowerbound and upperBound to 0 to make sure Jth reaction
    %is set to 0
    j
    iCel1315_model.lb = M_original.lb;
    iCel1315_model.lb(j) = 0;
    iCel1315_model.ub = M_original.ub;
    iCel1315_model.ub(j) = 0;
    rxnNamej = iCel1315_model.rxns{j}(1:end-1);
    rxnTypej = iCel1315_model.rxns{j}(end);
    if(rxnTypej == 'f')
        name=strcat(rxnNamej,'r');
    elseif(rxnTypej == 'r')
        name=strcat(rxnNamej,'f');
    end
    %iCel1315_model = changeRxnBounds(iCel1315_model, name, 0);
    iCel1315_model.lb(strcmp(name,iCel1315_model.rxns))=0;
    iCel1315_model.ub(strcmp(name,iCel1315_model.rxns))=0;
    %computing FluxBalanceAnalysis for all the reactions given Jth reaction
    %is set to 0
    submatrix = cell(length(iCel1315_model.rxns),1);
    %submatix2 = cell(length(iCel1315_model.rxns),1);
    for i = 1:length(iCel1315_model.rxns);
        fprintf("MAximizing i=%.4f and reaction=%s\n Deleting j=%.4f and reaction=%s",i,iCel1315_model.rxns{i},j,iCel1315_model.rxns{j})
        rxnName = iCel1315_model.rxns{i}(1:end-1);
        rxnType = iCel1315_model.rxns{i}(end);
        if(rxnType == 'f')
            name=strcat(rxnName,'r');
        elseif(rxnType == 'r')
            name=strcat(rxnName,'f');
        end
        %iCel1315_model = changeRxnBounds(iCel1315_model, name, 0);
        iCel1315_model.lb(strcmp(name,iCel1315_model.rxns))=0;
        iCel1315_model.ub(strcmp(name,iCel1315_model.rxns))=0;
        iCel1315_model=changeObjective(iCel1315_model,iCel1315_model.rxns(i));
        V=optimizeCbModel(iCel1315_model,'max');
        %F = fba(M, reactions{i}, 'max');
        if (isfield(V,"v"))
            submatrix{i} = (Vmax{i} - V.v(i))/Vmax{i};
        else
            submatrix{i} = (Vmax{i} - 0)/Vmax{i};
            fprintf("v not found")
        end%submatix2{i}=V.v;
        iCel1315_model.lb(strcmp(name,iCel1315_model.rxns))=model.lb(strcmp(name,model.rxns));
        iCel1315_model.ub(strcmp(name,iCel1315_model.rxns))=model.ub(strcmp(name,model.rxns));
    end
    correlation_matrix(:,j) = submatrix;
    %max_flux_matrix(:,j)=submatix2;
end
csvwrite('reactions_extended1315_fluxred.csv',correlation_matrix);
save('reaction_extended1315_coflux.mat','correlation_matrix');
reactions=iCel1315_model.rxns;
save('extended_reactions.mat','reactions');
%save('reactions_maxflux.mat','max_flux_matrix')


%% Gene Correlation matrix
geneCorrMatrix = {};
%geneCorrMatrix2={};
for i = 1:length(genesList);
    for j = 1:length(genesList);
        rxnList_i = values(G2RMap,{genesList{i}});
        rxnList_j = values(G2RMap,{genesList{j}});
        rxnCorrValues = {}
        rxnCorrValues2={}
        for k = 1:length(rxnList_i);
            for l = 1:length(rxnList_j);
                rxnCorrValues{end+1} = correlation_matrix{rxnList_i{k}, rxnList_j{l}};
                %rxnCorrValues2{end+1} = max_flux_matrix{rxnList_i{k}, rxnList_j{l}};
            end
        end
        geneCorrMatrix{i,j} = max(rxnCorrValues{:});
        %geneCorrMatrix2{i,j} = max(rxnCorrValues2{:});
    end
end

save('gene_extended_1315_coflux.mat','geneCorrMatrix');
save('genesList_iCel1315.mat','genesList');
%save('gene_extended_1315_maxflux.mat','geneCorrMatrix2');
csvwrite('gene_extended_1315_fluxred.csv',geneCorrMatrix);
%csvwrite('gene_extended_1315_maxflux.csv','geneCorrMatrix2');
csvwrite('reactions_extended1315_fluxred.csv',correlation_matrix);
%csvwrite('reactions_extended1315_maxflux.csv','max_flux_matrix')

