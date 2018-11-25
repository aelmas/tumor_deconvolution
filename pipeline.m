%% 1. Load data
load data_MC3.mat
load data_mRNA.mat
load data_PRO.mat

%% 2. Clean data
raw_values_mRNA = table2array(data_mRNA);
raw_values_PRO = table2array(data_PRO);
flag_normalize_mRNA = false;
flag_normalize_PRO = false;
% 2.a. Impute data
if ~isempty(find(isnan(raw_values_mRNA), 1)) 
    tic; values_mRNA = knnimpute(raw_values_mRNA); toc;
    flag_normalize_mRNA = true;
else
    values_mRNA = raw_values_mRNA;
end
if ~isempty(find(isnan(raw_values_PRO), 1))
    tic; values_PRO = knnimpute(raw_values_PRO); toc;
    flag_normalize_PRO = true;
else
    values_PRO = raw_values_PRO;
end

% 2.b. Remove outlier samples
c_dispersion = 3; %dispersion constant
non_outliers_mRNA = iqr(values_mRNA) <= mean(values_mRNA)+c_dispersion*std(values_mRNA);
values_mRNA = values_mRNA(:,non_outliers_mRNA);
non_outliers_PRO = iqr(values_PRO) <= mean(values_PRO)+c_dispersion*std(values_PRO);
values_PRO = values_PRO(:,non_outliers_PRO);
 
% 2.c. Normalize samples
% mRNA | normalized across samples at 75%% quantile | log2 normalized
if flag_normalize_mRNA
    values_mRNA = values_mRNA-repmat(mean(values_mRNA),size(values_mRNA,1),1); % 0-mean
    values_mRNA = values_mRNA./repmat(std(values_mRNA),size(values_mRNA,1),1); % 1-SD
end
% PRO | samples normalized | 0-mean, 1-SD | 
if flag_normalize_PRO
    values_PRO = values_PRO-repmat(mean(values_PRO),size(values_PRO,1),1); % 0-mean
    values_PRO = values_PRO./repmat(std(values_PRO),size(values_PRO,1),1); % 1-SD
end
% Check distributions
figure(1);
subplot(211); boxplot(values_mRNA); ylabel('BRCA Transcriptome');
subplot(212); boxplot(values_PRO); ylabel('BRCA Proteome');


%% 3. Deconvolution methods
%
%% 3.a. Tumor transcriptome
% xCell, performs cell type enrichment analysis from gene expression data for 64 
%    immune and stroma cell types by reducing associations between closely related cell types.
% xCell produces enrichment scores, not percentages.
% xCell uses the expression levels ranking and not the actual values, 
%    thus normalization does not have an effect.
% installation\\$ devtools::install_github('dviraran/xCell')
%
% pipe data to xCell
data_mRNA_cleaned = array2table(values_mRNA);
data_mRNA_cleaned.Properties.VariableNames = data_mRNA.Properties.VariableNames;
data_mRNA_cleaned.Properties.RowNames = data_mRNA.Properties.RowNames;
writetable(data_mRNA_cleaned,'BRCA_mRNA_formatted_normalized_cleaned.txt','Delimiter','\t','WriteRowNames',true);
% script to call xCell 
fid = fopen('call_xCell_mRNA.R','w');
fprintf(fid,[
    'library(xCell)\n' ...
    'expression_matrix = read.table("BRCA_mRNA_formatted_normalized_cleaned.txt",header=TRUE,row.names=1, as.is=TRUE)\n' ... 
    'xCell_result = xCellAnalysis(expression_matrix'...
    ', parallel.type = "FORK"' ... % faster in unix
    ', cell.types.use = NULL)\n' ... % use all cell types
    'write.table(xCell_result,file = "xCell_result_BRCA_mRNA_formatted_normalized_cleaned.txt",sep="\t",quote=FALSE)\n' ... 
    ]);
fclose(fid); 
% grand user permission of file executions
system('chmod u+x call_xCell_mRNA.R xCell-master/R/xCell.R') 
% call xCell
system('/usr/local/bin/Rscript call_xCell_mRNA.R');
% read xCell result
data_mRNA_cleaned_xCell_result = readtable('xCell_result_BRCA_mRNA_formatted_normalized_cleaned.txt','ReadRowNames',true);
% xCell has a bug that assigns the first variable to the row names | fix this
data_mRNA_cleaned_xCell_result.Properties.VariableNames = data_mRNA_cleaned.Properties.VariableNames; 
% display cell enrichments
figure(2); 
heatmap(data_mRNA_cleaned_xCell_result.Properties.VariableNames,data_mRNA_cleaned_xCell_result.Properties.RowNames,table2array(data_mRNA_cleaned_xCell_result))
title('BRCA Transcriptome cell enrichments')

%% 3.b. Tumor proteome
% a quick observation shows that there are ~10.5K common genetic elements between the analyzed transcriptome and proteome data sets:
length(intersect(data_mRNA.Properties.RowNames, data_PRO.Properties.RowNames))
% a deconvolution method designed for transcriptome data may worth studying
%    the proteome data that contains many common genetic signatures with transcriptome
%
% pipe data to xCell
data_PRO_cleaned = array2table(values_PRO);
data_PRO_cleaned.Properties.VariableNames = data_PRO.Properties.VariableNames;
data_PRO_cleaned.Properties.RowNames = data_PRO.Properties.RowNames;
writetable(data_PRO_cleaned,'BRCA_PRO_formatted_normalized_cleaned.txt','Delimiter','\t','WriteRowNames',true);
% script to call xCell 
fid = fopen('call_xCell_PRO.R','w');
fprintf(fid,[
    'library(xCell)\n' ... 
    'expression_matrix = read.table("BRCA_PRO_formatted_normalized_cleaned.txt",header=TRUE,row.names=1, as.is=TRUE)\n' ... 
    'xCell_result = xCellAnalysis(expression_matrix'...
    ', parallel.type = "FORK"' ... % faster in unix
    ', cell.types.use = NULL)\n' ... % use all cell types
    'write.table(xCell_result,file = "xCell_result_BRCA_PRO_formatted_normalized_cleaned.txt",sep="\t",quote=FALSE)\n' ... 
    ]);
fclose(fid); 
% grand user permission of file executions
system('chmod u+x call_xCell_PRO.R xCell-master/R/xCell.R') 
% call xCell
system('/usr/local/bin/Rscript call_xCell_PRO.R');
% read xCell result
data_PRO_cleaned_xCell_result = readtable('xCell_result_BRCA_PRO_formatted_normalized_cleaned.txt','ReadRowNames',true);
% xCell has a bug that assigns the first variable to the row names | fix this
data_PRO_cleaned_xCell_result.Properties.VariableNames = data_PRO_cleaned.Properties.VariableNames; 
% display cell enrichments
figure(3); 
heatmap(data_PRO_cleaned_xCell_result.Properties.VariableNames,data_PRO_cleaned_xCell_result.Properties.RowNames,table2array(data_PRO_cleaned_xCell_result))
title('BRCA Proteome cell enrichments')

%% 4. Analyze tumor deconvolution from transcriptome vs from proteome
% test the null hypothesis (H) that the two immune profiles inferred from
%    transcriptome and from proteome are not significantly different for the
%    same group of samples, are coming from distributions with equal means
% map sample IDs
[~, im_xv, ip_xv] = intersect(data_mRNA_cleaned_xCell_result.Properties.VariableNames,data_PRO_cleaned_xCell_result.Properties.VariableNames);
% H = 1 means that the null hypothesis can be rejected at 5% significance level
H_cell_types = ttest2(table2array(data_mRNA_cleaned_xCell_result(:,im_xv))', table2array(data_PRO_cleaned_xCell_result(:,ip_xv))')';
fprintf('|\n|\tFor the same group of samples, %u out of %u cell-types inferred differently.\n|\n', sum(H_cell_types), length(H_cell_types));
data_xCell_result_similar_mRNA_PRO = data_mRNA_cleaned_xCell_result.Properties.RowNames(H_cell_types==0);
% calculate correlations
cros_corr_cell_similar_mRNA_PRO = corr(table2array(data_PRO_cleaned_xCell_result(H_cell_types==0,ip_xv))',table2array(data_mRNA_cleaned_xCell_result(H_cell_types==0,im_xv))');
[corr_vals, ic_s] = sort(abs(diag(cros_corr_cell_similar_mRNA_PRO)),'descend');
ind_H0 = find(H_cell_types==0);
% The following 14 cell-types are inferred similarly when using transcriptome or proteome data
data_xCell_result_similar_mRNA_PRO(ic_s)

figure(4);
heatmap(data_mRNA_cleaned_xCell_result.Properties.VariableNames(im_xv), data_mRNA_cleaned_xCell_result.Properties.RowNames(ind_H0(ic_s)), table2array(data_mRNA_cleaned_xCell_result(ind_H0(ic_s),im_xv)))
title('BRCA Transcriptome, enrichments for 14 cells')
figure(5);
heatmap(data_PRO_cleaned_xCell_result.Properties.VariableNames(ip_xv), data_PRO_cleaned_xCell_result.Properties.RowNames(ind_H0(ic_s)), table2array(data_PRO_cleaned_xCell_result(ind_H0(ic_s),ip_xv)))
title('BRCA Proteome, enrichments for 14 cells')

%% 5. Proteome enrichment patterns due to mutation
% 5.a. Analyze proteomes of mutation carriers vs non-carriers in key driver genes
%    for each key driver gene: define carriers and non-carriers, compare their proteomes
%    then build protein profile for that mutation showing significance between carriers vs non-carriers
protein_profiles = cell(size(data_MC3,1),1);
signif_mut = false(1,size(data_MC3,1));
for k = 1:size(data_MC3,1) %for each mutation gene
    % define carriers and non-carriers of the mutation in mutation_gene(k) in MC3 data
    carrier_true = table2array(data_MC3(k,:))~='wt';
    % map those carriers/non-carriers in proteome data
    [~, ip_carrier, ~] = intersect(data_PRO_cleaned.Properties.VariableNames, data_MC3.Properties.VariableNames(carrier_true));
    [~, ip_noncarrier, ~] = intersect(data_PRO_cleaned.Properties.VariableNames, data_MC3.Properties.VariableNames(~carrier_true));
    if min(length(ip_carrier),length(ip_noncarrier)) > 1
        signif_mut(k) = true;
        % test the null hypothesis (H) that the proteomes of the mutation carriers and non-carriers in mutation_gene(k) are not significantly different
        % H = 1 means that the null hypothesis can be rejected at 5% significance level
        H_proteins = ttest2(table2array(data_PRO_cleaned(:,ip_carrier))', table2array(data_PRO_cleaned(:,ip_noncarrier))')';
        % any protein whose expression significantly differ (H=1) between mutation carriers and non-carriers
        protein_profiles{k}.proteins =  data_PRO_cleaned.Properties.RowNames(H_proteins==1);
        protein_profiles{k}.mutation_gene = data_MC3.Properties.RowNames(k);
    end
end
% 5.b. Analyze proteome-inferred cell types of mutation carriers vs non-carriers in key driver genes
%    for each key driver gene: define carriers and non-carriers, compare their
%    proteome-inferred cell enrichments then build cell-type profile for that mutation showing significance between carriers vs non-carriers
proteome_cell_type_profiles = cell(size(data_MC3,1),1);
for k = 1:size(data_MC3,1) %for each mutation gene
    % define carriers and non-carriers of the mutation in mutation_gene(k) in MC3 data
    carrier_true = table2array(data_MC3(k,:))~='wt';
    % map those carriers/non-carriers in proteome-inferred cell enrichments data
    [~, ip_carrier, ~] = intersect(data_PRO_cleaned_xCell_result.Properties.VariableNames, data_MC3.Properties.VariableNames(carrier_true));
    [~, ip_noncarrier, ~] = intersect(data_PRO_cleaned_xCell_result.Properties.VariableNames, data_MC3.Properties.VariableNames(~carrier_true));
    if min(length(ip_carrier),length(ip_noncarrier)) > 1
        % test the null hypothesis (H) that the proteome-inferred cell enrichments of the mutation carriers and non-carriers in mutation_gene(k) are not significantly different
        % H = 1 means that the null hypothesis can be rejected at 5% significance level
        H_cell_types = ttest2(table2array(data_PRO_cleaned_xCell_result(:,ip_carrier))', table2array(data_PRO_cleaned_xCell_result(:,ip_noncarrier))')';
        % any proteome-inferred cell-type whose enrichment scores significantly differ (H=1) between mutation carriers and non-carriers
        proteome_cell_type_profiles{k}.cell_types =  data_PRO_cleaned_xCell_result.Properties.RowNames(H_cell_types==1);
        proteome_cell_type_profiles{k}.mutation_gene = data_MC3.Properties.RowNames(k);
    end
end
% ------------------------------------------------------------------------
% Interpretation of analyses "5.a" and "5.b":
% ------------------------------------------------------------------------
k = 2;
% For the mutation in the "k-th key driver gene":
protein_profiles{k}.mutation_gene
%    the following list of "protein(s)" have different expression patterns between that mutation's carriers and non-carriers:
protein_profiles{k}.proteins
%    and the following inferred "cell-type(s)" are enriched differently as well:
proteome_cell_type_profiles{k}.cell_types
% That is, for k=2, the "ARID1A" gene's mutation-carriers and wild-types 
%    have different proteome signals on "374 proteins", which presumably
%    occurred in the cell-type "Myocytes".
%
% Similarly, for k=6, the "BRCA1" gene's mutation-carriers and wild-types 
%    have different proteome signals on "527 proteins", which presumably 
%    occurred in the cell-types "CD4+ Tem", "CMP", "MPP", "pro B-cells", "Th2 cells".
% This might suggest that these cell-types might be commonly exposed to mutations 
%    in the BRCA1 gene which may alter the production of those 527 proteins.
% ------------------------------------------------------------------------
% write profiles
fid = fopen('Supplementary_data.xls','w');
fprintf(fid,['Supplementary data. Each column represents the analysis of a ' ...
    'different gene mutation.\n The results include\n\t i) the group of cell types ' ...
    'that exhibit significanly different enrichment patterns between the mutation ' ...
    'carriers and non-carriers\n\t ii) the group of proteins that exhibit ' ...
    'significanly different expression patterns between the mutation carriers ' ...
    'and non-carriers \n\n\n']);
l_pctp_ct = zeros(length(find(signif_mut)),1);
l_pp_p = zeros(length(find(signif_mut)),1);
for k = find(signif_mut)
    fprintf(fid,'Mutation gene:\t');
    l_pctp_ct(k) = length(proteome_cell_type_profiles{k}.cell_types);
    l_pp_p(k) = length(protein_profiles{k}.proteins);
end
fprintf(fid,'\n');
for k = find(signif_mut)
    fprintf(fid,'%s\t',proteome_cell_type_profiles{k}.mutation_gene{:});
end
fprintf(fid,'\n\n');
for k = find(signif_mut)
    fprintf(fid,'Cell types:\t');
end
fprintf(fid,'\n');
for r = 1:max(l_pctp_ct)
    for k = find(signif_mut)
        if r <= l_pctp_ct(k)
            fprintf(fid,'%s\t',proteome_cell_type_profiles{k}.cell_types{r});
        else
            fprintf(fid,'\t');            
        end
    end
    fprintf(fid,'\n');
end    
fprintf(fid,'\n\n');
for k = find(signif_mut)
    fprintf(fid,'Proteins:\t');
end
fprintf(fid,'\n');
for r = 1:max(l_pp_p)
    for k = find(signif_mut)
        if r <= l_pp_p(k)
            fprintf(fid,'%s\t',protein_profiles{k}.proteins{r});
        else
            fprintf(fid,'\t');            
        end
    end
    fprintf(fid,'\n');
end    
fclose(fid);

%% 6. mRNA expression patterns due to mutation
% 6.a. Analyze mRNA-level gene expressions of mutation carriers vs non-carriers in key driver genes
%    for each key driver gene: define carriers and non-carriers, compare their mRNA expressions
%    then build mRNA profile for that mutation showing significance between carriers vs non-carriers
mRNA_profiles = cell(size(data_MC3,1),1);
for k = 1:size(data_MC3,1) %for each mutation gene
    % define carriers and non-carriers of the mutation in mutation_gene(k) in MC3 data
    carrier_true = table2array(data_MC3(k,:))~='wt';
    % map those carriers/non-carriers in mRNA data
    [~, im_carrier, ~] = intersect(data_mRNA_cleaned.Properties.VariableNames, data_MC3.Properties.VariableNames(carrier_true));
    [~, im_noncarrier, ~] = intersect(data_mRNA_cleaned.Properties.VariableNames, data_MC3.Properties.VariableNames(~carrier_true));
    if min(length(im_carrier),length(im_noncarrier)) > 1
        % test the null hypothesis (H) that the mRNA expressions of the mutation carriers and non-carriers in mutation_gene(k) are not significantly different
        % H = 1 means that the null hypothesis can be rejected at 5% significance level
        H_mRNAs = ttest2(table2array(data_mRNA_cleaned(:,im_carrier))', table2array(data_mRNA_cleaned(:,im_noncarrier))')';
        % any mRNA whose expression significantly differ (H=1) between mutation carriers and non-carriers
        mRNA_profiles{k}.mRNAs =  data_mRNA_cleaned.Properties.RowNames(H_mRNAs==1);
        mRNA_profiles{k}.mutation_gene = data_MC3.Properties.RowNames(k);
    end
end
% 6.b. Analyze mRNA-inferred cell types of mutation carriers vs non-carriers in key driver genes
%    for each key driver gene: define carriers and non-carriers, compare their
%    mRNA-inferred cell enrichments then build cell-type profile for that mutation showing significance between carriers vs non-carriers
mRNA_cell_type_profiles = cell(size(data_MC3,1),1);
for k = 1:size(data_MC3,1) %for each mutation gene
    % define carriers and non-carriers of the mutation in mutation_gene(k) in MC3 data
    carrier_true = table2array(data_MC3(k,:))~='wt';
    % map those carriers/non-carriers in mRNA-inferred cell enrichments data
    [~, im_carrier, ~] = intersect(data_mRNA_cleaned_xCell_result.Properties.VariableNames, data_MC3.Properties.VariableNames(carrier_true));
    [~, im_noncarrier, ~] = intersect(data_mRNA_cleaned_xCell_result.Properties.VariableNames, data_MC3.Properties.VariableNames(~carrier_true));
    if min(length(im_carrier),length(im_noncarrier)) > 1
        % test the null hypothesis (H) that the mRNA-inferred cell enrichments of the mutation carriers and non-carriers in mutation_gene(k) are not significantly different
        % H = 1 means that the null hypothesis can be rejected at 5% significance level
        H_cell_types = ttest2(table2array(data_mRNA_cleaned_xCell_result(:,im_carrier))', table2array(data_mRNA_cleaned_xCell_result(:,im_noncarrier))')';
        % any mRNA-inferred cell-type whose enrichment scores significantly differ (H=1) between mutation carriers and non-carriers
        mRNA_cell_type_profiles{k}.cell_types =  data_mRNA_cleaned_xCell_result.Properties.RowNames(H_cell_types==1);
        mRNA_cell_type_profiles{k}.mutation_gene = data_MC3.Properties.RowNames(k);
    end
end
% ------------------------------------------------------------------------
% Interpretation of analyses "6.a" and "6.b":
% ------------------------------------------------------------------------
% For k=6, the "BRCA1" gene's mutation-carriers and wild-types 
%    have different mRNA signals on "857 genes", which presumably 
%    occurred in the cell-types "CLP", "Macrophages", "Monocytes", "Neutrophils".
% This might suggest that these cell-types might be commonly exposed to mutations 
%    in the BRCA1 gene which may alter the mRNA production in those 857 genes.
% ------------------------------------------------------------------------
% Interpretation of analyses "5.a" and "6.a":
% ------------------------------------------------------------------------
% For k=6, the following "79 mRNA/protein" products are both altered 
%    due to mutation in BRCA1 gene (but not necessarily in the same cell-types):
intersect(mRNA_profiles{6}.mRNAs,protein_profiles{6}.proteins)
% ------------------------------------------------------------------------
 