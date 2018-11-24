clear; close all; clc
%% 1. load data
load data_MC3.mat
load data_mRNA.mat
load data_PRO.mat
% common samples between the transcriptome and proteome data sets
[common_samples_mRNA_PRO, is_mRNA, is_PRO] = intersect(data_mRNA.Properties.VariableNames, data_PRO.Properties.VariableNames);
% common genetic elements between the transcriptome and proteome data sets
[common_genes_mRNA_PRO, ig_mRNA, ig_PRO] = intersect(data_mRNA.Properties.RowNames, data_PRO.Properties.RowNames);

%% 2. clean data
raw_values_mRNA = table2array(data_mRNA);
raw_values_PRO = table2array(data_PRO);
flag_normalize_mRNA = false;
flag_normalize_PRO = false;
% 2.a. impute data
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

% 2.b. remove outlier samples
c_dispersion = 3; %dispersion constant
non_outliers_mRNA = iqr(values_mRNA) <= mean(values_mRNA)+c_dispersion*std(values_mRNA);
values_mRNA = values_mRNA(:,non_outliers_mRNA);
non_outliers_PRO = iqr(values_PRO) <= mean(values_PRO)+c_dispersion*std(values_PRO);
values_PRO = values_PRO(:,non_outliers_PRO);

% 2.c. normalize samples
% mRNA | normalized across samples at 75%% quantile | log2 normalized
if flag_normalize_mRNA
    values_mRNA = values_mRNA-repmat(mean(values_mRNA),size(values_mRNA,1),1); % 0-mean
    values_mRNA = values_mRNA./repmat(std(values_mRNA),size(values_mRNA,1),1); % 1-SD
    log_mean_norm_error_mRNA = log(mean(abs(mean(values_mRNA)-0)))
    log_std_norm_error_mRNA = log(mean(abs(std(values_mRNA)-1)))
end
% PRO | samples normalized | 0-mean, 1-SD | 
if flag_normalize_PRO
    values_PRO = values_PRO-repmat(mean(values_PRO),size(values_PRO,1),1); % 0-mean
    values_PRO = values_PRO./repmat(std(values_PRO),size(values_PRO,1),1); % 1-SD
    log_mean_norm_error_PRO = log(mean(abs(mean(values_PRO)-0)))
    log_std_norm_error_PRO = log(mean(abs(std(values_PRO)-1)))
end
% check distributions
figure(1);
subplot(211); boxplot(values_mRNA); ylabel('BRCA Transcriptome');
subplot(212); boxplot(values_PRO); ylabel('BRCA Proteome');


%% 3. deconvolution methods
%
%% 3.a. tumor transcriptome
% xCell, performs cell type enrichment analysis from gene expression data for 64 
%   immune and stroma cell types by reducing associations between closely related cell types.
% xCell produces enrichment scores, not percentages.
% xCell uses the expression levels ranking and not the actual values, 
%   thus normalization does not have an effect.
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

%% 3.b. tumor proteome
% a quick observation shows that there are ~10.5K common genetic elements between the analyzed transcriptome and proteome data sets:
length(common_genes_mRNA_PRO)
% a deconvolution method designed for transcriptome data may worth studying
%   the proteome data that contains many common genetic signatures with transcriptome
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

%% 4. Expression patterns of mutation carriers and non-carriers
% 4.a. analyze proteomes of mutation carriers vs non-carriers in key driver genes
% for each key driver gene: define carriers and non-carriers, compare their proteomes
%   then generate protein profile for that mutation showing significance bw carriers vs non-carriers
protein_profiles = cell(size(data_MC3,1),1);
for k = 1:size(data_MC3,1) %for each mutation gene
    % define carriers and non-carriers of the mutation in mutation_gene(k) in MC3 data
    carrier_true = table2array(data_MC3(k,:))~='wt';
    if sum(carrier_true) > 1
        % map those carriers/non-carriers in proteome data
        [~, ip_carrier, ~] = intersect(data_PRO_cleaned.Properties.VariableNames, data_MC3.Properties.VariableNames(carrier_true));
        [~, ip_noncarrier, ~] = intersect(data_PRO_cleaned.Properties.VariableNames, data_MC3.Properties.VariableNames(~carrier_true));
        % test the null hypothesis (H) that the proteomes of the mutation carriers and non-carriers in mutation_gene(k) are not significantly different
        % H = 1 means that the null hypothesis can be rejected at 5% significance level
        H_proteins = ttest2(table2array(data_PRO_cleaned(:,ip_carrier))', table2array(data_PRO_cleaned(:,ip_noncarrier))')';
        % any protein whose expression significantly differ (H=1) between mutation carriers and non-carriers
        protein_profiles{k}.proteins =  data_PRO_cleaned.Properties.RowNames(H_proteins==1);
        protein_profiles{k}.mutation_gene = data_MC3.Properties.RowNames(k);
    end
end
% 4.b. analyze inferred cell types of mutation carriers vs non-carriers in key driver genes
% for each key driver gene: define carriers and non-carriers, compare their
% inferred cell enrichments then generate cell-type profile for that mutation showing significance bw carriers vs non-carriers
cell_type_profiles = cell(size(data_MC3,1),1);
for k = 1:size(data_MC3,1) %for each mutation gene
    % define carriers and non-carriers of the mutation in mutation_gene(k) in MC3 data
    carrier_true = table2array(data_MC3(k,:))~='wt';
    if sum(carrier_true) > 1
        % map those carriers/non-carriers in cell enrichments data
        [~, ip_carrier, ~] = intersect(data_PRO_cleaned_xCell_result.Properties.VariableNames, data_MC3.Properties.VariableNames(carrier_true));
        [~, ip_noncarrier, ~] = intersect(data_PRO_cleaned_xCell_result.Properties.VariableNames, data_MC3.Properties.VariableNames(~carrier_true));
        % test the null hypothesis (H) that the inferred cell enrichments of the mutation carriers and non-carriers in mutation_gene(k) are not significantly different
        % H = 1 means that the null hypothesis can be rejected at 5% significance level
        H_cell_types = ttest2(table2array(data_PRO_cleaned_xCell_result(:,ip_carrier))', table2array(data_PRO_cleaned_xCell_result(:,ip_noncarrier))')';
        % any cell-type whose enrichment scores significantly differ (H=1) between mutation carriers and non-carriers
        cell_type_profiles{k}.cell_types =  data_PRO_cleaned_xCell_result.Properties.RowNames(H_cell_types==1);
        cell_type_profiles{k}.mutation_gene = data_MC3.Properties.RowNames(k);
    end
end
% ------------------------------------------------------------------------
% Interpretation of analyses "4.a" and "4.b":
% ------------------------------------------------------------------------
k = 2;
% For the mutation in the "k-th key driver gene":
protein_profiles{k}.mutation_gene
% the following list of "protein(s)" have different expression patterns between that mutation's carriers and non-carriers:
protein_profiles{k}.proteins
% and the following inferred "cell-type(s)" are enriched differently as well:
cell_type_profiles{k}.cell_types
%
% That is, for k=2, the "ARID1A" gene's mutation-carriers and wild-types 
% have different proteome signals on "374 proteins", which presumably
% occured in the cell-type "Myocytes".
%
% For k=6, the "BRCA1" gene's mutation-carriers and wild-types 
% have different proteome signals on "527 proteins", which presumably 
% occured in the cell-types "CD4+ Tem", "CMP", "MPP", "pro B-cells", "Th2 cells".
% ------------------------------------------------------------------------
 
