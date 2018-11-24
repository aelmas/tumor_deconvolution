
%% 1. load data
load data_MC3.mat
load data_mRNA.mat
load data_PRO.mat

%% 2. clean data
raw_values_mRNA = table2array(data_mRNA(2:end,2:end));
raw_values_PRO = table2array(data_PRO(2:end,2:end));
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

%% 3. deconvolution methods
%% 3.a. tumor transcriptome
% xCell, performs cell type enrichment analysis from gene expression data for 64 
% immune and stroma cell types by reducing associations between closely related cell types
% installation\\$ devtools::install_github('dviraran/xCell')
%
% pipe data to xCell
data_mRNA_cleaned = data_mRNA;
data_mRNA_cleaned(2:end,2:end) = array2table(values_mRNA);
writetable(data_mRNA_cleaned,'BRCA_mRNA_formatted_normalized_cleaned.txt','Delimiter','\t');
%
% script to call xCell 
fid = fopen('call_xCell.R','w');
fprintf(fid,[
    'library(xCell)\n' ...
    'expression_matrix = read.table("BRCA_mRNA_formatted_normalized_cleaned.txt",header=TRUE,row.names=1, as.is=TRUE)\n' ...
    'xCell_result = xCellAnalysis(expression_matrix)\n' ...
    'write.table(xCell_result,file = "xCell_result_BRCA_mRNA_formatted_normalized_cleaned.txt",sep="\t",quote=FALSE)\n' ...
    'return(xCell_result)' ...
    ]);
fclose(fid); 
% grand user permission of file executions
system('chmod u+x call_xCell.R xCell-master/R/xCell.R') 
% call xCell
system('/usr/local/bin/Rscript call_xCell.R');



