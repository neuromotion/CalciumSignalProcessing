%This Code should be run after extracting calcium signals using the
%FluoroSNNAP software. The output of the FluoroSNNAP analysis will be the
%"processed_analysis.mat" file. If you renamed that file, then you should
%change the pfile name. 

%This code will output a Corr_SampleName.mat file which contains the
%correlation matrix (newcorr) of each sample within the processed analysis
%file.  The code will also output a Data_Output table, which contains the
%sample name, the average correlation value, and the whole-tissue firing
%rate.

%Written by Elaina Atherton, 2021. 

pfile = 'processed_analysis.mat'; %processed file name

load(pfile);
startd = 1;
endd = size(processed_analysis,2); %number of samples within the experiment
Ca_Output = zeros(endd, 2) ;
for i = startd:endd
    %extract sample name from the current row of processed_analysis
    filepath = [processed_analysis(i).filename];
    lastslash_pos = find(filepath == '/', 1, 'last');
    filelength = size(filepath,2);
    filename = filepath((lastslash_pos +1):filelength);
    Sample = extractBefore(filename,'.'); %extracted sample name
    
    %make correlation matrix from single cell activity
    newcorr = corr(processed_analysis(i).dF_cell'); 
    newcorr(all(isnan(newcorr),2),:) = [];    %remove rows that are all nan
    newcorr(:, all(isnan(newcorr),1)) = [];   %remove cols that are all nan
    
    %save correlation matrix file
    file = strcat('Corr_', Sample, '.mat');
    save(file,'newcorr');
   
    %save sample name
    Data_output(i).Sample = Sample;
    
    %calculate average correlation
    colvect= newcorr(find(~tril(ones(size(newcorr)))));
    avgcorr = nanmean(colvect);
    Data_output(i).AvgCorrelation = avgcorr;
    
    %extract Whole Tissue Firing rate
    firingrate = size(processed_analysis(i).Spikes_whole,2);
    Data_output(i).FiringRate = firingrate;
    
   
end
file2 = strcat('DataOutput_', pfile);
 save(file2,'Data_output');
