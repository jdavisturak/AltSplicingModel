
% cd C:\Users\jeremy\Documents\Splicing\ExonsAnalysis

fid = fopen('hg19_refseq_072712.bed');
genesFile = textscan(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s');
fclose(fid);
genesData = [genesFile{5} genesFile{6} genesFile{10} genesFile{11} genesFile{4} genesFile{2} genesFile{13}];
genesData = genesData(2:end,:);

% Parameters:
NetElong = 2000; % kb/min
Ks       = 2;

% Loop through each gene 
myData = cell(size(genesData,1),4);

parfor g = 1:size(genesData,1)
% for g = 1:100
    gene=convertGene(struct('row',{genesData(g,:)}),true,1,'elong_splice1');
%    if(length(gene.introns)==8)
%        disp(g)
%    end

    try
        tic
        [FinalConcentration,OutputTimes,OutputX , FinalConcentration2,OutputTimes2,OutputX2, TranscriptNames,whichCompleteSplicers ] = ...
            runFullModel3_max9(gene,NetElong,repmat(Ks,1,length(gene.introns)*(length(gene.introns)+1)/2));
    
        myData(g,:) ={genesData{g,6}, genesData{g,7},FinalConcentration2(whichCompleteSplicers), TranscriptNames(whichCompleteSplicers)};
        toc
        
    catch err
    
    end
    
end