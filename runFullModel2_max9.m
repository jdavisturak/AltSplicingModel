%% JCDT 6/27/12
% This is the first attempt to use the alternative splicing model paired with calculating the mean
% elongation times.
% This version (currently) uses the pre-calculated S matrix of genes up to 9 introns long.
%clear; runVersion = 'Gillespie3.3'; LTa = {'31648071',	'31650077' ,'31648071,31648489,31648683,31649036',	'31648202,31648597,31648789,31650077','+'}; gene=convertGene(struct('row',{LTa}),true,1,'elong_splice1');

function [FinalConcentration,OutputTimes,OutputX , FinalConcentration2,OutputTimes2,OutputX2, TranscriptNames,whichCompleteSplicers ]=runFullModel2_max9(gene,NetElong,k_s_vec,options,verbose)
MAX_INTRONS = 9;

global GLOBAL_spliceGraph;
if isempty(GLOBAL_spliceGraph)    
    load('AltSplicing_1-9_structures.mat')
    GLOBAL_spliceGraph = spliceGraph;
end

if (nargin < 7)
    verbose = true;
end

% if(~isfield(options,'maxTranscripts'))
%     maxTranscripts = 10^10; % Default max.  Set to 33 introns
% else
%     maxTranscripts = options.maxTranscripts;
% end

NumIntrons=length(gene.introns);
if NumIntrons > MAX_INTRONS
    error('Too many introns! Max allowed is %d introns\n',MAX_INTRONS);
end


TranscriptNames = GLOBAL_spliceGraph(NumIntrons).TranscriptNames;
whichCompleteSplicers = completeSplicers(GLOBAL_spliceGraph(NumIntrons));

%% Calculate Time for each region
% region_elongation_times = calculateElongationMeans(gene,k_e_vec,k_p_vec,k_u_vec);
region_elongation_times = calculateElongationMeans2(gene,NetElong);



%% SPLICING:
% make sure that the K_s is of the correct dimensions
k_s_vec = reshape(k_s_vec,numel(k_s_vec),1);

% 1) Co-transcriptional splicing
[FinalConcentration,OutputTimes,OutputX] = SpliceDuringTxn(NumIntrons,region_elongation_times,k_s_vec,GLOBAL_spliceGraph);


% 2) Post-transcriptional splicing
[FinalConcentration2,OutputTimes2,OutputX2] = SpliceAllPostTxn(FinalConcentration,GLOBAL_spliceGraph(NumIntrons),region_elongation_times(end),k_s_vec);



end