%% This function calculates the mean elongation times for a given gene's regions.
% Jeremy Davis-Turak 5/16/13: Using simply NET elongation rate

function total_elongation=calculateElongationMeans2(gene,NetElong)

%% Calculate Time for each region
% Separate the regions of the gene into N1, N2, ... The beginning/end of a region is the end of the intron
regionStarts = gene.exonStarts(2:end);
regionEnds   = [regionStarts(2:end)-1 gene.exonEnds(end)];
numRegions = length(regionStarts);

total_elongation = zeros(1,numRegions);

% Loop through each region
for (i = 1:(numRegions))
    total_elongation(i) = abs(regionStarts(i) - regionEnds(i) + 1) / NetElong;
    
end

end