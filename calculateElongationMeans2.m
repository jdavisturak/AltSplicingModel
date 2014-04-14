%% This function calculates the mean elongation times for a given gene's regions.
% Jeremy Davis-Turak 5/16/13: Using simply NET elongation rate

function total_elongation=calculateElongationMeans2(gene,NetElong)

% Regions start at the first intron
total_elongation = ([gene.introns(2:end) 0] + gene.exons(2:end)) ./ NetElong;


end