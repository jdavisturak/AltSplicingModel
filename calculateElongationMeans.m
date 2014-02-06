%% From Lev's analytical solution to the original model
% This function calculates the mean elongation times for a given gene's regions.
% Jeremy Davis-Turak 6/27/12

function total_elongation=calculateElongationMeans(gene,k_e_vec,k_p_vec,k_u_vec)

%% Calculate Time for each region
% Separate the regions of the gene into N1, N2, ... The beginning/end of a region is the end of the intron
regionStarts = gene.exonStarts(2:end);
regionEnds   = [regionStarts(2:end)-1 gene.exonEnds(end)];
numRegions = length(regionStarts);

total_elongation = zeros(1,numRegions);

% Loop through each region
for(i = 1:(numRegions))
    these_bp = regionStarts(i):regionEnds(i);
    total_elongation(i)=expected_elongation(k_e_vec(these_bp),k_p_vec(these_bp),k_u_vec(these_bp));  % Get expected value of elongation
end

   function total_elongation=expected_elongation(k_e_vec,k_p_vec,k_u_vec)
        total_elongation = sum((1+k_p_vec./k_u_vec)./(k_p_vec+k_e_vec));
    end
end