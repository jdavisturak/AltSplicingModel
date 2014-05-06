function [ Knew ] = splicingRateLengthFunction( GENE, K0,SG,exponent, minimum)
%splicingRateLengthFunction Returns the K_splice vector for all moves based
%on a 1/D penalty for distance between splice sites
%  Distances less than 'minimum' get no penalty.
%  SG is the spliceGraph struct for the gene of interest
% Distance = sum over (i='from' ... 'to') of intron(i) + sum over (i='from+1' ... 'to') of exon (i)

if nargin < 5
    minimum = 170;
end

% Loop through 'Moves'
Distance = arrayfun(@(move) sum(GENE.introns(move.from:move.to)) + sum(GENE.exons((move.from+1):move.to)),SG.Moves);
Knew = K0.* (minimum ./max([Distance;repmat(minimum,1,length(SG.Moves))])).^exponent;

end

