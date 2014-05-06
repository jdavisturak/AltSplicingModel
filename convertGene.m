%% Function to convert UCSC-style data into exon or intron structure
function out=convertGene(row)

out = struct;
% Here convert raw UCSC row into numbered structure.
% Convert to num; row = [txStart,txEnd,exonStarts,exonStops]

exonStarts = str2num(row{3});  % converts comma-separated string into number array
exonEnds   = str2num(row{4});
strand = row{5};

out.exons = exonEnds - exonStarts;
out.introns = exonStarts(2:end) - exonEnds(1:(end-1));

if(strand=='-')
    out.exons = fliplr(out.exons);
    out.introns = fliplr(out.introns);
end

out.NumIntrons = length(out.introns);

end