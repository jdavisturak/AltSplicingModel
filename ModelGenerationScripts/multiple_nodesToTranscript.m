function [transcripts, exons, introns] = multiple_nodesToTranscript(nodesCell,spliceGraph)

% Convert multiple chars in nodesCell (each is comma-separated list of moves) to cell of exons
    %(each is a double or comma-separated list of exons included in that species;
    
exons = multiple_nodesToExons(nodesCell,spliceGraph,1);
introns = multiple_nodesToIntrons(nodesCell,spliceGraph,1);

transcripts = {};
for i = 1:length(exons)
    T = strcat(exons{i},introns{i});
    tCell = textscan(T,'%s','Delimiter',',');  % Splits the key into integers, comma-separated
    elements = sort(tCell{1});
    J = sprintf('%s,',elements{:});
    transcripts{i} = J(1:(end-1));
end