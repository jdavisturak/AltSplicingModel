function exons=multiple_nodesToExons(nodesCell,spliceGraph,makeChar)
if ~exist('makeChar','var')
   makeChar = false; % Get default option
end

% Convert multiple chars in nodesCell (each is comma-separated list of moves) to cell of exons
    %(each is a double or comma-separated list of exons included in that species;
    
exons = cellfun(@(x)nodeNameToExons(x,spliceGraph),nodesCell,'UniformOutput',false);

if makeChar 
    exons = cellfun(@(x)sprintf('%02de,',x),exons,'UniformOutput',false);
end

end