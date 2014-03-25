function introns=multiple_nodesToIntrons(nodesCell,spliceGraph,makeChar)
if ~exist('makeChar','var')
   makeChar = false; % Get default option
end

% Convert multiple chars in nodesCell (each is comma-separated list of moves) to cell of introns
    %(each is a double or comma-separated list of introns included in that species;
    
introns = cellfun(@(x)nodeNameToIntrons(x,spliceGraph),nodesCell,'UniformOutput',false);

if makeChar 
    for i = 1:length(introns)
        if ~isempty(introns{i})
            introns{i} = sprintf('%02di,',introns{i});
        else
            introns{i} = '';
        end
    end
end

end