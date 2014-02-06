function introns=movesToIntrons(Moves,spliceGraph)
% This function takes a set of Moves and returns a vector of introns that would be expressed in
% such a transcript.

% Most generic solution, hopefully fastest!
% Initially set ALL introns to included.
introns = 1:(spliceGraph.NumIntrons);

% Then, find those splices that remove introns
skipped = unique(cell2mat(arrayfun(@(x) x.from:x.to ,Moves,'UniformOutput',false)));

introns = setdiff(introns,skipped);


end

