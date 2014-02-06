function [Q] = makeQMatrixFromS(S,VmoveIndex)
% This function inputs a Stoichiometric matrix (S) and a vector VmoveIndex
% of length ncols(S).  VmoveIndex is an index into the 'Moves' variable,
% and each move corresponds to a unique rate constant. 

nNodes = size(S,1);
nEdges = size(S,2);


% First make off-diagonal entries
Q=zeros(nNodes);
for x =1:nEdges    % Loop through all COLUMNS of S
    Q( S(:,x) == -1, S(:,x) == 1 ) = VmoveIndex(x);
end

% Next make diagonal entries
for x =1:nNodes
    Q(x,x) = Q(x,x) - sum(Q(x,:));
end



end
