function [Q] = makeQMatrixFromS_noDiags(S,VmoveIndex)
% This function inputs a Stoichiometric matrix (S) and a vector VmoveIndex
% of length ncols(S).  VmoveIndex is an index into the 'Moves' variable,
% and each move corresponds to a unique rate constant. 

nNodes = size(S,1);

% First make off-diagonal entries
Q=zeros(nNodes);
for x =1:nNodes
    Q( S(:,x) == -1, S(:,x) == 1 ) = VmoveIndex(x);
end

end
