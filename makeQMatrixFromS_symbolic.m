function [Q2] = makeQMatrixFromS_symbolic(S,coefficients)
% This function inputs a Stoichiometric matrix (S) and a vector VmoveIndex
% of length ncols(S).  VmoveIndex is an index into the 'Moves' variable,
% and each move corresponds to a unique rate constant. 


nNodes = size(S,1);
nEdges = size(S,2);

% First make off-diagonal entries
Q2=sym(zeros(nNodes));
for x =1:nEdges
    Q2( S(:,x) == -1, S(:,x) == 1 ) = coefficients(x);
end

% Next make diagonal entries
for x =1:nNodes
    Q2(x,x) = Q2(x,x) - sum(Q2(x,:));
end



end
