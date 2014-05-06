function [ Mean ] = getPTS_time( SG, KKs, startingConcentration,Q )
%getPTS_time employ the phase-type distribution to calculate the expected
%value of comign to a fully splicec variant
%   SG is a spliceGraph and KKs is a vector of splicing rate constants (non-symbolic)
% See http://en.wikipedia.org/wiki/Phase-type_distribution, http://www.cs.cmu.edu/~osogami/thesis/html/node39.html

if nargin < 4
    % Obtain infinitesimal matrix of my markov chain
    Q = makeQMatrixFromS(SG.S,KKs(SG.VmoveIndex));
end

if nargin < 3
    startingConcentration = [1 zeros(1,SG.NumNodes-1)];
end


% make the subgenerator matrix S for the transient states
nonLeaf = unique(SG.VnodeIndex);
S=Q(nonLeaf,nonLeaf);
a = startingConcentration(nonLeaf);

% % Now make the Q matrix where node 0 is the absorbing node
% Q_ph = [zeros(1,length(nonLeaf)+1);-S*ones(length(nonLeaf),1) S];

% Obtain the mean:
Mean = -a * inv(S) * ones(size(S,2),1);

end

