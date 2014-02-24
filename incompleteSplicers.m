function nonLeafNodes=incompleteSplicers(spliceGraph)
nonLeafNodes = unique(spliceGraph.VnodeIndex);
% VnodeIndex is Vector of size NumEdges:
%   this is an index that tells me which node is supplying the concentration for each flux
%   therefore, any node in VnodeIndex cannot be a leaf node.
end