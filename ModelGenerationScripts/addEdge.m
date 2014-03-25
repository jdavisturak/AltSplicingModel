function addEdge(thisMove,NodeFrom,NodeTo)
global Edges NumEdges;

% Edges hash: 
%   Key=> Tuple of nodes (IN/OUT)
%   Value=> 'pointer' to the Move that it inherits

NumEdges = NumEdges + 1;

EdgeKey = sprintf('%d,%d,%d',NodeFrom,NodeTo,NumEdges);

Edges.put(EdgeKey,thisMove.index);

end