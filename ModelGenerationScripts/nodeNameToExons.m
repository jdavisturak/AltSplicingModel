function exons=nodeNameToExons(NodeName,spliceGraph)
exons = movesToExons(move_IdToStruct(keyToMoves(NodeName),spliceGraph),spliceGraph,1);

end