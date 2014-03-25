function exons=nodeNameToIntrons(NodeName,spliceGraph)
exons = movesToIntrons(move_IdToStruct(keyToMoves(NodeName),spliceGraph),spliceGraph);

end