function newNodeName = addToNodeName(ThisKey,thisMove)
    if(ThisKey(1) == 't')
        newNodeName = num2str(thisMove.index);
    else
        oldMoves = keyToMoves(ThisKey);
%         thisMove.index
%         oldMoves
        newNodeName = movesToKey([oldMoves; thisMove.index]);
    end

end