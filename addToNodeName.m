function newNodeName = addToNodeName(ThisKey,thisMove)
% Get the name of the next node give old node 'ThisKey'
    if(ThisKey(1) == 't')
        newNodeName = num2str(thisMove.index);
    else
        oldMoves = keyToMoves(ThisKey);
%         thisMove.index
%         oldMoves
        newNodeName = movesToKey([oldMoves; thisMove.index]);
    end

end