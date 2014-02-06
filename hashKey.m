function key=hashKey(Moves)
    % This function returns a unique name (key) of the node with moves Moves.
    
    % Moves is a list of all moves taken to get to this node.  The order is not important here
    Moves = sort(Moves);
    
    key = sprintf('%d,',Moves);
    key = key(1:(end-1)); % Get rid of trailing ','
end
