function key=movesToKey(movesArray)
    % This function returns a hash key, given a list of moves as input
    
    % right now, moves should be an array of doubles
    key = sprintf('%d,',sort(movesArray));
    key = key(1:(end-1));    
end
