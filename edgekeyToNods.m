function Moves=keyToMoves(key)
    % This function returns a list of moves given a hash key
    
    movesCell = textscan(key,'%u','Delimiter',',');  % Splits the key into integers, comma-separated
    Moves = movesCell{1};
end
