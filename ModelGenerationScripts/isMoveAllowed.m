function allowed = isMoveAllowed(newMove,previousMoves)
%% Function to test whether a new move is allowed for a node with moves 'previousMoves'

    if size(previousMoves,1) == 0
        allowed = 1;
        return;
    end

    froms = [previousMoves.from];
    tos   = [previousMoves.to];
    
    % Rules: NOT allowed if any, for all i(prev moves): 
    %           1) the new 'from' is between any to/from pair
    %           2) the new 'to' is between any to/from pair
    %           3) any 'to' is between the new to/from pair

    test1 = newMove.from <= tos & newMove.from >= froms;
    test2 = newMove.to <= tos & newMove.to >= froms;
    test3 = newMove.from <= tos & newMove.to >= tos;
    
    allowed = ~( any(test1)  | any(test2) | any(test3));
end