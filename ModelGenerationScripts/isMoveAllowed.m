function allowed = isMoveAllowed(newMove,previousMoves)
%% Function to test whether a new move is allowed for a node with moves 'previousMoves'


    if size(previousMoves,1) == 0
        allowed = 1;
        return;
    end

    % There are two sets of rules: one for splicing reactions, one for recruiting reactions
    
    froms = [previousMoves.from];
    tos   = [previousMoves.to];
    recruits   = [previousMoves.recruit];
    
        
     
    if newMove.to == 0
        %% Recruitment move: 
        % Rules: allowed if the intron still exists - in other words, NOT allowed for any i(prev moves):
        %      1) the new 'from' is between any to/from any Splicing pair
        %      2) This same Move has not alread been done
        
        test1 = newMove.from <= tos & newMove.from >= froms & recruits=='-';
        test2 = newMove.from == froms & recruits == newMove.recruit;
        allowed = ~( any(test1)  | any(test2));
        
        
    else
        %% Splicing move
        
        % Rules: NOT allowed if any, for any i(prev moves):
        %           1) the new 'from' is between any to/from Splicing pair
        %           2) the new 'to' is between any to/from pair
        %           3) any 'to' is between the new to/from pair
        %        Also, MUST have the following conditions apply, in any i (prev moves):
        %           4) previous Recruit move '5' with same 'from' must exist
        %           4) previous Recruit move '3' with to current 'to' must exist
        
        
        test1 = newMove.from <= tos & newMove.from >= froms & recruits=='-';
        test2 = newMove.to <= tos & newMove.to >= froms & recruits=='-';
        test3 = newMove.from <= tos & newMove.to >= tos & recruits=='-';
        test4 = newMove.from == froms & recruits=='5';
        test5 = newMove.to == froms & recruits=='3';
        
        allowed = any(test4) & any(test5) & ~( any(test1)  | any(test2) | any(test3));        
        %tests=[test1; test2; test3;test4; test5]
    end
end




