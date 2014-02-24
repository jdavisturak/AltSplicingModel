function moves=addNewRecruitMoves(n,Nmoves)
    
    % n = the current intron
    
    % Nmoves = the # of total moves already in existence
    
    % The 'from' is the current intron: the 'to' is left as 0, as a flag
    
    % 'recruit' tells us which spliceosome is being recruited
    
    moves = struct('from',{n,n},'to',{0,0},'recruit',{'5','3'},'index',num2cell(Nmoves+(1:2)));
end