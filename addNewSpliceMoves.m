function moves=addNewSpliceMoves(n,Nmoves)
    % n = the newest intron to add
    
    % Nmoves = the # of total moves already in existence
    
    % The'from' and 'to' fields indicate which introns' splice sites are
    % joined in this splicing reaction
    
    moves = struct('from',num2cell(1:n),'to',n,'recruit','-','index',num2cell(Nmoves+(1:n)));
end