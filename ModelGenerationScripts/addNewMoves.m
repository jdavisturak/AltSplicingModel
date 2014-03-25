function moves=addNewMoves(n,Nmoves)
    % n = the newest intron to add
    % Nmoves = the # of total moves already in existence
    moves = struct('from',num2cell(1:n),'to',n,'index',num2cell(Nmoves+(1:n)));
end