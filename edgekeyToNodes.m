function nodes=edgekeyToNodes(key)
    % This function returns integers, given a Edge-hash key
    % These elements are the from, to, and ID of the Edge.
    
    
    nodesCell = textscan(key,'%u','Delimiter',',');  % Splits the key into integers, comma-separated
    nodes = nodesCell{1};
end
