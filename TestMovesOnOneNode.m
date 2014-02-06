function TestMovesOnOneNode(ThisNodeName,NewMoves)
global Nodes NumNodes AllMoves Edges;
ThisNodeID = Nodes.get(ThisNodeName);

% Get previous moves from this node 
previousMoves = AllMoves(keyToMoves(ThisNodeName));

%% Loop through each possible new move
% NewMoves, Nodes, NodeString

for (moveIndex = 1:length(NewMoves))
    thisMove = NewMoves(moveIndex);
    
    % Is this move legal?
    if isMoveAllowed(thisMove,previousMoves)
        % If allowed, construct the string of a new node, otherwise skip
        targetNodeName = addToNodeName(ThisNodeName,thisMove);
        targetNode = Nodes.get(targetNodeName);        % Retrieve this node from the hashtable
        
        % If this node does exist, simply add an edge
        % Otherwise, A) Add the node; B) Add the edge C)Try ALL moves with this
        % node.
        
        % Test if that node exists already
        if isempty(targetNode)
            
            % Add this node to the hashtable
            NumNodes = NumNodes+1;
            targetNode = NumNodes;
            Nodes.put(targetNodeName,targetNode);
          
            % Add EDGE to the record of this node
            addEdge(thisMove,ThisNodeID,targetNode); % Reports that this move was used between these two nodes
            
            %fprintf('Added new node %d\n',targetNode);
            %disp(Edges);                        
            
            % test this NEW node, with all OTHER moves!
            otherMovesIndex = ~ismember([AllMoves.index],keyToMoves(targetNodeName));             
            TestMovesOnOneNode(targetNodeName,AllMoves(otherMovesIndex));
            
            
        else
            % Add EDGE
            addEdge(thisMove,ThisNodeID,targetNode); % Reports that this move was used between these two nodes
           
            %fprintf('Added new edge only %d\n',thisMove.index);
            %disp(Edges);
        end        
    end   
end