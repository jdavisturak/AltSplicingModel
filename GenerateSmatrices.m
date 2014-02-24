%% This script creates the S matrix and the data structures describing the splicing graph

%% IMPORTANT STRUCTURES TO KEEP TRACK OF

% Possible Moves: 1*M struct

% Nodes:  Hashtable, where keys uniquely identify the node; value is index?

% Edges hash
%   Key=> Tuple of nodes (IN/OUT) AND EdgeID
%   Value=> 'pointer' to the Move that it inherits

clear Nodes NumNodes AllMoves Edges NumEdges; clear global
global Nodes NumNodes AllMoves Edges NumEdges;

%% Nodes Hash
Nodes = java.util.Hashtable;
Nodes.put('t',1); % First node: pre-mRNA
NumNodes = 1;

%% Edges 
% For each node,keep track of WHICH MOVE it is; also, WHICH NODES use this move?
Edges = java.util.Hashtable;
NumEdges = 0;

%% Main Program:
% Loop through all N introns.
AllMoves = struct([]); 
for N = 1:3 % number of introns to consider
   
    %fprintf('Beginning Test of intron %d\n',N)
    %%  Get new moves that this intron introduces
    
    %NewMoves = addNewMoves(N,size(AllMoves,2));
    NewRecruitMoves = addNewRecruitMoves(N,size(AllMoves,2));
    AllMoves = [AllMoves NewRecruitMoves];
    
    NewSpliceMoves = addNewSpliceMoves(N,size(AllMoves,2));
    AllMoves = [AllMoves NewSpliceMoves];
    
    
    
    %% For each EXISTING node, add these moves, if possible
    
    % Get names from hash
    myNodesTempNames = cell(NumNodes,1);
    enum = Nodes.keys();
    count = 0;
    while enum.hasMoreElements   
        count = count+1;
        myNodesTempNames{count} = enum.nextElement;
    end
    
    % Next, loop through all EXISTING nodes and add their legal moves
    % This ALSO adds new nodes
    for(i = 1:length(myNodesTempNames))
        ThisNodeName = myNodesTempNames{i};
        %ThisNodeID = Nodes.get(ThisNodeName);
        %fprintf('Testing nodes afresh\n')
        TestMovesOnOneNode(ThisNodeName,[ NewRecruitMoves NewSpliceMoves])
    end
    
    %% Make S matrix: sparse Matrix notation
    % in S, rows are reactants (nodes), columns are reactions/fluxes (edges)
    S = sparse(NumNodes,NumEdges);
    % Also, Make VmoveIndex.  This maps individual nodes into the moves they belong to
    VmoveIndex = zeros(NumEdges,1);
    
    % loop through fluxes
    enum = Edges.keys();
    while enum.hasMoreElements   % Loop through all EXISTING nodes
        key = enum.nextElement;
        [nodes] = edgekeyToNodes(key);
        move = Edges.get(key);
        S(nodes(1),nodes(3)) = -1; %from
        S(nodes(2),nodes(3)) = 1;
        VmoveIndex(nodes(3)) = move;
    end
    
    %% Need to make a way of easily getting the fluxes during the simulation.  
    % First thing is to get the node's concentration: that comes from the -1 entry.
    % (wrong)Second thing is to divide by the number of edges leaving that node.
    
    S2 = S;
    S2(S2==1) = 0;
    % 1)
    VnodeIndex = -1*(1:NumNodes)*S2;  % Vector of size NumEdges; this is an index that tells me which node is supplying the concentration for each flux
%     % 2)
%     OutPerNode = full(sum(S2,2));
%     OutPerNode(OutPerNode==0) = Inf;
%     VnodeMultiplier = -1./OutPerNode; % Vector of size NumNodes; in the simulation, I will need to multiply this times the vector of node 'conentrations' to determine to flux leaving each node
%    
%     % In the simulation, I will use this Code to determine each flux:
%     % FluxByNode = (Node_Concentraions.*VnodeMultiplier); 
    % Flux = FluxByNode(VnodeIndex);
        
    save(sprintf('AltSpliceRecruitMatrices%2d.mat',N), 'Nodes', 'NumNodes', 'AllMoves', 'Edges', 'NumEdges','S','VnodeIndex','VmoveIndex');

end

%% Save data together
clear spliceGraph;
for  N =1:3
    disp(N)
    load(sprintf('AltSpliceRecruitMatrices%2d.mat',N))

    % Reverse the hash: map node ID's to node keys:
    NodesTable = cell(NumNodes,1);
    enum = Nodes.keys();
    count = 0;
    while enum.hasMoreElements   % Loop through all EXISTING nodes
        key = enum.nextElement;
        NodesTable{Nodes.get(key)} = key;
        count = count + 1;
    end
    
    % Final Save
    spliceGraph(N) = struct('S',S,'VnodeIndex',VnodeIndex,'VmoveIndex',VmoveIndex','Nodes',Nodes,'NodesTable',{NodesTable},'Edges',Edges,'Moves',AllMoves,'NumIntrons',N,'NumNodes',NumNodes,'NumEdges',NumEdges,'TranscriptNames','');
    spliceGraph(N).TranscriptNames = multiple_nodesToTranscript(NodesTable,spliceGraph(N));
      
end
save(sprintf('AltSplicingRecruit_1-%d_structures.mat',N),'spliceGraph');


%% Report number of edges:
%[AllMoves.from; AllMoves.to; AllMoves.index]
enum = Nodes.keys();
count = 0;
while enum.hasMoreElements   % Loop through all EXISTING nodes
    enum.nextElement;
    count = count + 1;
end
enum = Edges.keys();
count2 = 0;
dat=[];
while enum.hasMoreElements   % Loop through all EXISTING nodes
    edge = enum.nextElement;
    count2 = count2 + 1;
    dat(end+1) = Edges.get(edge);
end
fprintf('%d Nodes, %d Edges\n',count,count2);
