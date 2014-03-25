function exons=movesToExons(Moves,spliceGraph,includeOthers)
    % This function takes a set of Moves and returns a vector of exons that would be expressed in
    % such a transcript.
    
    %% First convert from intron splice sites  to exons included.
    % For example, 1-1 means INTRON 1 5' SS to INTRON 1 3' SS
    %        1-1 means INTRON 1 5' SS to INTRON 2 3' SS, which SKIPS exon 2!!!
    % Therefore the first intron tells us that exon 'FROM' is included.
    %           the second intron tells that exon 'TO' + 1 is included
    
    starts = [Moves.from]; % These are now EXON coordinates
    ends = [Moves.to] + 1;
    
    
    %% But what the exons before and after these guys?  And in the middle?
    if exist('includeOthers','var')
        if includeOthers          
            % Most generic solution, hopefully fastest!
            % Initially set ALL exons to included.            
            exons = 1:(spliceGraph.NumIntrons+1);
            
            % Then, find those splices that skip exons.
                % skips are those where the ending EXON is > 1 from the starting exon.
            skipped = unique(cell2mat(arrayfun(@(x) (x.from+1):x.to ,Moves,'UniformOutput',false)));
            
            exons = setdiff(exons,skipped);            
            return
            
        end
    end
    
    %% Default: simple way
    
    exons = unique([starts; ends]);  % This initally get only the ones that have been spliced together: doesn't include the untranscribed ones
    
    
end

