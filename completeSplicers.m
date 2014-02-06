function leafNodes = completeSplicers(spliceGraph)
    leafNodes = setdiff(1:size(spliceGraph.S,1),incompleteSplicers(spliceGraph));
end