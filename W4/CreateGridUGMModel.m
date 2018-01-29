function [edgePot,edgeStruct]=CreateGridUGMModel(nRows, nCols, nStates, lambda)
%
%
% NumFils, NumCols: image dimension
% K: number of states
% lambda: smoothing factor

tic

nNodes = nRows * nCols;

adj = sparse(nNodes, nNodes);

% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;
 
% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;
 
% Add Up/Left Edges
adj = adj+adj';

% Edge struct
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

% Edge pot
potts_val = exp(lambda);
potts_mtrx = potts_val(2)*eye(nStates) + potts_val(1)*(ones(nStates) - eye(nStates));

edgePot = zeros(nStates, nStates, edgeStruct.nEdges);
for e = 1:edgeStruct.nEdges
    edgePot(:,:,e) = potts_mtrx;
end


toc;