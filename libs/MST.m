function [G, pred] = MST(A, c, R)
%% Get the minimal spanning tree from the locations of the sensors and 
%  the sink.
%
% Args:
%   A: list of locations of the sensors
%   c: location of the sink
%   R: communication range
%
% Return:
%   G: matrix showing the MST
%   pred: predecessor array
addpath('./lldistkm/');

nodes = vertcat(A, c);  % add the sink to the nodes array
n = size(nodes, 1);     % number of nodes
G = zeros(n);           % init a symmetric matrix for connection graph
                        
% create the undirected graph and fill the matrix
for p = 1:n
    for q = p+1:n
        % check the distance between nodes
        [d1km, d2km] = lldistkm(nodes(p, :), nodes(q, :));
        if d1km < R
            % update the undirected graph
            G(p, q) = d1km;
            G(q, p) = d1km;
        end
    end
end

% find the minimal spanning tree in graph, with the sink as the root
[Tree, pred] = graphminspantree(sparse(G), n);
end

