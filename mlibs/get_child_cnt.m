function [child_cnt] = get_child_cnt(pred)
%% Get the number of child of each node in the MST from the predecessor 
%  vector.
%
% Args:
%   pred: the predecessor vector of each node. NaN means non-connected.
%         the last node represents the sink. pred(-1) = 0 (root).
%
% Return:
%   child_cnt: the number of child at each node. 0 means non-connected or
%              leaf node

% init result vector to zero
n = length(pred);         % include the last sink node
child_cnt = zeros(1, n);

for i = 1:n
    if ~isnan(pred(i))    % if this node is connected
        % add one to all the predecessor, all the way to the sink
        cur_pred = pred(i);
        while cur_pred ~= 0
            child_cnt(cur_pred) = child_cnt(cur_pred) + 1;
            cur_pred = pred(cur_pred); % iterate
        end
    end
end
%fprintf('get child cnt:\n');
%disp(pred);
%disp(child_cnt);
end

