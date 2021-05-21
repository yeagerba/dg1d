function [B] = spblkdiag(A)
% Determine the sizes of the blocks (assumed to be the same size)
[n,m] = size([A{1}]);
% Determine the total number of rows and columns
N = n*numel(A); M = m*numel(A);
% Create the sparse block diagonal matrix
i = reshape(repmat(reshape((1:N)',n,[]),m,1),[],1);
j = reshape(repmat(1:M,[n,1]),[],1);
s = reshape([A{:}],[],1);
B = sparse(i,j,s);
end

