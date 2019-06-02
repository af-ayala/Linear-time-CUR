%% Select target points using Nearest-Neighbors
% Require:
% Y: set of n target points, Y is an mxd matrix, d: geometric dimension
% t: number of selected points Y
% Returns:
% J: indices of target points closest to the source domain
function [P] = NN_Sampling(Y,X,t)

% Finding the distance
    DX = bsxfun(@minus,Y(:,1),X(:,1)');
    DY = bsxfun(@minus,Y(:,2),X(:,2)');
    DZ = bsxfun(@minus,Y(:,3),X(:,3)');
    D = sqrt(DX.^2+DY.^2+DZ.^2); % The i-th line of D is the distance from
    % the i-th target point to X
    d = min(D(:));

for i=1:size(D,1)
    [dist(i),~] = min(D(i,:));
end
[~,P] = mink(dist,t); % Find t points on Y closest to X
end
