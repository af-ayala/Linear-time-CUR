%% Select target points using geometrically balanced partition
% Require:
% Y: set of n target points, Y is an mxd matrix, d: geometric dimension
% t: number of subclusters to obtain from Y
% Returns:
% J: indices of target points closest to the gravity centers of the t subclusters

function [J] = GC_Sampling(Y,t)

% Sanity check
l = log(t)/log(2);
if(floor(l) ~= l)
   error('t must be a power of 2!');
end

if(t==1); l=1; end

% Get 1st generation of sons
[G{1},S{1},GGC] = Geo_Bal_Partition(Y);
% GGC: index of target point closest to the gravity center of Y

% Sons of further generations
for i=2:l
S{i} = {};
G{i} = {};
    for j = 1:size(S{i-1},2)   % number of clusters at previous generation
            [g,s] = Geo_Bal_Partition(S{i-1}{j});
            S{i} = cat(2,S{i},s);
            G{i} = cat(2,G{i},g);
    end
end

% Getting the indices of target points closest to gravity centers
for j=1:size(G{l},2)
    for i=1:size(Y,1)
    if(G{l}{j}'==Y(i,:))
    J(j)=i;
    end
    end
end

if(t==1)
J=GGC;
end

end


%% Function Geo_Bal_Partition
% Performs geometrically ballanced partition to divide a cluster into two clusters son
% Requires:
% S_y: cluster of points
% Returns:
% Son: list of two cluster sons
% G: contains the gravity centers of cluster sons
% GCC: index of target point closest to the gravity center of S_y

function [G,Sons,GGC] = Geo_Bal_Partition(S_y)

[n,~] = size(S_y);

g = S_y'*ones(n,1)/n;
Cov = S_y - g';
[~,~,v] = svd(Cov); v=v(:,1);

L = Cov*v;
b_1 = (L>0);
b_2 = (L<0);
Sons{1} = S_y(b_1,:);
Sons{2} = S_y(b_2,:);

% Getting the index of target point closest to the gravity center of S_y
v = (S_y - g.'); z = zeros(n,1);
    for i=1:n
        z(i)=norm(v(i,:));
    end
[~,GGC]=min(z);

% Getting indices of target points closest to the gravity centers of clusters son
for j=1:2
n_s = size(Sons{j},1);  G{j} = Sons{j}'*ones(n_s,1)/n_s;
v = (Sons{j} - G{j}.');
z = zeros(n_s,1);
    for i=1:n_s
        z(i)=norm(v(i,:));
    end
[~,f]=min(z);
G{j}=Sons{j}(f,:)';
end
end
