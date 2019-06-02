%% CUR approximation with Gravity points criterion
% Requires:
% X,Y: Source and target points, given as matrices of size (mxd) and (nxd)
% respectively, where d is the geometric dimension
% k: fixed approximation rank.
% fun: kernel function, e.g. Laplacian kernel: fun = @(x,y) -1/(2*pi)*log(norm(x-y);
% Exponential kernel: fun = @(x,y) exp(1i*norm(x-y))/norm(x-y);
% Gravitation kernel:  fun = @(x,y) 1/(4*pi*norm(x-y));
% Returns:
% CUR: a rank-k approximation of matrix A(i,j)=fun(X(i,:),Y(j,:)).
% A \approx CxUxR, where C, R, U are complex matrices of size (mxk), (kxn), (kxk)

function [CUR] = CUR_GCS(fun,X,Y,k)

m = size(X,1);  n = size(Y,1);

% Finding t: number of sampling columns
l = nextpow2(k);
% t = pow2(l);
if(k > pow2(l-1) && k>2 )
    t = pow2(l+1);
else
    t = pow2(l);
end
if(k==1); t = 1; end

C=zeros(m,t);
R=zeros(k,n);

% Decompose target domain into t subdomains
[J] = GC_Sampling(Y,t);
% [J] = NN_Sampling(Y,X,t); % Alternatively use Nearest-Neighbors sampling

% Form matrix C of sampling columns, C is of size mxt
for i=1:size(X,1)
    for j=1:t
        C(i,j) =  fun(X(i,:),Y(J(j),:));
    end
end

[Q,~,p_c]=qr(C,'vector');
Q=Q(:,1:k);
C=C(:,p_c(1:k));

% Get column indices
[~,~,p_r]=qr(Q.','vector');
I=p_r(1:k);

% Form Matrix R
for i=1:k
    for j=1:size(Y,1)
        R(i,j) =  fun(X(I(i),:),Y(j,:));
    end
end

% Construct the CUR rank-k approximation
G=C(I,:);
CUR=C*(G\R); % Use Algorithm 1 for computing this skeleton approximation in order to better handle and control the selected indices

return
