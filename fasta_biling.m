%   Solve the L1-penalized non-negative least-squares problem
%           minimize_{X,Y}  mu|X|+.5||S-XY'||^2
%           subject to  norm(Y)<=1, X>=0, Y>=0
%   using the solver FASTA.  Note that choosing mu=0 yields standard
%   non-negative least squares without an L1 penality.
%
%  Inputs:
%    S_EN, S_FR		: A matrix of data for English and French Words, size V_EN X N and V_FR X N, respectively
%    X0_EN, X0_FR  	: Initial guess of solutions (sparse codings), size V_EN X K and V_FR X K, respectively
%    Y0_EN, Y0_FR	: Initial guess of solutions (dictionary), size N X K
%    AL			: The alignment-score matrix, size V_EN X V_FR - has to be sparse!!
%    mu_en, mu_fr	: Scalar regularization parameters
%    opts		: Optional inputs to FASTA
%
%   For this code to run, the solver "fasta.m" must be in your path.
%
%   For more details, see the FASTA user guide, or the paper "A field guide
%   to forward-backward splitting with a FASTA implementation."
%
%   Copyright: Tom Goldstein, 2014.



function [ X_EN_sol, Y_EN_sol, X_FR_sol, Y_FR_sol, outs ] = ...
    fasta_biling(S_EN, X0_EN, Y0_EN, S_FR, X0_FR, Y0_FR, AL, mu_en, mu_fr, mu_xling, opts )

%%  Check that we have matrices
assert(isnumeric(S_EN) & isnumeric(X0_EN) & isnumeric(Y0_EN),'The inputs must be matrices (EN).')
assert(isnumeric(S_FR) & isnumeric(X0_FR) & isnumeric(Y0_FR),'The inputs must be matrices (FR).')

%  Check for 'opts'  struct
if ~exist('opts','var') % if user didn't pass this arg, then create it
    opts = [];
end

% Make unknowns into a single large matrix so that FASTA can handle them
Z0 = [X0_EN;Y0_EN;X0_FR;Y0_FR];
X_EN_rows = 1:size(X0_EN,1);
Y_EN_rows = size(X0_EN,1)+1:size(X0_EN,1)+size(Y0_EN,1);
X_FR_rows = size(X0_EN,1)+size(Y0_EN,1)+1:size(X0_EN,1)+size(Y0_EN,1)+size(X0_FR,1);
Y_FR_rows = size(X0_EN,1)+size(Y0_EN,1)+size(X0_FR,1)+1:...
    size(X0_EN,1)+size(Y0_EN,1)+size(X0_FR,1)+size(Y0_FR,1); 

%%  Define ingredients for FASTA
%  Note: fasta solves min f(Ax)+g(x).
%  f(Z) = .5 ||XY' - S||^2
f = @(Z)     .5*norm(Z(X_EN_rows,:)*Z(Y_EN_rows,:)' - S_EN,'fro')^2 + .5*norm(Z(X_FR_rows,:)*Z(Y_FR_rows,:)' - S_FR,'fro')^2 + .5*mu_xling*xling(Z, AL, X_EN_rows, X_FR_rows);
A = @(x) x;
At = @(x) x;
grad = @(Z) gradZ(Z,S_EN,S_FR,X_EN_rows,Y_EN_rows,X_FR_rows,Y_FR_rows, AL, mu_xling);


% g(z) = mu*|z|
g = @(Z) matrix1norm(Z(X_EN_rows,:))*mu_en + matrix1norm(Z(X_FR_rows,:))*mu_fr;
% proxg(Z,t):  Perform shrink on X, and ensure Z is positive.  Ensure Y
% lies in [0,1]

prox = @(Z,t) [max(Z(X_EN_rows,:)-t*mu_en,0) ; proxY(Z,Y_EN_rows); ...
    max(Z(X_FR_rows,:)-t*mu_fr,0); proxY(Z,Y_FR_rows)];

%% Call solver
[solution, outs] = fasta(A,At,f,grad,g,prox,Z0,opts);

X_EN_sol = solution(X_EN_rows,:);
Y_EN_sol = solution(Y_EN_rows,:);
X_FR_sol = solution(X_FR_rows,:);
Y_FR_sol = solution(Y_FR_rows,:);

end


%% Gradient function
function grad = gradZ(Z,S_EN,S_FR,X_EN_rows,Y_EN_rows,X_FR_rows,Y_FR_rows, AL, mu_xling)
X_EN = Z(X_EN_rows,:);
Y_EN = Z(Y_EN_rows,:);
X_FR = Z(X_FR_rows,:);
Y_FR = Z(Y_FR_rows,:);

[i,j,v] = find(AL);
[v_en,n_en] = size(X_EN);
[v_fr,n_fr] = size(X_FR);

% Calculate gradient from cross lingual term
dxling_en = zeros(v_en,n_en);
dxling_fr = zeros(v_fr,n_fr);

for counter = 1:length(v)

	curr_en = i(counter);
	curr_fr = j(counter);
	curr_al = v(counter);

	X_EN_i = X_EN(curr_en,:);
	X_FR_j = X_FR(curr_fr,:);
	dxling_en(curr_en,:) = dxling_en(curr_en,:) + curr_al*(X_EN_i - X_FR_j);
	dxling_fr(curr_fr,:) = dxling_fr(curr_fr,:) + curr_al*(X_FR_j - X_EN_i);

end

% Calculate gradient from monolingual term
diff_EN = X_EN*Y_EN'-S_EN;
dX_EN = diff_EN*Y_EN+dxling_en*mu_xling;
dY_EN = diff_EN'*X_EN;

diff_FR = X_FR*Y_FR'-S_FR;
dX_FR = diff_FR*Y_FR+dxling_fr*mu_xling;
dY_FR = diff_FR'*X_FR;

% Combine into the gradient for Z
grad = [dX_EN;dY_EN;dX_FR;dY_FR];
end 

%% 1-norm of a matrix
function norm1 = matrix1norm(Z)
norm1 = sum(abs(Z(:)));
end 

%% The cross lingual term in the function
function val = xling(Z, AL, X_EN_rows, X_FR_rows)
val = 0;
[i,j,v] = find(AL);

X_EN = Z(X_EN_rows,:);
X_FR = Z(X_FR_rows,:);

for counter = 1:length(v)
	curr_en = i(counter);
	curr_fr = j(counter);
	curr_al = v(counter);
	
	X_EN_i = X_EN(curr_en,:);
	X_FR_j = X_FR(curr_fr,:);

	val = val + curr_al*norm(X_EN_i-X_FR_j)^2;

end
end

%prox operator for characteristic function of the l2 ball
function proxYval = proxY(Z, Yrows)

Y = Z(Yrows,:);
YT = Y';
[r,c] = size(YT);
proxYvalT = zeros(r,c);
    

for i = 1:r
    proxYvalT(i,:) = YT(i,:)/max(norm(YT(i,:)), 1);
end

proxYval = proxYvalT';
end
