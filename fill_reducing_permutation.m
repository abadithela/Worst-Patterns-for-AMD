%% Finding the best fill reducing permutation by brute-force for small N
% Apurva Badithela
% July 17th, 2017

% We desire to find the best fill-in reducing permutation for a symmetric
% positive definite sparse square matrix A of size n. For a given n, the
% number of possible sparsity patterns are 2^{n(n+1)/2}.

% Given an A, finding a P for which fill-in is minimum:
% A permutation matrix of size n can be represented as a permutation of the
% sequence of first n natural numbers. This permuted sequence represents
% the order in which the rows of an identity matrix are swapped, givign us
% the permutation P. Even for small n, the total number of permutations
% will grow quickly, therefore we desire to find a method in which all
% permutations will processed in a lexicographic order. For this purpose,
% we are implementing the Heap's algorithm to go through all possible P's:

% Initialization:
function [g_amd, g_colamd, g_symrcm, MIN_FILL, fill_reduce] = fill_reducing_permutation(A);
n = size(A,1); % Size of symmetric matrix A
% TODO: See if there is a way to keep A sparse and perform LDL on A.
% A = full(generatesparseSPDmatrix(n, 0.5)); 

nonzeros = nnz(tril(A,-1));
A = A + n*speye(n); % Make it positive definite
MAX_FILL = n*(n-1)/2; % Maximum possible fill-in
MIN_FILL = MAX_FILL;
p = 1:n; % Initial permutation
c = zeros(1, n);
index = 0;
I = speye(n);
fill_reduce = []; % List of matrices with different permutations but still minimum fill

% Check fill-in for initial permutation -
P0 = I;
p_opt = p;
L = chol(P0*A*P0', 'lower');
fill = nnz(L) - nonzeros;
if(fill < MIN_FILL)
    MIN_FILL = fill;
    p_opt = p;
    fill_reduce = [p_opt]; % Optimal Pattern
end

% Calling Heap function to generate next lexicographic function
while(index<n)
    [p, c, index] = heap(n, p, c, index);
    P = I(p,:);
    L = chol(P*A*P', 'lower'); % [L, D] = ldl(X) returns L such that X = L*D*L'
    fill = nnz(tril(L,-1))-nonzeros;
    if(fill <= MIN_FILL)
        MIN_FILL = fill;
        p_opt = p;
        fill_reduce = [p_opt];
%     elseif(fill == MIN_FILL)
%         fill_reduce = [fill_reduce; p];
    end
end

% Comparison with AMD:
p_amd = amd(A);
L_amd = chol(A(p_amd, p_amd));
fill_amd = nnz(L_amd) - nonzeros;
g_amd = fill_amd - MIN_FILL;

% Comparison with Reverse Cuthill McKee:
p_symrcm = symrcm(A);
L_symrcm = chol(A(p_symrcm, p_symrcm));
% [L_symrcm,~, ~] = ldl(A(p_symrcm, p_symrcm));
fill_symrcm = nnz(L_symrcm) - nonzeros;
g_symrcm = fill_symrcm - MIN_FILL;

% Comparison with COLAMD:
p_colamd = colamd(sparse(A));
L_colamd = chol(sparse(A(p_colamd, p_colamd)));
% [L_colamd, ~] = lu(sparse(A(p_colamd, p_colamd)));
fill_colamd = nnz(L_colamd) - nonzeros;
g_colamd = fill_amd - MIN_FILL;

end

%% Some details - brute forcing is feasible only for size(P)<=10. For n = 11 onwards, the average time required to compute P increases to 50 minutes. For n = 20, we will require about a day to brute-force.
% Finding the most evil matrix A will be even harder.
% Average number of matrices processed per second ~ 6600.

% size = 10:17;
% projected_sol_time = factorial(size)./86400./6600e8; % In minutes. Solving a matrix of size 20 would take 11 hours on regular Matlab on PC