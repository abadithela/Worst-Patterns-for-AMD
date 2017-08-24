%% Finding most evil sparse matrix A for which AMD does poorly:
% Apurva Badithela
% July 18th, 2017

close all
clear all

% We go through all possible 2^n(n-1)/2 sparsity patterns of A:
tic
sz = 5;
nCk = 1:sz*(sz-1)/2;
N = length(nCk);
Z = zeros(1, N);
GMAX = 0; % Min. possible fill-in. AMD is doing as well as minimum fill-in
Evil_Patterns = [];
Fill_Reduce = [];


for K = 1:N
    m = 0;
    h = K;
    iteration = 1;
    nonzeros = K;
    pattern = 1:K;
    vec = Z;
    vec(pattern) = 2;
    A = sparse(fillMatrix(sz, vec)); % First Combination
    [GMAX, min_fill, Evil_Patterns, Fill_Reduce] = AnalyzeEvilPattern(A, GMAX, Evil_Patterns, Fill_Reduce);
    
    while(iteration > 0)
        [pattern, m, h, iteration] = GetNextCombination(N, K, pattern, m, h, iteration);
        vec = Z;
        vec(pattern) = 2;
        A = sparse(fillMatrix(sz, vec)); % First Combination
        [GMAX, min_fill, Evil_Patterns, Fill_Reduce] = AnalyzeEvilPattern(A, GMAX, Evil_Patterns, Fill_Reduce);
    end
end
toc

for ii = 1:1
    figure(ii)
    %caption = sprintf('Worst matrices for AMD (n = 5, g_{amd} = %d)', GMAX);    
    for jj = 1:9
        subplot(3,3,jj);
        spy(Evil_Patterns(:,:,(ii-1)*9 + jj));
    end
    %suptitle(caption);
end
