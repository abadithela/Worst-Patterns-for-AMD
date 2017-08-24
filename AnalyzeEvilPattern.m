function [GMAX, min_fill, Evil_Patterns, Fill_Reduce] = AnalyzeEvilPattern(A, GMAX, Evil_Patterns, Fill_Reduce)
%% Function to check if input A is the most evil pattern:

[g_amd, g_colamd, g_symrcm, MIN_FILL, fill_pattern] = fill_reducing_permutation(A);
min_fill = [];
if(g_amd > GMAX)
    GMAX = g_amd;
    Evil_Patterns = [full(A)];
    Fill_Reduce = [fill_pattern];
    min_fill = MIN_FILL;
elseif(g_amd == GMAX)
    Evil_Patterns = cat(3,Evil_Patterns, full(A));
    % Fill_Reduce = cat(3, Fill_Reduce, fill_pattern);
end

end