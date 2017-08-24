function A  = generatesparseSPDmatrix(n, density)
% Generate sparse random symmetric matrix:
A = sprandsym(n, density);
% Adding identity to make A positive definite
A = A + n*speye(n); 
end