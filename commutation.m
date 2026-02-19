function k=commutation(n,m)
%The function computes the commutation matrix
%A is a matrix nxm
%k(n,m)*vec(A)=vec(A')

k = reshape(kron(vec(eye(n)), eye(m)), n*m, n*m);

