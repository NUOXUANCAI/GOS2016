function [C,Ja]=Cmatrix(nu,p,q,K,d1,d2)
%This function computes matrix C_{\nu} defined in section 3.2

D=duplication(p);       %(p*p)x(p*(p+1)/2)
Dplus=(inv(D'*D))*D';   %Moore-Penrose inverse of D, (p*(p+1)/2)x(p*p)
W1=commutation(p*(p+1)/2,p*K);  %commutation matrices
W2=commutation(p,p);
W3=commutation(p*q,p*K);
W4=commutation(p,q);
E1=[eye(d1);zeros(d2,d1)]; %selection matrix, dxd1
E2=[zeros(d1,d2); eye(d2)];%selection matrix, dxd2

J11=W1*(kron(eye(K),(kron(eye(p),Dplus)*kron(W2,eye(p))*kron(eye(p),vec(eye(p))))));
%J11 is a (p^2*K*(p+1)/2)xp*K matrix
J22=W3*(kron(eye(K),(kron(eye(p),W4)*kron(W4,eye(p))*kron(eye(q),vec(eye(p))))));
%J22 is a (p^2*q*K)xq*K matrix
Ja=[J11, zeros(p^2*K*(p+1)/2,q*K); zeros(p^2*q*K,p*K), J22];

C=(E1'-kron(eye(d1),nu')*Ja*E2')';  %C is a dxd1 matrix


