 function V=varcovNW(T,L,y)
% This function computes the  Newey-West estimator of the var-cov matrix
% T=number of time serie observations
% L=number of lags
% y=vector of observations

m=(T-L)^(1/4);    % finite "bandwidth" parameter (ref. Hayashi (2000))
m=m-rem(m,1);     % integer part of m
y=y(1:T-L,:);     % (T-L)xK (lagged)
ybar=y-repmat(mean(y),T-L,1); % vector y-E(y)

t=1;              
gamma0=ybar(t:T-L,:)'*ybar(t:T-L,:)/(T-L);  % Cov(yt,yt), T-LxT-L matrix

sum1=0; 
for j=1:m;    
    gamma=(ybar(j+1:T-L,:)'*ybar(1:T-L-j,:))/(T-L); %Cov(st,s(t-j))
    sum1= sum1 + (1-j/(m+1))*(gamma+gamma');       
end;

% The Newey-West estimator of the var-cov matrix
V=gamma0+sum1;
