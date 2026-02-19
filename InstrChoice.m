function [Z,Zi]=InstrChoice(q,T,n,Ti,BM,I)
% This function loads and standardizes common instruments and
% asset-specific instruments.

% Common Instruments
load Wspace_Instrum 
ds=(DefSpread-mean(DefSpread))./std(DefSpread);       
ts=(TermSpread-mean(TermSpread))./std(TermSpread);    
Z=[ones(T,1),ds,ts];  	

% Asset-specific instruments
Zi=zeros(T,q,n);   
for i=1:n;
    %Standardiz. of BM (book-to-mkt)
    m=sum(BM(:,i))/Ti(i);
    s=sqrt(sum((BM(:,i)-I(:,i).*repmat(m,T,1)).^2)/Ti(i));
    z1=I(:,i).*((BM(:,i)-m)/s);
    Zi(:,:,i)=z1;
end;

