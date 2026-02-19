
%% ========================================================================
%%  TIME-VARYING RISK PREMIUM IN LARGE CROSS-SECTIONAL EQUITY DATASETS 
%% P.Gagliardini, E.Ossola, O.Scaillet
%% June 2014
%% ========================================================================
clc;
clear;

%% INPUT
% Choice of dataset
dataset=1;          % 1=25FF portfolios, 2=44 Indu. portfolios, 3=Individual stocks
% Number of factors (exclude the constant)
K=4;
% Model 
Mod=1;              % Mod=1 --> the code estimates time-invariant and time-varying models 
                    % Mod=0 --> the code estimates only time-invariant model

%% Load Datasets
% Each dataset contains:
% R=returns, I=indicator for unbalanced characteristic, BM=book-to-mkt
% Ti= time-series observations for each asset, T=TS size, n=CS size
if dataset==1;     load Wspace_25FF;            %25 Fama-French portfolio Jul64-Dec09
elseif dataset==2; load Wspace_44Indu;          %44 Indu. portfolio Jul64-Dec09
elseif dataset==3; load Wspace_CRSPCMST_ret;    %Dataset after matching CRSP and Compustat Jul64-Dec09
end;
% Factors matrix
load Wspace_Fact F
F=F(:,1:K).*10;                                 %Rescaled factors    
% Excess Returns
load RiskFree;                                  %Free-risk rate
ExcR=(I==1).*(R-repmat(Rf,1,n))+0;              %Txn matrix of excess returns

%%  OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		TIME-INVARIANT MODEL		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Estimation of no-time varying risk premia, nu, and statistic (see Section 3) 
    [n,nchi,nchiSTAT,Table1,Table2,Table4]=EstTimeInvariantMod(n,T,K,ExcR,[],F,I,Ti,dataset);
    %Display annualized results
    display(n);      %cross-sectional dimension
    display(nchi);   %trimmed cross-sectional dimension (Estimation)
    display(Table1); %no-time varying risk premia and its confidence intervals (CI): [riskpremia,CIbelow,CIup]
    display(Table2); %no-time varying nu and its confidence intervals: [nu,CIbelow,CIup]
    if dataset==3; display(nchiSTAT); end;      %trimmed cross-sectional dimension (Statistics - Table 4)
    display(Table4); %statistics and p-values for invariant case: [stat(Qe),pv;stat(Qa),pv]
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		TIME-VARYING MODEL		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Mod==1;
    %Number of instruments
    p=3;                %number of common instruments Z
    q=1;                %number of asset-specific instruments Zi
    %load Instruments
    [Z,Zi]=InstrChoice(q,T,n,Ti,BM,I);            
    %Estimation of time varying risk premia, nu, and statistic (see Section 3) 
    [n,nchi,nchiSTAT,Table3,Table5,RiskPremia,CIbRiskPremia,CIuRiskPremia,nuT,CIbnuT,CIunuT,TimeVarTest]=EstTimeVariantMod(K,p,q,n,T,ExcR,[],I,F,Z,Zi,dataset);
    %Display annualized results
    display(n);         %cross-sectional dimension
    display(nchi);      %trimmed cross-sectional dimension (Estimation)
    display(Table3);    %estimated components and CI: [vec(Fhat'),CIvecFbelow,CIvecFup,nu,CInubelow,CInuup]
    if dataset==3; display(nchiSTAT); end;         %trimmed cross-sectional dimension (Statistics - Table 5)
    display(Table5);    %statistics and p-values for variant case: [stat(Qe),pv;stat(Qa),pv]
    %Plots of time-varying risk premia and the paths of nu_{t}
    plotsRiskPremia(T,K,RiskPremia,CIbRiskPremia,CIuRiskPremia,Table1(:,1));
    plotsNu(T,K,nuT,CIunuT,CIbnuT,Table2(:,1));
    display(TimeVarTest);%Time variation test statistics and p-value, [XiF,pvF;XiNU,pvNU]
end;
