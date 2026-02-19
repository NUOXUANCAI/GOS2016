function plotsRiskPremia(T,K,lambda,CIu,CIb,lUncMod)
%This function produces four figures: 
%Figure 1: plot of estimated annualized MKT risk premia
%Figure 2: plot of estimated annualized SMB risk premia
%Figure 3: plot of estimated annualized HML risk premia
%Figure 4: plot of estimated annualized MOM risk premia
%These figures correspond to the Figures 1,2 and 3 in the paper.

Dates=(1:1:T);

ly=lambda;      %Re-definition
lCondMod=mean(ly,2);    %average over time of time-varying risk premia

%Economic crises
x2=[Dates(1,66),Dates(1,77)];   %Dic 1969 - Nov 1970
x3=[Dates(1,113),Dates(1,129)]; %Nov 1973 - Mar 1975
x4=[Dates(1,187),Dates(1,193)]; %Jan 1980 - Jul 1980
x5=[Dates(1,205),Dates(1,221)]; %Jul 1981 - Nov 1982
x6=[Dates(1,313),Dates(1,321)]; %Jul 1990 - Mar 1991
x7=[Dates(1,441),Dates(1,448)]; %Mar 2001 - Ott 2001
x8=[Dates(1,522),Dates(1,540)]; %Dic 2007 - Jun 2009

MINy=-20; 
MAXy=40;
y=MAXy.*ones(2,1);
ylim([MINy MAXy]);
xlim([Dates(1,1) Dates(1,T)]);
XT=[(Dates(1,7):60:Dates(1,T))';Dates(1,T)];
XTL=['65';'70';'75';'80';'85';'90';'95';'00';'05';'10'];

for k=1:K;
    figure(k);
    hold on;
    h2=area(x2,y,MINy); h3=area(x3,y,MINy); h4=area(x4,y,MINy);
    h5=area(x5,y,MINy); h6=area(x6,y,MINy); h7=area(x7,y,MINy); h8=area(x8,y,MINy);
    h=[h2;h3;h4;h5;h6;h7;h8];
    set(h,'FaceColor',[.9,0.9,0.9]);
    set(h,'LineStyle','none')
    plot(Dates(1,2:T),ly(k,:), 'LineWidth',2.2);
    plot(Dates(1,2:T),CIb(k,:), '.k','MarkerSize',3.7);
    plot(Dates(1,2:T),CIu(k,:),'.k','MarkerSize',3.7);
    plot(Dates(1,2:T),lUncMod(k,1).*ones(T-1,1),'--r','LineWidth',1.5);
    plot(Dates(1,2:T),lCondMod(k,1).*ones(T-1,1),'b','LineWidth',1.5);
    set(gca,'XTick',XT,'XTickLabel',XTL);
    if k==1;
        title('$\hat{\lambda}_{m,t}$','FontSize',18,'Interpreter','latex');
        hold off;
    elseif k==2;
        title('$\hat{\lambda}_{smb,t}$','FontSize',18,'Interpreter','latex');
        hold off;
    elseif k==3;
        title('$\hat{\lambda}_{hml,t}$','FontSize',18,'Interpreter','latex');
        hold off;
    elseif k==4;
        title('$\hat{\lambda}_{mom,t}$','FontSize',18,'Interpreter','latex');
        hold off;
    end;
end;

    

