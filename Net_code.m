clear all; clc;
load('sptimes.mat')
dt=0.005;
mint=-0.1900;
maxt=0.3;
N=45;
STTC=zeros(N,N);

%calculating correlation matrix using STTC------------
for Node_i= 1:N
    for Node_j=1:N
        TA=0;TB=0;
        PA=0;PB=0;
        %calculating T_A------------------------------
        HA=zeros(2*length(sptimes{Node_j}),1);
        for k=1:length(sptimes{Node_j})
            HA(2*k-1)=sptimes{Node_j}(k)-dt;
            HA(2*k)=sptimes{Node_j}(k)+dt;
            if k>1
                if HA(2*k-1)<HA(2*(k-1))
                TA=TA+2*dt-(HA(2*(k-1))-HA(2*k-1));
                else
                 TA=TA+2*dt;
                end
            else
                TA=TA+2*dt;
            end
        end
        TA=TA/(maxt-mint);
        %calculating T_B-------------------------------
        HB=zeros(2*length(sptimes{Node_i}),1);
        for k=1:length(sptimes{Node_i})
            HB(2*k-1)=sptimes{Node_i}(k)-dt;
            HB(2*k)=sptimes{Node_i}(k)+dt;
            if k>1
                if HB(2*k-1)<HB(2*(k-1)) 
                    TB=TB+2*dt-(HB(2*(k-1))-HB(2*k-1));
                else
                    TB=TB+2*dt;
                end
            else
                TB=TB+2*dt;
            end
        end
        TB=TB/(maxt-mint);
        %calculating P_A-------------------------------     
        for k=1:length(sptimes{Node_j})
            for l=1:length(sptimes{Node_i})
                if sptimes{Node_j}(k)>HB(2*l-1) && sptimes{Node_j}(k)<=HB(2*l)
                    PA=PA+1;
                    break
                end
            end
        end
        PA=PA/length(sptimes{Node_j});    
        %calculating P_B-------------------------------
        for k=1:length(sptimes{Node_i})
            for l=1:length(sptimes{Node_j})
                if sptimes{Node_i}(k)>HA(2*l-1) && sptimes{Node_i}(k)<=HA(2*l)
                    PB=PB+1;
                    break
                end
            end
        end
        PB=PB/length(sptimes{Node_i});
        STTC(Node_i,Node_j)=0.5*((PA-TB)/(1-PA*TB)+(PB-TA)/(1-PB*TA));
    end
end
STTC=STTC-eye(N,N);
STTC_M=STTC;

%degress distribution in weighted graph---------------
% degree_k_wg=sum(STTC,2);
% meanV_degree_k_wg=mean(degree_k_wg);
% histogram(degree_k_wg,'BinWidth',0.5)
% annotation('textbox',[0.7,0.7,0.1,0.1],'string',sprintf('k_{mean} = %0.2f',meanV_degree_k_wg))
%making unweighted network----------------------------
epsilon=0.65;
for i=1:N
    for j=1:N
        if STTC_M(i,j)>=epsilon
            STTC(i,j)=1;
        else
            STTC(i,j)=0;
        end
    end
end
%calculating degree distribution/mean valu------------ 
degree_k=sum(STTC,2);
meanV_degree_k=mean(degree_k);
edge=unique(degree_k);
for i=1:length(edge)
    edge(i)=sum(edge(i)==degree_k);
end
[counts,binedge]=histcounts(degree_k,unique(degree_k));

%graph/plot/histogram----------------------------------
%plot(unique(degree_k),edge)
figure;
histogram(degree_k,'BinWidth',0.5)
annotation('textbox',[0.7,0.7,0.1,0.1],'string',sprintf('<k> = %0.2f',meanV_degree_k))
xlabel('k')
ylabel('Number of Nodes')
title('Degree Distribution (\epsilon=0.65)');
figure;
plot(graph(STTC))
title('Network Visualization (\epsilon=0.65)');

%calculating giant component---------------------------
Epsilon=(0:0.05:1);
giant_c=zeros(length(Epsilon),1);
for e=1:length(Epsilon)
for i=1:N
    for j=1:N
        if STTC_M(i,j)>=Epsilon(e)
            STTC(i,j)=1;
        else
            STTC(i,j)=0;
        end
    end
end
bins=conncomp(graph(STTC));
giant_c(e)=max(histcounts(bins,'BinMethod','integer'));
end
figure;
plot((1:-0.05:0),giant_c/N,'O')
xlabel('1-\epsilon,(connectivity)')
ylabel('g,(giant component)')
title('');
%calculcating critical exponent-------------------------
c_critical=0.25;
figure;
% to determine the critical exponent,should only use the data points near critical
% C0 which is 0.25. for this special case we consider n=7 to 13.
C=Epsilon(7:length(Epsilon))./(c_critical*ones(1,length(Epsilon)-6))-ones(1,length(Epsilon)-6);
G=flipud(giant_c(1:length(Epsilon)-6))/N;
plot(log10(C(1:6)),log10(G(1:6)),'o','DisplayName','data')
hold on
exponent=polyfit(log10(C(1:6)),log10(G(1:6)),1);
fittedcurve=polyval(exponent,log10(C(1:6)));
plot(log10(C(1:6)),fittedcurve,'DisplayName','Fitted Curve')
annotation('textbox',[0.2,0.7,0.1,0.1],'string',sprintf('\\beta = %0.2f',exponent(1)))
legend('show','Location','NorthWest')
xlabel('log(|1-c/c_{0}|)')
ylabel('log(g)')