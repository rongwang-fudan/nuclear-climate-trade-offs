% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
tic
clear;

EndSav=1;

load('files\codelist.dat','-mat');
load('..\nuclear\NEI_sentest.dat','-mat');
anslist_mean=zeros(5,4);
anslist_median=zeros(5,4);
anslist_SE=zeros(5,4);
anslist_Skew=zeros(5,4);

pollist=[1 2 3 4];% 1.5C 2C NDC CUR
testlist=[1 4];%

for t=1
    for p=1:4
        if t==2
            pollist=[1 2 3 4];
        end
        for k=1:2
            if t>1
                NEI_1_0=nvmc(4*(testlist(t)-2)+4+pollist(p),1:39);
                NEI_2_0=nvmc(4*(testlist(t)-2)+4+pollist(p),39+1:39+39*5);
                NEI_3_0=nvmc(4*(testlist(t)-2)+4+pollist(p),39+39*5+1:39+39*5+39*5*39);
            else
                NEI_1_0=nvmc(4*(testlist(t)-1)+pollist(p),1:39);
                NEI_2_0=nvmc(4*(testlist(t)-1)+pollist(p),39+1:39+39*5);
                NEI_3_0=nvmc(4*(testlist(t)-1)+pollist(p),39+39*5+1:39+39*5+39*5*39);
            end
            NEI_1_1=[];
            NEI_2_1=[];
            NEI_3_1=[];
            NEI_1_2=[];
            NEI_2_2=[];
            NEI_3_2=[];
            n11=1;n21=1;
            for n=1:39
                if NEI_1_0(n)<-100
                    NEI_1_1(n11)=NEI_1_0(n);
                    n11=n11+1;
                else
                    NEI_1_2(n21)=NEI_1_0(n);
                    n21=n21+1;
                end
            end
            n12=1;n22=1;
            for n=1:39*5
                if NEI_2_0(n)<-100
                    NEI_2_1(n12)=NEI_2_0(n);
                    n12=n12+1;
                else
                    NEI_2_2(n22)=NEI_2_0(n);
                    n22=n22+1;
                end
            end
            n1=1;n2=1;
            for n=1:39*5*39
                if NEI_3_0(n)<-100
                    NEI_3_1(n1)=NEI_3_0(n);
                    n1=n1+1;
                else
                    NEI_3_2(n2)=NEI_3_0(n);
                    n2=n2+1;
                end
            end

            if k==1
                NEI_1=NEI_1_1;
                NEI_2=NEI_2_1;
                NEI_3=NEI_3_1;
                edges = [-2100:100:-100];
                edgesplot0=[-290:10:-100]-5;
                kn=10;
            else
                NEI_1=NEI_1_2;
                NEI_2=NEI_2_2;
                NEI_3=NEI_3_2;
                edges = [-100:10:100];
                edgesplot0=edges(2:end)-5;
                kn=1;
            end



            [N1] = histcounts(NEI_1,edges);
            N21=N1./(39*kn);


            [N2] = histcounts(NEI_2,edges);
            N22=N2./(39*5*kn);

            [N3] = histcounts(NEI_3,edges);
            N23=N3./(39*5*39*kn);

            N31(1+20*(k-1):20+20*(k-1))=N21;
            N32(1+20*(k-1):20+20*(k-1))=N22;
            N33(1+20*(k-1):20+20*(k-1))=N23;
            edgesplot(1+20*(k-1):20+20*(k-1))=edgesplot0;
            if t==1
                anslist_mean(1,p)=mean(NEI_1_0);
                anslist_mean(2,p)=mean(NEI_2_0);
                anslist_mean(3,p)=mean(NEI_3_0);
                anslist_median(1,p)=median(NEI_1_0);
                anslist_median(2,p)=median(NEI_2_0);
                anslist_median(3,p)=median(NEI_3_0);
                anslist_SE(1,p)=std(NEI_1_0);
                anslist_SE(2,p)=std(NEI_2_0);
                anslist_SE(3,p)=std(NEI_3_0);
                anslist_Skew(1,p)=sum((NEI_1_0-mean(NEI_1_0)).^3)/(std(NEI_1_0)^3)/size(NEI_1_0,2);
                anslist_Skew(2,p)=sum((NEI_2_0-mean(NEI_2_0)).^3)/(std(NEI_2_0)^3)/size(NEI_2_0,2);
                anslist_Skew(3,p)=sum((NEI_3_0-mean(NEI_3_0)).^3)/(std(NEI_3_0)^3)/size(NEI_3_0,2);
                [a2,b2]=find(NEI_3_0<-10);
                ans0(p)=size(a2,2)/7605*100;
            else
                anslist_mean(4,p)=mean(NEI_3_0);
                anslist_median(4,p)=median(NEI_3_0);
                anslist_SE(4,p)=std(NEI_3_0);
                anslist_Skew(4,p)=sum((NEI_3_0-mean(NEI_3_0)).^3)/(std(NEI_3_0)^3)/size(NEI_3_0,2);
            end
        end
        if t==1
            subplot(4,2,p)
            save_bar(p,:)=log10(N33)+5;
            b=bar(edgesplot,log10(N33)+5);hold on;
            b.FaceColor=[191 191 191]./255;
            ylim([0 5]);
            xlim([-300 100]);
        end
    end

end

load('files\codelist2.dat','-mat');
nucroulist=zeros(5,39*39);
mlist=ones(5);
for i=39+39*5+1:39*39*5+39+39*5
    if codelist2(i,7)==0
        nucroulist(1,mlist(1))=i;
        mlist(1)=mlist(1)+1;
    elseif codelist2(i,7)==0.25
        nucroulist(2,mlist(2))=i;
        mlist(2)=mlist(2)+1;
    elseif codelist2(i,7)==0.5
        nucroulist(3,mlist(3))=i;
        mlist(3)=mlist(3)+1;
    elseif codelist2(i,7)==0.75
        nucroulist(4,mlist(4))=i;
        mlist(4)=mlist(4)+1;
    elseif codelist2(i,7)==1
        nucroulist(5,mlist(5))=i;
        mlist(5)=mlist(5)+1;
    end
end

Pines7=0.00159;
pollist=[1 2 3 4];% 1.5C 2C NDC CUR6t
testlist=[2 3 5];%1-2for fig3a  305 for fig3b
for t=1:3
    for p=1:4
        for m=1:5
            nvnei=nvmc(4*(testlist(t)-2)+4+pollist(p),nucroulist(m,:));
        if t==1
        NEI(1,m,p)=mean(nvnei);
        NEI(2,m,p)=median(nvnei);
        NEI(3:13,m,p)=prctile(nvnei,[5 10 20 30 40 50 60 70 80 90 95]);
        else
        NEI(1,m,p)=mean(nvnei);
        NEI(2,m,p)=median(nvnei);
        end
        end


subplot(4,2,p+4)
colorrb2=[255 88 137
   255 140 195
   255 188 225
   255 223 248
   234 215 255
   238 243 248
   212 212 255
   186 186 255
   156 156 255
   138 138 255
   83 83 255]./255;
if t==1
    colors=[17 88 212]./255;
elseif t==2
    colors=[238 88 96]./255;
else
    colors=[155 187 89]./255;
end
if t==1
for i=11:-1:1
    plot([1:5],NEI(i+2,:,p),'-','LineWidth',1,'Color',colorrb2(12-i,1:3)); hold on;
    ylim([-150 75]);

end
    plot([1:5],NEI(1,:,p),'LineStyle','--','LineWidth',1,'Color',colors); hold on;
    plot([1:5],NEI(2,:,p),'LineStyle','-','LineWidth',1,'Color',colors); hold on;
    ylim([-150 75]);
else
    plot([1:5],NEI(1,:,p),'LineStyle','--','LineWidth',1,'Color',colors); hold on;
plot([1:5],NEI(2,:,p),'LineStyle','-','LineWidth',1,'Color',colors); hold on;
ylim([-150 75]);
end

    end

end

