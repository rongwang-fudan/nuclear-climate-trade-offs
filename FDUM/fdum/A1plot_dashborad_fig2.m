% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2023.7.5
tic
clear;

load('files\cndata134.dat','-mat'); cndata=cndata2; clear cndata2;
load('files\cou_iform.dat','-mat'); % 1: id for 222 countries; 2: 2 developing/ 1 developed; 3: 12 region id; 4 OECD; 5 id for 112 countries; 6 pi temperature
cou_iform=cou_iform2; clear cou_iform2;
WGR=load('files\WGR_GHG.txt');
load('..\nuclear\output_pop.dat','-mat'); % output_pop=zeros(400,cn_num-1);
cn_num=size(cndata,1);
EndSav=1;
cmap2=jet(21);
realtime=[45:48 49:12:121 122:196];% 2015-2100
realtime2=[45:48 49:12:121 122:296];% 2015-2200


ans0=zeros(5,4);
for km=1:2
    if km==1
        TW=1;%1for 0TW and 11 for 1TW
        color=jet(6);
    else
        TW=11;
        color=jet(6);
    end
    explist=[1:4];
    for ex=1:4
        T=explist(ex);
        load(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_var_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_emi_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');

        dW=zeros(16,2);
        for s=1:16
            dW(s,1)=dW(s,1)+sum(output_nuclear(1:(cn_num-1),s,2),1); % global Gt CO2 abated by nuclear power
            for cn=1:(cn_num-1)
                dW(s,2)=dW(s,2)+(output_nuclear(cn,s,3)-output_nuclear(cn,1,3))*output_pop(400,cn); % consumption trillion $
            end
        end
        nv(:,ex)=dW(:,2);


        subplot(4,4,1) %co2 emission
        if ex==1
            emi_neg=15;
        elseif ex==2
            emi_neg=5;
        else
            emi_neg=0;
        end
        for t=1:86
            if realtime(t)<121
                emi_neg0(t) = 0; % deploying negative emissions since 2025
            else
                emi_neg0(t) = emi_neg-emi_neg*exp(-(realtime(t)-121)^2/(-(2050-2025)^2/log(0.1))); % deploying emi_neg after 2025
            end
        end
        for k=2:6
            m=1;
            clear ghgplot
            for q=1:36
                if WGR(k,q)==-999
                else
                    if k==6
                        ghgplot2(q,1)=2014+q;
                        ghgplot2(q,2)=WGR(k,q)-(WGR(k,q)*(1-0.685)*(1-0.5/36*q));
                        m=m+1;
                    elseif k==5
                        ghgplot2(q,1)=2014+q;
                        ghgplot2(q,2)=WGR(k,q)-(WGR(k,q)*(1-0.685)*(1-0.5/36*q));
                        m=m+1;
                    elseif k==4
                        ghgplot2(q,1)=2014+q;
                        ghgplot2(q,2)=WGR(k,q)-(WGR(k,q)*(1-0.685)*(1-0.6/36*q));
                        m=m+1;
                    elseif k==3
                        ghgplot2(q,1)=2014+q;
                        ghgplot2(q,2)=WGR(k,q)-(WGR(k,q)*(1-0.685)*(1-0.6/36*q));
                        m=m+1;
                    else
                        ghgplot2(q,1)=2014+q;
                        ghgplot2(q,2)=WGR(k,q)-(WGR(k,q)*(1-0.685)*(1-0.7/36*q));
                        m=m+1;
                    end
                end
            end
            plot(ghgplot2(36,1),ghgplot2(36,2),'o','MarkerSize',5,'MarkerEdgeColor','none','MarkerFaceColor',cmap2(k,:)); hold on;
            policynum(k-1)=ghgplot2(36,2);
        end
        for i=1:max(size(realtime))
            co2plot(i)=output_var(realtime(i),TW,5)-output_var(realtime(i),TW,27)+emi_neg0(i);
        end
        if km==1
            plot([2015:2100],co2plot,'LineStyle','-','LineWidth',1,'Color',color(ex,:)); hold on; % 2025
        elseif km==2
            plot([2015:2100],co2plot,'LineStyle','--','LineWidth',1,'Color',color(ex,:)); hold on; % 2025
        end
        axis([2015 2100 0 50]);
        set(gca,'xtick',2015:10:2100);


        tt=[2015:2200];
        subplot(4,4,2)% warming
        if km==1
            plot(tt,output_var(realtime2,TW,1),'LineStyle','-','LineWidth',1,'Color',color(ex,1:3)); hold on; % warming
        elseif km==2
            plot(tt,output_var(realtime2,TW,1),'LineStyle','--','LineWidth',1,'Color',color(ex,1:3)); hold on; % warming
        end
        axis([2015 2200 0 3]);

        subplot(4,4,3)% nuclear capacity factor
        if km==1

        elseif km==2
            plot([2025:2100],output_var(122:197,TW,25),'LineStyle','-','LineWidth',1,'Color',color(ex,1:3)); hold on;
        end
        axis([2025 2100 0 1]);

        for t=1:400
            percapcons(t,TW)=output_var(t,TW,4)*(1-output_var(t,TW,6))/sum(output_pop(t,1:(cn_num-1)),2)*1000;
            prob=1;
            for tip_i=1:5
                prob=prob*output_var(t,TW,tip_i+18);
            end
            tip(t,TW)=1-prob;
        end
        subplot(4,4,4)% probability of tipping
        if km==1
            plot(tt,tip(realtime2,TW),'LineStyle','-','LineWidth',1,'Color',color(ex,1:3)); hold on;
        elseif km==2
            plot(tt,tip(realtime2,TW),'LineStyle','--','LineWidth',1,'Color',color(ex,1:3)); hold on;
        end
        axis([2015 2200 0 0.6]);

        subplot(4,4,5)% per capita consumption
        if km==1
            plot(tt,percapcons(realtime2,TW),'LineStyle','-','LineWidth',1,'Color',color(ex,1:3)); hold on;
        elseif km==2
            plot(tt,percapcons(realtime2,TW),'LineStyle','--','LineWidth',1,'Color',color(ex,1:3)); hold on;
        end
        realtime3=[121:396];
        if km==1
            for i=1:max(size(realtime3))
                util_0TW(i,ex)=sum(output_var(realtime3(i),TW,15:18),3);
            end
        else
            for i=1:max(size(realtime3))
                util_1TW(i,ex)=sum(output_var(realtime3(i),TW,15:18),3);
            end
        end



    end
end

subplot(4,4,6) %dU 0TW
for i=1:3
    utilplot_1(:,i)=util_0TW(:,i)-util_0TW(:,i+1);
    plot([121:396],utilplot_1(:,i),'LineStyle','-','LineWidth',1,'Color',color(i,1:3)); hold on;
end

subplot(4,4,7) %dU 1TW
for i=1:4
    utilplot_2(:,i)=util_1TW(:,i)-util_0TW(:,i);
    plot([121:396],utilplot_2(:,i),'LineStyle','-','LineWidth',1,'Color',color(i,1:3)); hold on;
end
for km=1:2
    if km==1
        TW=1;%1for 0TW and 11 for 1TW
    else
        TW=11;
    end
    explist=[8:11];
    for ex=1:4
        T=explist(ex);
        load(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_var_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_emi_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
       if km==1
            for i=1:max(size(realtime3))
                util_0TW(i,ex)=sum(output_var(realtime3(i),TW,15:18),3);
            end
        else
            for i=1:max(size(realtime3))
                util_1TW(i,ex)=sum(output_var(realtime3(i),TW,15:18),3);
            end
        end
    end
end
for i=1:4
    utilplot_2(:,i)=util_1TW(:,i)-util_0TW(:,i);
    plot([121:396],utilplot_2(:,i),'LineStyle','--','LineWidth',1,'Color',color(i,1:3)); hold on;
end


for k=1:2
    if k==1
        TW=1;%1for 0TW and 11 for 1TW
        color=gray(6);
    else
        TW=11;
        color=jet(6);
    end

    explist=[1:4];
    for ex=1:4

        T=explist(ex);


        load(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_var_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_emi_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');

        dW=zeros(16,2);
        for s=1:16
            dW(s,1)=dW(s,1)+sum(output_nuclear(1:(cn_num-1),s,2),1); % global Gt CO2 abated by nuclear power
            for cn=1:(cn_num-1)
                dW(s,2)=dW(s,2)+(output_nuclear(cn,s,3)-output_nuclear(cn,1,3))*output_pop(400,cn); % consumption trillion $
            end
        end
        nv(:,ex)=dW(:,2);
        for t=121:396
            cost_0(t,1)=output_var(t,TW,13);% renew E cost
            cost_0(t,2)=(output_var(t,TW,37)+output_var(t,TW,36)-output_var(t-1,TW,36));% warming cost
            cost_0(t,3)=(output_var(t,TW,36));% tipping cost
            cost_0(t,4)=output_var(t,TW,39)./output_var(t,TW,4);% nuc install cost
            cost_0(t,5)=output_var(t,TW,40)./output_var(t,TW,4);% nuc operation cost
            cost_0(t,6)=output_var(t,TW,28)./output_var(t,TW,4);% nuc dam cost
        end
        cost_0(121,4)=cost_0(122,4);% nuc install cost
        if ex==2 && k==1 %2C 0TW
            cost_2C=cost_0;
        elseif ex==4 && k==1 % CUR 0TW
            cost_CUR=cost_0;
        elseif ex==4 && k==2 % CUR 1TW
            cost_CUR_2=cost_0;
        elseif ex==2 && k==2 %2C 1TW
            cost_2C_2=cost_0;
        elseif ex==3 && k==1 % NDC 0TW
            cost_NDC=cost_0;
        end



    end
end
color=[238 181 197
    16 89 210
    34 173 47
    80 80 80
    178 161 199]./255;
for x=1:3
    costplot1(:,x)=cost_CUR(:,x)-cost_2C(:,x);
end
lowline=min(costplot1,[],'all');

DAC2300=zeros(3,5);
for x=1:5
    costplot2(:,x)=cost_2C(:,x+1)-cost_2C_2(:,x+1);
end
lowline2=min(costplot2,[],'all');

for x=1:5
    costplot3(:,x)=cost_CUR(:,x+1)-cost_CUR_2(:,x+1);
end
lowline3=min(costplot3,[],'all');

Fcolor=[180 180 180
    238 181 197
    16 89 210
    34 173 47
    229 223 236]./255;
for p=1:3

    if p==1
        costplot=costplot1;
        knum=[1:3];
    elseif p==2
        costplot=costplot2;
        knum=[1:5];
    else
        costplot=costplot3;
        knum=[1:5];
    end
    lowline=min(costplot,[],'all');
    partplot1=zeros(size(costplot,1),size(costplot,2));
    partplot2=zeros(size(costplot,1),size(costplot,2));
    m=1;
    l=1;
    for k0=1:max(size(knum))
        k=knum(k0);
        for n=1:size(costplot,1)
            if costplot(n,k)<0
                partplot1(n,k0)=costplot(n,k);

            elseif costplot(n,k)>0
                partplot2(n,k0)=costplot(n,k);
            end
        end
    end

    for k2=1:k0
        partplot1_sum(:,k2)=sum(partplot1(:,1:k2),2);
        partplot2_sum(:,k2)=sum(partplot2(:,1:k2),2);
        subplot(4,4,p+7)
        plot([121:396],partplot1_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color',Fcolor(k2,:));hold on
        plot([121:396],partplot2_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color',Fcolor(k2,:));hold on
    end
    subplot(4,4,p+7)
    plot([121:396],sum(partplot1(121:396,:),2)+sum(partplot2(121:396,:),2),'LineStyle','-','LineWidth',1,'Color',[0 0 0]); hold on;
    plot([121:396],partplot1_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color','b'); hold on;
    plot([121:396],partplot2_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color','r'); hold on;
    if p==3
        ylim([-0.005 0.01]);
    elseif p==2
        ylim([-0.005 0.005]);
    end
    for d=1:size(costplot,2)
        DAC2300(p,d)=sum(costplot(:,d))./275.*100;
    end
end

for ex=1:3
for k=1:2
    if k==1
        TW=1;%1for 0TW and 11 for 1TW
    else
        TW=11;
    end
        T=4+ex;
        load(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_var_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
        load(strcat('..\nuclearmonte\output_emi_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');

        dW=zeros(16,2);
        for s=1:16
            dW(s,1)=dW(s,1)+sum(output_nuclear(1:(cn_num-1),s,2),1); % global Gt CO2 abated by nuclear power
            for cn=1:(cn_num-1)
                dW(s,2)=dW(s,2)+(output_nuclear(cn,s,3)-output_nuclear(cn,1,3))*output_pop(400,cn); % consumption trillion $
            end
        end
        nv(:,ex)=dW(:,2);
            for t=121:396
                cost_0(t,1)=output_var(t,TW,13);% renew E cost
                cost_0(t,2)=(output_var(t,TW,37)+output_var(t,TW,36)-output_var(t-1,TW,36));% warming cost
                cost_0(t,3)=(output_var(t,TW,36));% tipping cost
                cost_0(t,4)=output_var(t,TW,39)./output_var(t,TW,4);% nuc install cost
                cost_0(t,5)=output_var(t,TW,40)./output_var(t,TW,4);% nuc operation cost
                cost_0(t,6)=output_var(t,TW,28)./output_var(t,TW,4);% nuc dam cost
            end
               cost_0(121,4)=cost_0(122,4);% nuc install cost
        if k==1
            cost_CUR=cost_0;
        elseif  k==2
            cost_CUR_2=cost_0;
        end
end
color=[238 181 197
    16 89 210
    34 173 47
    178 161 199]./255;


for x=1:5
    costplot3(:,x)=cost_CUR(:,x+1)-cost_CUR_2(:,x+1);
end

Fcolor=[180 180 180
    238 181 197
    16 89 210
    34 173 47
    229 223 236]./255;
p=ex;
        costplot=costplot3;
        knum=[1:5];

    partplot1=zeros(size(costplot,1),size(costplot,2));
    partplot2=zeros(size(costplot,1),size(costplot,2));
    m=1;
    l=1;
    for k0=1:max(size(knum))
        k=knum(k0);
    for n=1:size(costplot,1)
        if costplot(n,k)<0
            partplot1(n,k0)=costplot(n,k);

        elseif costplot(n,k)>0
            partplot2(n,k0)=costplot(n,k);
        end
    end
    end

    for k2=1:k0
        partplot1_sum(:,k2)=sum(partplot1(:,1:k2),2);
        partplot2_sum(:,k2)=sum(partplot2(:,1:k2),2);
        subplot(4,4,p+10)
        plot([121:396],partplot1_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color',Fcolor(k2,:));hold on
        plot([121:396],partplot2_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color',Fcolor(k2,:));hold on
    end
      lowline=min(partplot1_sum,[],'all');
    subplot(4,4,p+10)
    plot([121:396],sum(partplot1(121:396,:),2)+sum(partplot2(121:396,:),2),'LineStyle','-','LineWidth',1,'Color',[0 0 0]); hold on;
    plot([121:396],partplot1_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color','b'); hold on;
    plot([121:396],partplot2_sum(121:396,k2),'LineStyle','-','LineWidth',1,'Color','r'); hold on;
    for d=1:5
    DAC2300_2(ex,d)=sum(costplot(:,d))./275.*100;
    end

end