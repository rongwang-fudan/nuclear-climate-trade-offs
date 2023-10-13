% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
tic
clear;

load('files\cndata134.dat','-mat'); cndata=cndata2; clear cndata2;
load('files\cou_iform.dat','-mat'); % 1: id for 222 countries; 2: 2 developing/ 1 developed; 3: 12 region id; 4 OECD; 5 id for 112 countries; 6 pi temperature
cou_iform=cou_iform2; clear cou_iform2;
WGR=load('files\WGR_GHG.txt');
load('..\nuclear\output_pop.dat','-mat'); % output_pop=zeros(400,cn_num-1);
cn_num=size(cndata,1);
EndSav=1;
ncp=[0:0.1:1.5]; ncp(12:16)=[2 3 4 5 6]; ntw=ncp';
% ncp=[0:0.5:7.5]; ncp(12:16)=[5.2 5.4 5.6 5.8 6];ntw=ncp';
load('files\codelist.dat','-mat');
ntw2=[0.5 1 1.5 2 3 4]; ntw2=ntw2';
tip_para=[0.063 50 15 0; 0.188 1500 10 0.067; 0.104 500 5 0.2; 0.163 50 5 1; 0.053 50 10 0.2];
cmap2=jet(21);
cmap3=jet(6);

%FudanCCM
events=zeros(2000,7+5); % 7 level of accident + 5 tipping elements
i2=0;
color=jet(6);
nmc=1;

% nv=zeros(16,14*4);
% sccmc=zeros(15,nmc,4);
% scc=zeros(15,14*4);
% explist=[10 15 20 25 30]+67;

nvmc=zeros(30,39+39*5+39*39*5);
ans0=zeros(5,4);

for ex=1:30
for mc=1:(39+39*5+39*39*5)

T=(39+39*5+39*39*5)*(ex-1)+mc+1;

load(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
load(strcat('..\nuclearmonte\output_var_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
% load(strcat('..\nuclearmonte\output_emi_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');

dW=zeros(16,2);
test=0;
for s=1:16
    dW(s,1)=dW(s,1)+sum(output_nuclear(1:(cn_num-1),s,2),1); % global Gt CO2 abated by nuclear power
    for cn=1:(cn_num-1)
        dW(s,2)=dW(s,2)+(output_nuclear(cn,s,3)-output_nuclear(cn,1,3))*output_pop(400,cn); % consumption trillion $
%       test=test+output_nuclear(cn,11,3)*output_pop(400,cn);
    end
end
nvmc(ex,mc)=dW(11,2);

end   
end

save(strcat('..\nuclear\NEI_sentest.dat'),'nvmc');

