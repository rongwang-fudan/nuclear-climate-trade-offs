% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2023
tic
clear;

load('files\cndata134.dat','-mat'); cndata=cndata2; clear cndata2;
load('..\nuclear\output_pop.dat','-mat'); 
load('..\nuclear\output_consum.dat','-mat');
cn_num=size(cndata,1);
EndSav=1;
load('files\codelist.dat','-mat');


nvmc=zeros(21,39+39*5+39*39*5);
ans0=zeros(5,4);
exlist=[1:20];

for ex0=1:20
    ex=exlist(ex0);
for mc=1:(39+39*5+39*39*5)

T=(39+39*5+39*39*5)*(ex-1)+mc+1;

load(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');
load(strcat('..\nuclearmonte\output_var_EndSav',num2str(EndSav),'_MC',num2str(T),'.dat'),'-mat');

dW=zeros(16,2);
test=0;
for s=1:16
    dW(s,1)=dW(s,1)+sum(output_nuclear(1:(cn_num-1),s,2),1); % global Gt CO2 abated by nuclear power
    for cn=1:(cn_num-1)
        dW(s,2)=dW(s,2)+(output_nuclear(cn,s,3)-output_nuclear(cn,1,3))/(output_pop(121,cn)/output_consum(cn,1))^1.45*(1000^1.45); % consumption trillion $
    end
end
nvmc(ex0,mc)=dW(11,2);

end   
end

save(strcat('..\nuclear\NEI_sentest.dat'),'nvmc');
