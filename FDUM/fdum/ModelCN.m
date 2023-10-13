% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.31
% Establishing the energy system by zone

function gdpcn = ModelCN( FFlux, tempz, L, calrsav, output_cap, deffs, covidyear, endsav, mcmode )

global elas carbonbudget20 Eland EFco2 Fex cndata deltarf fudanccm_exo cou_iform burkecoef tempzone nuclear_tip nuclear_tech nuclear_acc randnumc invest_trans energy_growth mac_curves neg_year nuclear_D_GDP abty INES_sen NPG_sen sen_open

cn_num=size(cndata,1);

if mcmode==1
    mc2=36;
elseif mcmode>=2 && mcmode<=6
    mc2=1000;
elseif mcmode==7
    mc2=0;
elseif mcmode==8
    mc2=2;
elseif mcmode==9
    xra=[-1.96000000000000;-1.64813271378742;-1.28613372071428;-1.03696234321609;-0.844492884450751;-0.676925320813962;-0.527894367333415;-0.402920558568819;-0.265930233109539;-0.149229844668551;-0.0205257114860915;0.113957817061441;0.242432349790607;0.376949571488937;0.508044017944601;0.658551851844843;0.833469925620126;1.03016998482561;1.28025683193193;1.63438823247971;1.96000000000000];
    mc2=4000;
elseif mcmode==10
    mc2=2;
end


for mc=8:11

    display(mc);

    dpo=2; elas0=0.4; esc2=3.1; nonCO2=1;
    elasmu=1.45; % elasticity of marginal utility of consumption 1.45 (default)
    tempdiff=3; % 1 using the same temperature everywhere; 2 temperature by country, but the same warming; 3 temperature and warming by country (default)
    fudanccm_exo=1; % 1 fudanccm endo (default); 2 fudanccm exo; 3 dice
    weitzman=0; % 0 for normal damage function (default); 1 for weitzman damage function ;2 for Global_function aT^2
    tipping=0; % 0 for no tipping; 1 for tipping without interaction; 2 for tipping with interaction (default)
    neg=0; % negative emissions Gt CO2 y-1
    neg_year=1970; % default year to deploy negative emissions
    rour=15; % social time preference of consumption 15 (default)
    LR=0.2; % learning rate 0.2 (default)
    tempopt=13; % optimal temperature for productivity 13 (default)
    gc=[0.014 0.04]; % initital per capita GDP growth rate for rich / poor countries
    qlen=10; % length of transition 10 (default)
    dicetfp=0; % 0 for exogeneous; 1 for Burke2015 with free TFP; 2 for Burke2015 (default); 3 for Burke2015 constrained + learning
    nuclear_tip=0; % 0 for no nuclear; 1 for nuclear tipping; 2 for stochastic nuclear tipping
    nuclear_tech=1; % 0 no tech advances; 1 national tech advances; 2 global cooperation; 3 OECD/REF/ALM/ASIA cooperation; 4 no accident; 5 low nuclear costs ($2/We); 6 lifetime 20 yrs
    nuclear_D_GDP=1;%1 damage use normal GDP 2 damage use 2020 GDP
    nuclear_acc=5; % 1 damage as a % of GDP; 2 damage as a % of GDP forced (calibrated to level 7); 3 damage as absolute $; 4 damage as absolute $ (calibrated to level 7)
    burkecoef=[0.0005 0.0005]; % coefficient for poor and rich countries 0.0005 (default)
    invest_trans=[0 0 0 0 0 0 ]; % 1 fraction of investment transfer; 2 threshold of per capita GDP k$; 3 fraction of investment; 4 block(1)/open(2) transfer from region 5 to region 6
    energy_growth=0; % 0 without (default) / 1 with fast energy production growth in developing countries
    mac_curves=0; % 0 using the global MACC (default); 1 using the MACC by project
    techdiff=0; % 0 for technological diffusion by country (default); 1 for global technological diffusion; 2 for regional technological diffusion; 3 for no technological diffusion
    aff_pop=0.1;%Percentage of affected population
    lim_t=1.5;%tem of limit for tip
    CONSTR_nuc=10;
    load('files\nuccnlist.dat','-mat');
    load('files\GDPcorrect_Pforacc.dat','-mat');
    load('files\popcorrect_Pforacc.dat','-mat');
    abty=2020;
INES_sen=0;  NPG_sen=0;
sen_open=0;
    unc=ones(1,19);
    % 1 for sens; 2 MC-FUDAM; 3 MC-FUDAM (nonCO2); 4 MC-FUDAM (nonCO2+neg); 5 MC-DICE; 6 MC-DICE (burke); 7 SCC; 8 SCC-nuclear; 9 Stochastic-SCC-nuclear
    if mcmode==1
        if mc==1
            gc(2)=gc(2)-0.01;
        elseif mc==2
            gc(2)=gc(2)+0.01;
        elseif mc==3
            esc2=2.6;
        elseif mc==4
            esc2=3.9;
        elseif mc==5
            LR=0.1;
        elseif mc==6
            LR=0.3;
        elseif mc==7
            tempopt=12;
        elseif mc==8
            tempopt=14;
        elseif mc==9
            rour=10;
        elseif mc==10
            rour=20;
        elseif mc==11
            elasmu=0.9;
        elseif mc==12
            nonCO2=0.5;
        elseif mc==13
            weitzman=1;
        elseif mc==14
            tipping=2;
        elseif mc==15
            neg=10;
        elseif mc>=16
            fudanccm_exo=3;
            if mc==17
                rour=10;
            elseif mc==18
                rour=20;
            elseif mc==19
                elasmu=0.9;
            elseif mc==20
                weitzman=1;
            elseif mc==21
                tipping=2;
            elseif mc==22
                dicetfp=3;
            elseif mc==23
                dicetfp=2;
            elseif mc==24
                dicetfp=2; rour=10;
            elseif mc==25
                dicetfp=2; elasmu=0.9;
            elseif mc==26
                dicetfp=2; weitzman=1;
            elseif mc==27
                dicetfp=2; tipping=2;
            elseif mc==28
                dicetfp=2; neg=10;
            elseif mc==29
                neg=10;
            elseif mc==30
                nonCO2=0.5;
            elseif mc==31
                dicetfp=2; rour=20;
            elseif mc==32
                dicetfp=2; nonCO2=0.5;
            elseif mc==33
                tipping=2; neg=10;
            elseif mc==34
                fudanccm_exo=1; tipping=2; neg=10;
            elseif mc==35
                dicetfp=2; neg=10; tipping=2;
            elseif mc==36
                fudanccm_exo=1; neg=10; nonCO2=0.5;
            end
        end
    elseif mcmode>=2 && mcmode<=6
        if mc>0
            esc2 = min(3.9,max(2.6,3.1 + 1.3/4 * randn)); % (3.9-2.6)/2/2=0.32 as sigma
            LR = 0.1+(0.3-0.1)*rand; % uniform distribution
            unc(16) = min(2,max(0.2,1+randn)); % uncertainty for the slope of EUE to omega
            unc(17) = min(2,max(0.2,1+randn)); % uncertainty for the slope of ENE to omega
            tempopt = min(14,max(12,13 + 0.5 * randn)); % +/-0.5 as sigma
            if mcmode==3
                nonCO2 = 0.5;
            elseif mcmode==4
                nonCO2 = 0.5; neg=10;
            elseif mcmode==5
                fudanccm_exo=3;
            elseif mcmode==6
                fudanccm_exo=3; dicetfp=2;
                burkecoef(1)=burkecoef(1)*(1+0.8 * randn); % 80% standard error for poor countries
                burkecoef(2)=burkecoef(2)*(1+0.2 * randn); % 20% standard error for poor countries
            end
        end
    elseif mcmode==7 % scc

    elseif mcmode==8 % scc-nuclear-average
        tipping=2;
        neg=14.3; % negative emissions Gt CO2 y-1
        neg_year=2050;
        nuclear_tip=1;  % 0 for no nuclear; 1 for nuclear tipping; 2 for stochastic nuclear tipping
        nuclear_tech=1; % 0 no tech advances; 1 national tech advances; 2 global cooperation; 3 OECD/REF/ALM/ASIA cooperation; 4 no accident; 5 low nuclear costs ($2/We); 6 lifetime 20 yrs
        nuclear_acc=4;  % 1 damage as a % of GDP; 2 damage as a % of GDP forced (calibrated to level 7); 3 damage as absolute $; 4 damage as absolute $ (calibrated to level 7)
        randnumc=randnu(1,1:3290);
        randnumc(1,1:5)=0.5;
        if mc==0
            nzlen=20;
        elseif mc==1
            nzlen=10;
        elseif mc==2
            nzlen=20;
            nuclear_tech=4; % no nuclear accident
        elseif mc==3
            nzlen=20;
            tipping=0; % no tipping
        elseif mc==4
            nzlen=20;
            fudanccm_exo=3; dicetfp=0; % DICE
        elseif mc==5
            nzlen=20;
            fudanccm_exo=3; dicetfp=2; % DICE-Burke
        end
    elseif mcmode==9 % scc-nuclear-stochastic
        tipping=2;
        neg=10; % negative emissions Gt CO2 y-1
        neg_year=2050;
        nuclear_tip=2; % 0 for no nuclear; 1 for nuclear tipping; 2 for stochastic nuclear tipping
        nuclear_tech=1; % 0 no tech advances; 1 national tech advances; 2 global cooperation; 3 OECD/REF/ALM/ASIA cooperation; 4 no accident; 5 low nuclear costs ($2/We); 6 lifetime 20 yrs
        nuclear_acc=5; % 1 damage as a % of GDP; 2 damage as a % of GDP forced (calibrated to level 7); 3 damage as absolute $; 4 damage as absolute $ (calibrated to level 7)
        nzlenlist=[7 14 23 33];
        if mc==1
            nzlen=nzlenlist(mc); % 1.5C pathway
            nuclear_tech=1;
            esc2=3;nonCO2=0.5;rour=20;
            neg=15;
        elseif mc==2
            nzlen=nzlenlist(mc); % 2C pathway
            nuclear_tech=1;
            esc2=3;nonCO2=0.5;rour=20;
            neg=5;
        elseif mc==3
            nzlen=nzlenlist(mc); % NDC
            nuclear_tech=1; %
            esc2=3;nonCO2=0.6;rour=20;neg=0;
        elseif mc==4
            nzlen=nzlenlist(mc); % CUR
            nuclear_tech=1;
            esc2=3;nonCO2=0.7;rour=20;neg=0;
        elseif mc==5
            nzlen=nzlenlist(4); % INSE rou test
            nuclear_tech=1; 
            esc2=3;nonCO2=0.7;rour=20;neg=0;
            INES_sen=1;  
        elseif mc==6
            nzlen=nzlenlist(4); % construction year test
            nuclear_tech=1;
            esc2=3;nonCO2=0.7;rour=20;neg=0;
            NPG_sen=1;
        elseif mc==7
            nzlen=nzlenlist(4); % tipping test
            nuclear_tech=1;
            esc2=3;nonCO2=0.7;rour=20;neg=0;
            sen_open=1;
        elseif mc==8
            nzlen=nzlenlist(1); % 1.5C pathway & no NPG improvement
            nuclear_tech=1;
            esc2=3;nonCO2=0.5;rour=20;
            neg=15;
            NPG_sen=2;
        elseif mc==9
            nzlen=nzlenlist(2); % 2C pathway & no NPG improvement
            nuclear_tech=1;
            esc2=3;nonCO2=0.5;rour=20;
            neg=5;
            NPG_sen=2;
        elseif mc==10
            nzlen=nzlenlist(3); % NDC & no NPG improvement
            nuclear_tech=1; %
            esc2=3;nonCO2=0.6;rour=20;neg=0;
            NPG_sen=2;
        elseif mc==11
            nzlen=nzlenlist(4); % CUR & no NPG improvement
            nuclear_tech=1;
            esc2=3;nonCO2=0.7;rour=20;neg=0;
            NPG_sen=2;
        end
    elseif mcmode==10
        LR=0.2;
        techdiff=0; % 0 for technological diffusion by country (default); 1 for global technological diffusion; 2 for regional technological diffusion; 3 for no technological diffusion
        % 1 fraction of investment transfer completed; 2 threshold of per capita GDP k$; 3 fraction of investment
        invest_trans=[0 10 0.05 0 0 0]; % 1 fraction of investment transfer; 2 threshold of per capita GDP k$; 3 fraction of investment; 4 block(1)/open(2) transfer from region 5 to region 6
        energy_growth=10; % final per capita energy MWh
        mac_curves=1; % 0 using the global MACC (default); 1 using the MACC by project
    end

    %Regression of regional temperature with global temperature
    tempzone=zeros(11,2);
    for regionid=1:11
        tempzone(regionid,1) = max(0,min(2.5,tempz(regionid,1) + tempz(regionid,2) * (unc(19)-1))); % slope of the varied curve
        tempzone(regionid,2) = tempz(regionid,3) - tempz(regionid,4) * tempzone(regionid,1); % intercept of the varied curve
    end

    %Equilibrium sensitivity of climate
    FFlux0=FFlux;
    FFlux0(1) = esc2/deltarf; % K / (W m2)

    %Non-CO2 GHG
    AerosolsLUC;
    Fex0=Fex*nonCO2;

    %Calibration of climate damage function
    [dcoef, xy_damage] = damage( dpo, 2 );

    %Parameters
    EFco20=EFco2;
    carbonbudget20s=carbonbudget20;
    Eland0 = Eland;
    dcoef0=dcoef;
    LR0=LR;

    %Coefficient in the damage function
    dcoef = dcoef0 * unc(1);
    %Elasticity of substitution (avoid a zero-like elas)
    elas = elas0 * unc(15);
    %Learning rate on the cost curve
    LR = LR0 * unc(18);
    %Equilibrium sensitivity of climate
    FFlux(1) = FFlux0(1) * unc(2);
    %Time inertia of climate system to reach equilibirum (year)
    FFlux(2) = FFlux0(2) * unc(3);
    %air to land biosphere GtC/yr
    FFlux(4) = FFlux0(4) * unc(6);
    %air to surface ocean GtC/yr
    FFlux(5) = FFlux0(5) * unc(7);
    %land biosphere to soil GtC/yr
    FFlux(6) = FFlux0(6) * unc(8);
    %surface soil to deep soil GtC/yr
    FFlux(7) = FFlux0(7) / unc(9);
    %surface ocean to deep ocean GtC/yr
    FFlux(8) = FFlux0(8) / unc(10);
    %CO2 emission factors for fossil fuel only tCO2 / MJ
    EFco2 = EFco20 * unc(4);
    carbonbudget20(:,2) = carbonbudget20s(:,2) * unc(4);
    %CO2 emissions from land use change
    Eland = Eland0 * unc(5);
    %Radiative forcing by 1 CH4, 2 N2O, 3 CFCs, 4 aerosol
    for i=1:4
        Fex(:,i) = Fex0(:,i) * unc(i+10);
    end

    %Slope of induced efficiency changes
    iec_cn = Calibration_IEC( cndata, 0 );
    for cn=1:cn_num
        iec_cn(cn,2) = iec_cn(cn,2) * (1 + (unc(16)-1)*iec_cn(cn,3));
        iec_cn(cn,4) = iec_cn(cn,4) * (1 + (unc(16)-1)*iec_cn(cn,5));
        iec_cn(cn,12) = iec_cn(cn,12) * (1 + (unc(17)-1)*iec_cn(cn,13));
        iec_cn(cn,14) = iec_cn(cn,14) * (1 + (unc(17)-1)*iec_cn(cn,15));
    end

    parameters=[esc2, LR, tempopt, iec_cn(25,2), iec_cn(25,12), iec_cn(25,14)];

    %Switcher for  C1	C2	S1	S2	S3	S4	S5	T1	T2	T3	T4
    switcher = ones(1,12);

    %Input of data by country: cndata(1+119,1+49*4) KEYL+Et+Ez
    [econcn] = countrydata( output_cap, cndata );
    % save('files\econcn134.dat','econcn'); % econcn=zeros(49,23,cn_num)

    %Exogeneous FudanCCM
    if fudanccm_exo==2
        switcher(2)=0; % deactivating the impact of climate change on efficieney of ENE
        switcher(4)=0; % deactivating the impact of climate change on efficieney of EUE
    end

    % Abatement by country
    % 1 start of mitigation; 2 final share of clean energy (%); 3 length of transition; 4 growth rate of per capita GDP in 2025
    % 5 fraction of emission abatements (0 for net-zero emissions; >0 for additional negative emissions); 6 national negative emissions (Gt CO2)
    % 7 start of nuclear transition; 8 global new nuclear power TW; 9 length of nuclear transition

    abtcn=zeros(cn_num,9);
    abtcn(1:cn_num,1)=abty;
    abtcn(1:cn_num,2)=1;
    abtcn(1:cn_num,3)=qlen;
    for cn=1:cn_num
        if iec_cn(cn,29)==1
            abtcn(cn,4)=gc(2); % poor countries
        else
            abtcn(cn,4)=gc(1); % rich countries
        end
    end
    tt=[0:299]; ett=exp(-tt/50);
    optlrsav=zeros(300,1); % 1 for 2025
    scc=[2025 0]; % 1 year of scc; 2 estimate scc by adding a pulse emissions

    % Output variables
    output_var=zeros(400,16,40);
    % 1 warming C, 2 energy PWh, 3 net-zero energy PWh, 4 net gdp t$, 5 emi Gt CO2,
    % 6 saving rate, 7 eue, 8 ene, 9 capital non-energy, 10 capital energy,
    % 11 capital green energy, 12 energy cost share, 13 cost of mitigation as
    % a percentage of GDP, 14 carbon price, 15-18 util for OECD REF ASIA and ALM
    % 19-23 tipping case AMOC GIS WAIS AMAZ and ENSO; 24 total gdp t$
    % 25 nuclear capacity factor, 26 nuclear power TW, 27 % nuclear abated emissions Gt CO2
    % 28 nuclear accident damage t$, 29 nuclear cost t$, 30 number of nuclear accidents
    % 31-36 damage of tipping AMOC GIS WAIS AMAZ, ENSO and total as % to GDP
    % 37 damage of climate warming as % to GDP
    % 39 nuclear cost of install 40 nuclear cost of operation 

    output_util=zeros(cn_num-1,16,50);
    output_util2=zeros(cn_num-1,16,50);
    output_scc=zeros(cn_num,16,50);
    output_pop=zeros(400,cn_num-1);
    output_gdp=zeros(400,cn_num-1,3);
    output_cost=zeros(276,4*3,cn_num-1);
    output_temp=zeros(36,16); % warming in 2025, 2030, ... 2200
    output_nuclear=zeros(cn_num-1,16,3+100*3); % 1 cumul energy PWh; 2 cumul abated Gt CO2; 3 utility 2025-2300; 4-6 nuclear accident 1 (time, level, damage); 7-9 nuclear accident 2
    output_invest=zeros(76,cn_num-1,16*2); % 1 climate investment contributed by a country; 2 climate investment used in this country


    for s=16:-1:1
        abtcn(1:cn_num,1)=2020+s*5;
        % Saving rate
        if endsav==1
            optlrsav(1:300,1)=0.258; % default value in DICE
        elseif endsav==2
            optlrsav=0.26+0.001*(ett'-1); % endogeneous saving rate
            if s<16 && mcmode==1
                optlrsav(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
            end
        end

        % Nuclear energy
        if mcmode==8 || mcmode==9
            ncp=[0:0.1:1.5]; ncp(12:16)=[2 3 4 5 6]; 
            abtcn(1:cn_num,1)=abty;
            abtcn(1:cn_num,2)=1;
            abtcn(1:cn_num,3)=nzlen;
            abtcn(1:cn_num,7)=2025;   % the time starting to deploy nuclear power
            abtcn(1:cn_num,8)=ncp(s); % global nuclear power TW
            abtcn(1:cn_num,9)=10;%CONSTR_nuc; % the time needed to construct a nuclear power plant
        end

        % Investment transfers
        if mcmode==10
            abtcn(1:cn_num,1)=abty;
            abtcn(1:cn_num,2)=1;
            abtcn(1:cn_num,3)=200;
            invest_trans(1)=1;
            invest_trans(2)=20;
            if s<=11
                invest_trans(3)=(s-1)/2; % fraction of investment
            else
                continue;
            end
        end

        % Simulation from 2025 to 2300
        [npv, tip, S, emico2_0,costadd] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav, elasmu, rour, scc,  lim_t, nuccnlist,CONSTR_nuc,GDPcorrect_Pforacc,popcorrect_Pforacc);

        %GDP by country
        gdpcn=zeros(136,2);
        gdpcn(1:136,1)=S(49,7,1:136);
        gdpcn(1:136,2)=S(121+2200-2025,7,1:136);
        output_pop(1:400,1:(cn_num-1))=S(1:400,35,3:(cn_num+1));

        % Monte Carlo for MIC
        if mcmode>=2 && mcmode<=6
            display(npv);
            for ri=1:50
                for cn=1:(cn_num-1)
                    for ti=121:396
                        output_util(cn,s,ri) = output_util(cn,s,ri) + S(ti,36,cn+2) * ((1-ri*0.001)/(1-rour*0.001))^(ti-121);
                    end
                end
            end
            output_temp(1:36,s)=S(121:5:296,31,1); % warming in 2025, 2030, ... 2200
            continue;
        end

        % Endogeneous saving rate
        if endsav==2
            % Optimizing step 1
            adj_a=0;
            adj_b=0;
            optlrsav2=0.26-0.01+0.001*(ett'-1);
            if s<16
                optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
            end
            [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
            if npv2>=npv
                adj_a=1;
                npv=npv2;
                S=S2;
                tip=tip2;
                optlrsav=optlrsav2;
            else
                optlrsav2=0.26+0.01+0.001*(ett'-1);
                if s<16
                    optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
                end
                [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
                if npv2>=npv
                    adj_a=-1;
                    npv=npv2;
                    S=S2;
                    tip=tip2;
                    optlrsav=optlrsav2;
                end
            end
            if adj_a~=0
                for i_sav=2:10
                    optlrsav2=0.26-0.01*adj_a*i_sav+0.001*(ett'-1);
                    if s<16
                        optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
                    end
                    [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
                    if npv2>=npv
                        npv=npv2;
                        S=S2;
                        tip=tip2;
                        optlrsav=optlrsav2;
                        adj_b=i_sav;
                    else
                        break;
                    end
                end
            end
            % Optimizing step 2
            adj_c=0;
            adj_d=0;
            optlrsav2=0.26-0.01*adj_a*adj_b-0.001+0.001*(ett'-1);
            if s<16
                optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
            end
            [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
            if npv2>=npv
                adj_c=1;
                npv=npv2;
                S=S2;
                tip=tip2;
                optlrsav=optlrsav2;
            else
                optlrsav2=0.26-0.01*adj_a*adj_b+0.001+0.001*(ett'-1);
                if s<16
                    optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
                end
                [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
                if npv2>=npv
                    adj_c=-1;
                    npv=npv2;
                    S=S2;
                    tip=tip2;
                    optlrsav=optlrsav2;
                end
            end
            if adj_c~=0
                for i_sav=2:9
                    optlrsav2=0.26-0.01*adj_a*adj_b-0.001*adj_c*i_sav+0.001*(ett'-1);
                    if s<16
                        optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
                    end
                    [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
                    if npv2>=npv
                        npv=npv2;
                        S=S2;
                        tip=tip2;
                        optlrsav=optlrsav2;
                        adj_d=i_sav;
                    else
                        break;
                    end
                end
            end
            % Optimizing step 3
            adj_e=0;
            for j_sav=1:10
                optlrsav2=0.26-0.01*adj_a*adj_b-0.001*adj_c*adj_d+(0.001+0.01*j_sav)*(ett'-1);
                if s<16
                    optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
                end
                [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
                if npv2>=npv
                    npv=npv2;
                    S=S2;
                    tip=tip2;
                    optlrsav=optlrsav2;
                    adj_e=j_sav;
                else
                    break;
                end
            end
            % Optimizing step 4
            for j_sav=1:9
                optlrsav2=0.26-0.01*adj_a*adj_b-0.001*adj_c*adj_d+(0.001+0.01*adj_e+0.001*j_sav)*(ett'-1);
                if s<16
                    optlrsav2(1:(s-1)*5,1)=output_var(121:(120+(s-1)*5),16,6); % baseline
                end
                [npv2, tip2, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav2, elasmu, rour, scc, aff_pop );
                if npv2>=npv
                    npv=npv2;
                    S=S2;
                    tip=tip2;
                    optlrsav=optlrsav2;
                else
                    break;
                end
            end
        end

        % Output variables
        %     output_emi(1:400,:)=S(1:400,20,:);
        output_emi(:,s)=emico2_0(1:400);
        output_var(1:400,s,1)=S(1:400,31,1);
        output_var(1:400,s,2)=sum(S(1:400,12,3:(cn_num+1)),3);
        output_var(1:400,s,3)=sum(S(1:400,21,3:(cn_num+1)),3);
        output_var(1:400,s,4)=sum(S(1:400,7,3:(cn_num+1)),3);
        output_var(1:400,s,5)=sum(S(1:400,20,3:(cn_num+1)),3);
        output_var(1:32,s,6)=calrsav(3);
        output_var(33:37,s,6)=calrsav(4);
        output_var(38:120,s,6)=calrsav(5);
        output_var(121:400,s,6)=optlrsav(1:280,1);
        output_var(1:400,s,7)=S(1:400,1,26);
        output_var(1:400,s,8)=S(1:400,3,26);
        output_var(1:400,s,10)=sum(S(1:400,11,3:(cn_num+1)),3);
        output_var(1:400,s,9)=sum(S(1:400,10,3:(cn_num+1)),3)-output_var(1:400,s,10);
        output_var(1:400,s,11)=sum(S(1:400,17,3:(cn_num+1)),3);
        output_var(1:400,s,14)=mean(S(1:400,19,3:(cn_num+1)),3);

        output_var(1:400,s,39)=sum(costadd(1:400,1,3:(cn_num+1)),3);
        output_var(1:400,s,40)=sum(costadd(1:400,2,3:(cn_num+1)),3);
        for tipp_i=1:5
            output_var(1:400,s,tipp_i+18)=tip(1:400,tipp_i);
            output_var(1:400,s,tipp_i+30)=tip(1:400,tipp_i+5);
        end
        for cn=1:(cn_num-1)
            zz=cou_iform(cndata(cn+1,1),4);
            output_var(1:400,s,zz+14)=output_var(1:400,s,zz+14)+S(1:400,36,cn+2); % Sum of NPV of utility in 2025
            output_var(121:396,s,12)=output_var(121:396,s,12)+S(121:396,13,cn+2).*S(121:396,15,cn+2); % energy cost
            output_var(121:396,s,13)=output_var(121:396,s,13)+S(121:396,13,cn+2).*(S(121:396,5,cn+2)-S(121:396,43,cn+2)./S(120:395,13,cn+2));  % mitigation cost
            output_var(121:396,s,36)=output_var(121:396,s,36)+S(121:396,13,cn+2).*tip(121:396,11);    % tipping damage
            output_var(121:396,s,37)=output_var(121:396,s,37)+S(121:396,13,cn+2).*(S(121:396,6,cn+2)-S(121:396,42,cn+2)-tip(121:396,11));  % warming damage
            output_var(121:396,s,24)=output_var(121:396,s,24)+S(121:396,13,cn+2);

        end
        output_var(121:396,s,12)=output_var(121:396,s,12)./output_var(121:396,s,24); % energy cost share
        output_var(121:396,s,13)=output_var(121:396,s,13)./output_var(121:396,s,24); % mitigation cost as % to GDP
        output_var(121:396,s,36)=output_var(121:396,s,36)./output_var(121:396,s,24); % tipping damage as % to GDP
        output_var(121:396,s,37)=output_var(121:396,s,37)./output_var(121:396,s,24); % warming damage as % to GDP
        if s==16
            output_gdp(1:400,1:(cn_num-1),1)=S(1:400,13,3:(cn_num+1));
            output_gdp(1:400,1:(cn_num-1),2)=S(1:400,15,3:(cn_num+1));
            output_gdp(1:400,1:(cn_num-1),3)=S(1:400,12,3:(cn_num+1));
        end

        % NPV of utility in 2025
        for cn=1:(cn_num-1)
            for ri=1:50
                for ti=121:396
                    output_util(cn,s,ri) = output_util(cn,s,ri) + S(ti,36,cn+2) * ((1-ri*0.001)/(1-rour*0.001))^(ti-121);
                end
            end
            for i=1:276
                if mod(s-1,6)==0
                    si=floor((s-1)/6); % for 2025 2055 2085
                    output_cost(i,si*4+1,cn)=S(i+120,13,cn+2)*S(i+120,5,cn+2); % abatement cost
                    output_cost(i,si*4+2,cn)=S(i+120,13,cn+2)*S(i+120,6,cn+2); % climate cost
                    output_cost(i,si*4+3,cn)=S(i+120,13,cn+2)*(1-S(i+120,5,cn+2)-S(i+120,6,cn+2))*S(i+120,15,cn+2); % energy cost
                    output_cost(i,si*4+4,cn)=S(i+120,13,cn+2)*(1-S(i+120,5,cn+2)-S(i+120,6,cn+2))*(1-S(i+120,15,cn+2)); % non-energy cost
                end
            end
        end

        % Climate Investment Transfer
        if mcmode==10
            for cn=1:(cn_num-1)
                output_invest(1:(invest_trans(4)-2025+1),cn,s)=S(121:(invest_trans(4)-2025+121),45,cn+2); % climate investment contributed by a country
                output_invest(1:(invest_trans(4)-2025+1),cn,s+16)=S(121:(invest_trans(4)-2025+121),46,cn+2); % climate investment used in this country
            end
        end

        % Social cost of carbon
        if mcmode==7
            scc2=[2025 0.1];
            [npv, tip, S2] = Abatementcn( FFlux, L, dpo, dcoef, tempopt, LR, covidyear, iec_cn, calrsav, deffs, techdiff, switcher, econcn, abtcn, tempdiff, weitzman, dicetfp, tipping, neg, optlrsav, elasmu, rour, scc2, aff_pop );
            for cn=1:(cn_num-1)
                % sum of NPV of utility in 2025
                for ri=1:50
                    for ti=121:396
                        output_util2(cn,s,ri) = output_util2(cn,s,ri) + (S(ti,36,cn+2) - S2(ti,36,cn+2)) * ((1-ri*0.001)/(1-rour*0.001))^(ti-121);
                    end
                end
            end
        end

        % nuclear
        if mcmode==8 || mcmode==9
            for cn=1:(cn_num-1)
                output_var(121:396,s,25)=output_var(121:396,s,25)+S(121:396,37,cn+2).*S(121:396,40,cn+2); % nuclear capacity factor
                output_var(121:396,s,26)=output_var(121:396,s,26)+S(121:396,37,cn+2); % nuclear power TW
                output_var(121:396,s,27)=output_var(121:396,s,27)+S(121:396,39,cn+2); % nuclear abated emissions Gt CO2
                output_var(121:396,s,28)=output_var(121:396,s,28)+S(121:396,42,cn+2).*S(120:395,13,cn+2); % nuclear damage trillion $
                output_var(121:396,s,29)=output_var(121:396,s,29)+S(121:396,43,cn+2); % nuclear cost trillion $
                output_nuclear(cn,s,1)=sum(S(121:396,38,cn+2),1); % nuclear energy PWh
                output_nuclear(cn,s,2)=sum(S(121:396,39,cn+2),1); % nuclear abated CO2 emissions Gt CO2
                output_nuclear(cn,s,3)=sum(S(121:396,36,cn+2),1); % nuclear utility discounted to 2025
                idx_nuclear=find(S(121:396,41,cn+2)>0);
                if size(idx_nuclear,1)>0
                    for acci=1:min(100,size(idx_nuclear,1))
                        output_nuclear(cn,s,acci*3+1)=idx_nuclear(acci)+2024; % year of accident
                        output_nuclear(cn,s,acci*3+2)=S(120+idx_nuclear(acci),41,cn+2)+S(120+idx_nuclear(acci),44,cn+2)/10000; % number of accident
                        output_nuclear(cn,s,acci*3+3)=S(120+idx_nuclear(acci),42,cn+2)*S(119+idx_nuclear(acci),13,cn+2); % damage of accident trillion $
                    end
                end
            end
            for tttt=121:396
                if output_var(tttt,s,26)>0
                    output_var(tttt,s,25)=output_var(tttt,s,25)./output_var(tttt,s,26); % nuclear capacity factor
                end
                output_var(tttt,s,30)=sum(S(tttt+1,41,3:(cn_num+1)),3); % number of nuclear accidents
                output_var(tttt,s,38)=sum(S(tttt+1,44,3:(cn_num+1)),3); % number of nuclear accidents level 5-7
            end
        end
    end

    %MIC
    t=[5:5:25];
    zoneid8=[2 2 2 6 5 3 7 1 1 8 4 1];
    % 1. Total Eastern and Southern Africa; 2. Total Northern Africa; 3. Total Western and Central Africa;
    % 4. Total East Asia; 5. Total South and South-east Asia; 6. Total Western and Central Asia; 7. Total Europe;
    % 8. Total Caribbean; 9. Total Central America; 10. Total North America; 11. Total Oceania; 12. Total South America.
    unit1=((sum(S(49,7,3:(cn_num+1)),3)/sum(S(49,35,3:(cn_num+1)),3)*1000)^(1-elasmu)-(sum(S(1,7,3:(cn_num+1)),3)/sum(S(1,35,3:(cn_num+1)),3)*1000)^(1-elasmu))/(1-elasmu);

    if mcmode>=2 && mcmode<=6
        output_mic=zeros(cn_num,12,50);
        for ri=1:50
            for s=1:12
                ms0=0;
                for cn=1:(cn_num-1)
                    A=zeros(1,5);
                    A(1:5)=output_util(cn,s:(s+4),ri);
                    [rs,ms,bs] = regression(t,A);
                    output_mic(cn,s,ri)=-ms/(1-ri*0.001)^(s*5+10)/(S(49,35,cn+2)/1000*unit1); % discounting to the mitigation year (i.e. 2035)
                    ms0=ms0+ms;
                end
                output_mic(cn_num,s,ri)=-ms0/(1-ri*0.001)^(s*5+10)/(sum(S(49,35,3:(cn_num+1)),3)/1000*unit1);
            end
        end
    elseif mcmode==1
        output_mic=zeros(cn_num+8,12); % 8 for eight regions in fig 4a
        for s=1:12
            totalmic=zeros(9,2);
            for cn=1:(cn_num-1)
                A=zeros(1,5);
                A(1:5)=output_util(cn,s:(s+4),rour);
                [rs,ms,bs] = regression(t,A);
                output_mic(cn,s)=-ms/(1-rour*0.001)^(s*5+10)/(S(49,35,cn+2)/1000*unit1); % discounting to the mitigation year (i.e. 2035)
                zone12=zoneid8(cou_iform(cndata(cn+1,1),3));
                totalmic(1,1)=totalmic(1,1)+ms;
                totalmic(1,2)=totalmic(1,2)+S(49,35,cn+2);
                totalmic(zone12+1,1)=totalmic(zone12+1,1)+ms;
                totalmic(zone12+1,2)=totalmic(zone12+1,2)+S(49,35,cn+2);
            end
            for zz9=1:9
                output_mic(cn_num+zz9-1,s)=-totalmic(zz9,1)/(1-rour*0.001)^(s*5+10)/(totalmic(zz9,2)/1000*unit1);
            end
        end
        % mic 50
        output_mic50=zeros(cn_num+8,12,50);
        for ri=1:50
            for s=1:12
                totalmic=zeros(9,2);
                for cn=1:(cn_num-1)
                    A=zeros(1,5);
                    A(1:5)=output_util(cn,s:(s+4),ri);
                    [rs,ms,bs] = regression(t,A);
                    output_mic50(cn,s,ri)=-ms/(1-ri*0.001)^(s*5+10)/(S(49,35,cn+2)/1000*unit1); % discounting to the mitigation year (i.e. 2035)
                    zone12=zoneid8(cou_iform(cndata(cn+1,1),3));
                    totalmic(1,1)=totalmic(1,1)+ms;
                    totalmic(1,2)=totalmic(1,2)+S(49,35,cn+2);
                    totalmic(zone12+1,1)=totalmic(zone12+1,1)+ms;
                    totalmic(zone12+1,2)=totalmic(zone12+1,2)+S(49,35,cn+2);
                end
                for zz9=1:9
                    output_mic50(cn_num+zz9-1,s,ri)=-totalmic(zz9,1)/(1-ri*0.001)^(s*5+10)/(totalmic(zz9,2)/1000*unit1);
                end
            end
        end
    end

    % Social cost of carbon
    if mcmode==7
        demi=0.1/1000; % trillion tonne
        for s=1:16
            for cn=1:(cn_num-1)
                for ri=1:50
                    output_scc(cn,s,ri)=output_util2(cn,s,ri) * ((1-optlrsav(1,1))*S(121,7,cn+2)*1000/S(121,35,cn+2))^elasmu /demi; % scc $/tonne CO2
                    output_scc(cn_num,s,ri)=output_scc(cn_num,s,ri)+output_scc(cn,s,ri);
                end
            end
        end
    end

    % Social cost of carbon for nuclear power
    if mcmode==8 || mcmode==9
        for s=1:15
            demi=(sum(output_nuclear(1:(cn_num-1),s+1,2),1)-sum(output_nuclear(1:(cn_num-1),s,2),1))/1000; % Gt CO2 abated by nuclear power
            for cn=1:(cn_num-1)
                output_pop(400,cn)=((1-optlrsav(1,1))*S(121,7,cn+2)*1000/S(121,35,cn+2))^elasmu;
                for ri=1:50
                    output_scc(cn,s,ri)=(output_util(cn,s+1,ri)-output_util(cn,s,ri)) * output_pop(400,cn) /demi; % scc as USD/tonne CO2 (positive for higher benefits of abated CO2)
                    output_scc(cn_num,s,ri)=output_scc(cn_num,s,ri)+output_scc(cn,s,ri);
                end
            end
        end
    end

    if mcmode==1
        save(strcat('..\mic\output_mic_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_mic');
        save(strcat('..\mic\output_mic50_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_mic50');
        save(strcat('..\mic\output_var_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_var');
        save(strcat('..\mic\output_gdp_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_gdp');
        save(strcat('..\mic\output_pop.dat'),'output_pop');
        if mc==0 || mc==16 || mc==23
            save(strcat('..\mic\output_cost_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_cost');
        end
    elseif mcmode>=2 && mcmode<=6
        save(strcat('..\micmonte\output_mic_sce',num2str(mcmode),'_MC',num2str(mc),'.dat'),'output_mic');
        save(strcat('..\micmonte\output_temp_sce',num2str(mcmode),'_MC',num2str(mc),'.dat'),'output_temp');
        save(strcat('..\micmonte\parameters_sce',num2str(mcmode),'_MC',num2str(mc),'.dat'),'parameters');
        if mc==0
            save(strcat('..\mic\output_pop.dat'),'output_pop');
        end
    elseif mcmode==8
        save(strcat('..\nuclear\output_var_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_var');
        save(strcat('..\nuclear\output_scc_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_scc');
        save(strcat('..\nuclear\output_nuclear_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_nuclear');
        save(strcat('..\nuclear\output_pop.dat'),'output_pop');
    elseif mcmode==9
        save(strcat('..\nuclearmonte\output_var_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_var');
        save(strcat('..\nuclearmonte\output_scc_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_scc');
        save(strcat('..\nuclearmonte\output_nuclear_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_nuclear');
        save(strcat('..\nuclear\output_pop.dat'),'output_pop');
        save(strcat('..\nuclearmonte\output_emi_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_emi');

        %     save(strcat('..\nuclearmonte2\output_emi.dat'),'output_emi');
    elseif mcmode==10
        output_var2=output_var(121:(invest_trans(4)-2025+121),:,:); output_var=output_var2;
        save(strcat('..\invest\output_var_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_var');
        save(strcat('..\invest\output_invest_EndSav',num2str(endsav),'_MC',num2str(mc),'.dat'),'output_invest');
        save(strcat('..\invest\output_pop.dat'),'output_pop');
    end

end

end
