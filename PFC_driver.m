% Path to mechanism files, get_PFC_cell.m, and get_PFC_1layer.m
model_dir='/home/jason/models/dynasim/mechanisms/JSS_PFC';
cd(model_dir);

% generic input and state equations
input_def={'input(V)=iAMPA(V); monitor input; onset=50; offset=inf;';
           'iAMPA(V)=-gAMPA.*sAMPA(k,:).*(V-EAMPA); EAMPA=0; gAMPA=0;';
           'sAMPA=getPoissonGating(baseline,dcAMPA,0,0,0,onset,offset,tauAMPA,T,Npop); dcAMPA=0; tauAMPA=2; baseline=0;';
          };
state_equations=['dV/dt=(@current+input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;' input_def{:}];

% PFC models
spec=get_PFC_cell('DS02PYjs',1);
spec=get_PFC_1layer('DS02PYjs',1,'DS02FSjs',1,[],0);
spec=get_PFC_1layer('DS02PYjs',2,'DS02FSjs',1,[],0);
spec=get_PFC_1layer('DS02PYjs',2,'DS02FSjs',2,'DS02RSNPjs',1);
spec=ApplyModifications(spec,{'Es','equations',state_equations});
vary={'Es','Iapp',0;'Es','gAMPA',[1 2]*1e-5;'Es','dcAMPA',20e3;'Es','offset',1500};
solver_options={'tspan',[0 200],'solver','rk1','dt',.01,'compile_flag',1,'verbose_flag',1};
data=SimulateModel(spec,'vary',vary,solver_options{:});
PlotData(data);

%% PFC network with persistence
% input_def={'input(V)=iAMPA(V)+iNMDA(V)+iGABA(V); monitor input,iAMPA,iNMDA,iGABA; onset=50; offset=inf;';
input_def={'input(V)=iAMPA(V)+iNMDA(V)+iGABA(V); monitor input,iAMPA,iNMDA,iGABA; onset=50; offset=inf;';
           'iAMPA(V)=-gAMPA.*sAMPA(k,:).*(V-EAMPA); EAMPA=0; gAMPA=0;';
           'iNMDA(V)=-gNMDA.*sNMDA(k,:).*(1.50265./(1+0.33*exp(V./(-16)))).*(V-ENMDA); ENMDA=0; gNMDA=0;';
           'iGABA(V)=-gGABA.*sGABA(k,:).*(V-EGABA); EGABA=-75; gGABA=0;';
           'sAMPA=getPoissonGating(bgAMPA,dcAMPA,0,0,0,onset,offset,tauAMPA,T,Npop); dcAMPA=0; bgAMPA=0; tauAMPA=2;';
           'sNMDA=getPoissonGating(bgNMDA,dcNMDA,0,0,0,onset,offset,tauNMDA,T,Npop); dcNMDA=0; bgNMDA=0; tauNMDA=150;';
           'sGABA=getPoissonGating(bgGABA,dcGABA,0,0,0,onset,offset,tauGABA,T,Npop); dcGABA=0; bgGABA=0; tauGABA=5';
           };
    EL='Ed'; IL='FS'; Ne=20; Ni=5;
    % input parameters for pyramidal cells
    gAMPAe=1e-3;        % uS, AMPA->PY
    gNMDAe=gAMPAe/50;   % uS, NMDA->PY
    gGABAe=.6e-3;       % uS, GABA->PY
    bgAMPAe=20000;      % Hz, (4000 input neurons firing rate <FR>=5Hz)
    bgNMDAe=bgAMPAe;    % Hz
    bgGABAe=13125;      % Hz
    dcAMPAe=0;      % Hz
    dcNMDAe=0;    % Hz
    dcGABAe=0;      % Hz
    gAMPAi=.74e-3;      % uS, AMPA->FS
    gNMDAi=gAMPAi/50;   % uS, NMDA->FS
    gGABAi=.6e-3;       % uS, GABA->FS
    bgAMPAi=bgAMPAe/2;  % Hz
    bgNMDAi=bgNMDAe/2;  % Hz
    bgGABAi=bgGABAe/2;  % Hz
    dcAMPAi=0;      % Hz
    dcNMDAi=0;    % Hz
    dcGABAi=0;      % Hz
    gAMPAee=3e-3*(100/Ne);   % 0.0150
    gNMDAee=gAMPAee/50;      % 0.0003
    gGABAie=.2e-3*(37/Ni);   % 0.0015
    gAMPAei=.74e-3*(100/Ne); % 0.0037
    gNMDAei=gAMPAei/50;      % .000074
    gGABAii=.6e-3*(37/Ni);   % 0.0044
    tauNMDAr=2.3; tauNMDA=95;
tauNMDAinp=150; 
gAMPAe=1e-3;   gNMDAe=gAMPAe/50; gGABAe=.6e-3; bgAMPAe=20000;     bgNMDAe=bgAMPAe;   bgGABAe=13125;     dcAMPAe=0; dcNMDAe=0; dcGABAe=0;
gAMPAi=.74e-3; gNMDAi=gAMPAi/50; gGABAi=.6e-3; bgAMPAi=bgAMPAe/2; bgNMDAi=bgNMDAe/2; bgGABAi=bgGABAe/2; dcAMPAi=0; dcNMDAi=0; dcGABAi=0;

gAMPAi=0; gNMDAi=0; gGABAi=0; gAMPAe=[.25 .5 .75]*1e-3; gNMDAe=gAMPAe(1)/50; gGABAe=.6e-3;
gAMPAi=0; gNMDAi=0; gGABAi=0; gAMPAe=.5e-3; gNMDAe=gAMPAe(1)/50; gGABAe=.6e-3;
gAMPAi=0; gNMDAi=0; gGABAi=0; gAMPAe=.65e-3; gNMDAe=gAMPAe(1)/50; gGABAe=.6e-3;
gAMPAi=.1e-3; gNMDAi=gAMPAi/50; gGABAi=.1e-3; gAMPAe=.55e-3; gNMDAe=gAMPAe(1)/50; gGABAe=.6e-3;
dcAMPAe=20e3; dcNMDAe=dcAMPAe;

spec=get_PFC_1layer('DS02PYjs',Ne,'DS02FSjs',Ni,[],0);
state_equations=['dV/dt=(@current+input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;' input_def{:}];
spec=ApplyModifications(spec,{EL,'equations',state_equations});
E_input_parameters={'gAMPA',gAMPAe,'gNMDA',gNMDAe,'gGABA',gGABAe,'dcAMPA',dcAMPAe,'dcNMDA',dcNMDAe,'dcGABA',dcGABAe,'bgAMPA',bgAMPAe,'bgNMDA',bgNMDAe,'bgGABA',bgGABAe};
spec.populations(2).parameters=cat(2,spec.populations(2).parameters,E_input_parameters);
spec=ApplyModifications(spec,{IL,'equations',state_equations});
I_input_parameters={'gAMPA',gAMPAi,'gNMDA',gNMDAi,'gGABA',gGABAi,'dcAMPA',dcAMPAi,'dcNMDA',dcNMDAi,'dcGABA',dcGABAi,'bgAMPA',bgAMPAi,'bgNMDA',bgNMDAi,'bgGABA',bgGABAi};
spec.populations(3).parameters=cat(2,spec.populations(3).parameters,I_input_parameters);
Kee=ones(Ne)-eye(Ne); 
spec=ApplyModifications(spec,{'Es->Ed','netcon',Kee});
Kii=ones(Ni)-eye(Ni); 
spec=ApplyModifications(spec,{'FS->FS','netcon',Kii});

gAMPAee=.01; gNMDAee=1; onset=4000; offset=5000; tauNMDA=150; 
gAMPAee=.002; gNMDAee=.0001; gGABAie=.002; onset=4000; offset=5000; tauNMDA=95;
gAMPAee=.002/4; gNMDAee=.0001*10; gGABAie=.002; gAMPAei=.00074; tauNMDA=150; tauNMDAr=10.6; gGABAii=.006; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
gAMPAee=.002/4; gNMDAee=.0001*20; gGABAie=.002; gAMPAei=gAMPAee/2; tauNMDA=95; tauNMDAr=10.6; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
  % need slightly larger E->I drive (maybe gAMPAei=.0005)
  % need larger E->E drive to achieve persistence
% NO INHIBITION:
gAMPAee=.002/4; gNMDAee=.0001*50; gGABAie=.00; gAMPAei=gAMPAee/2; tauNMDA=95; tauNMDAr=10.6; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=6000;
gAMPAee=.002/4; gNMDAee=.0001*50; gGABAie=.00; gAMPAei=gAMPAee/2; tauNMDA=95; tauNMDAr=20; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=6000;
% ** HAS PERSISTENCE WITHOUT INHIBITION -----------------
gAMPAee=.003/4; gNMDAee=.0001*150; gGABAie=.00; gAMPAei=gAMPAee/2; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=15e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=6000;
gAMPAee=.003/4; gNMDAee=.0001*150; gGABAie=.00; gAMPAei=gAMPAee/2; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=15.5e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
% **
% WITH INHIBITION:
gAMPAee=.003/4; gNMDAee=.0001*150; gGABAie=.2e-3; gAMPAei=gAMPAee/2; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=15.5e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
gAMPAee=.003/4; gNMDAee=.0001*150; gGABAie=.4e-3; gAMPAei=gAMPAee/2; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=15.5e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
gAMPAee=.00075; gNMDAee=.015; gGABAie=.0010; gAMPAei=.000555; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=16.5e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
dcAMPAe=20e3; dcNMDAe=dcAMPAe; gAMPAee=.00075; gNMDAee=.0175; gGABAie=.0010; gAMPAei=.001; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=16e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.00075; gNMDAee=.0175; gGABAie=.0010; gAMPAei=.001; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=16e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.0009; gNMDAee=.0175; gGABAie=.0010; gAMPAei=.001; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=16e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.001; gNMDAee=.0175; gGABAie=.0025; gAMPAei=.001; bgAMPAe=16.3e3; tauNMDA=95; tauNMDAr=10.6; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.001; gNMDAee=.0175; gGABAie=.0015; gAMPAei=.0037; bgAMPAe=16.3e3; tauNMDA=95; tauNMDAr=10.6; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
% gAMPAi=.2e-3; gNMDAi=gAMPAi/50; gGABAi=gAMPAi*(.6/.74); dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.001; gNMDAee=.0175; gGABAie=.0015; gAMPAei=.0037; bgAMPAe=16.3e3; tauNMDA=95; tauNMDAr=10.6; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;

% % ------------------------- TESTING
% % orig: gAMPAe=1e-3; gAMPAi=.74e-3; bgAMPAe=20e3; bgGABAe=13125; gAMPAee=0.015; gNMDAee=.0003; gGABAie=.0015; gAMPAei=.0037;
% 
% gAMPAe=.55e-3;                               gAMPAi=.74*gAMPAe;             gGABAi=gAMPAi*(.6/.74); gNMDAi=gAMPAi/50; gNMDAe=gAMPAe/50;
% bgAMPAe=17.1e3; bgGABAe=13125;                              dcAMPAe=15e3;                  dcNMDAe=dcAMPAe; bgNMDAe=bgAMPAe; bgAMPAi=bgAMPAe/2; bgNMDAi=bgNMDAe/2; bgGABAi=bgGABAe/2;
% gAMPAee=.001; gNMDAee=.0175;                  gGABAie=.0015; gAMPAei=.0037; tauNMDA=95; tauNMDAr=10.6; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
% tauNMDAinp=150; tauNMDA=150;
% % gAMPAe=1e-3; gAMPAi=.74e-3;                  gGABAi=gAMPAi*(.6/.74); gNMDAi=gAMPAi/50; gNMDAe=gAMPAe/50;
% % bgAMPAe=20e3; bgGABAe=13125; dcAMPAe=15e3;  dcNMDAe=dcAMPAe; bgNMDAe=bgAMPAe; bgAMPAi=bgAMPAe/2; bgNMDAi=bgNMDAe/2; bgGABAi=bgGABAe/2;
% % gAMPAee=.015; gNMDAee=.0003; gGABAie=.0015; gAMPAei=.0037; tauNMDA=95; tauNMDAr=10.6; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
% dcAMPAe=0;dcNMDAe=0;
% -------------------------

vary={EL,'dcAMPA',dcAMPAe;EL,'dcNMDA',dcNMDAe; IL,'gAMPA',gAMPAi;IL,'gNMDA',gNMDAi;IL,'gGABA',gGABAi;IL,'bgAMPA',bgAMPAi;IL,'bgNMDA',bgNMDAi;IL,'bgGABA',bgGABAi;EL,'gAMPA',gAMPAe;EL,'gNMDA',gNMDAe;EL,'gGABA',gGABAe;EL,'bgAMPA',bgAMPAe;EL,'bgNMDA',bgNMDAe;EL,'bgGABA',bgGABAe;'FS->Es','gGABA',gGABAie;'Es->Ed','gNMDA',gNMDAee;'(Es->Ed,Es->FS)','tauNMDA',tauNMDA;'(Es->Ed,Es->FS)','tauNMDAr',tauNMDAr;'Es->Ed','gAMPA',gAMPAee;'(Es,Ed,FS)','onset',onset;'(Es,Ed,FS)','offset',offset;'(Es,Ed,FS)','Iapp',0;'Es->FS','gAMPA',gAMPAei;'Es->FS','gNMDA',gNMDAei;'FS->FS','gGABA',gGABAii;'(Es,Ed,FS)','tauNMDA',tauNMDAinp};%;
% vary={EL,'dcAMPA',dcAMPAe;EL,'dcNMDA',dcNMDAe; IL,'gAMPA',gAMPAi;IL,'gNMDA',gNMDAi;IL,'gGABA',gGABAi;IL,'bgAMPA',bgAMPAi;IL,'bgNMDA',bgNMDAi;IL,'bgGABA',bgGABAi;EL,'gAMPA',gAMPAe;EL,'gNMDA',gNMDAe;EL,'gGABA',gGABAe;EL,'bgAMPA',bgAMPAe;EL,'bgNMDA',bgNMDAe;EL,'bgGABA',bgGABAe;'FS->Es','gGABA',gGABAie;'Es->Ed','gNMDA',gNMDAee;'(Es->Ed,Es->FS)','tauNMDA',tauNMDA;'(Es->Ed,Es->FS)','tauNMDAr',tauNMDAr;'Es->Ed','gAMPA',gAMPAee;'(Es,Ed,FS)','onset',onset;'(Es,Ed,FS)','offset',offset;'(Es,Ed,FS)','Iapp',0;'Es->FS','gAMPA',gAMPAei;'Es->FS','gNMDA',gNMDAei;'FS->FS','gGABA',gGABAii};%;'(Es,Ed,FS)','tauNMDA',tauNMDA
solver_options={'tspan',[0 10000],'solver','rk1','dt',.01,'compile_flag',1,'verbose_flag',1};
data=SimulateModel(spec,'vary',vary,solver_options{:});
PlotData(data,'plot_type','rastergram','threshold',-20);

% Need: higher spontaneous FR and lower delay FR in E-cells; a more prominant role for inhibition
% Q: does decreasing stim (dc) decrease delay FR and/or kill delay activity?
%   A: it decreased the stimulus-period FR but did not affect the delay-period FR; --> delay FR relates to recurrent synaptic properties
% Q: can very slightly increasing bgAMPA&bgNMDA causes spont FR .5-3Hz and still get state transition?
% Q: how to maintain effect w/ the addition of inhibition?

% PlotData(data,'variable',{'Ed_Es_iNMDA_INMDA','Es_V'});
% PlotData(data,'variable',{'Ed_V'});
data.varied={'Es_Ed_gNMDA'}; PlotFR(data,'threshold',-20);

t=data(1).time; dur=offset-onset;
tix1=(t>=onset-dur)&(t<onset);
tix2=(t>onset)&(t<=offset);
tix3=(t>offset)&(t<=offset+dur);
isyn=[[sum(sum(data.Ed_Es_iAMPA_IAMPA(tix1,:))) sum(sum(data.Ed_Es_iNMDA_INMDA(tix1,:))) sum(sum(data.Es_FS_iGABA_IGABA(tix1,:))) sum(sum(data.FS_Es_iAMPA_IAMPA(tix1,:))) sum(sum(data.FS_Es_iNMDA_INMDA(tix1,:))) sum(sum(data.FS_FS_iGABA_IGABA(tix1,:)))]' ...
[sum(sum(data.Ed_Es_iAMPA_IAMPA(tix2,:))) sum(sum(data.Ed_Es_iNMDA_INMDA(tix2,:))) sum(sum(data.Es_FS_iGABA_IGABA(tix2,:))) sum(sum(data.FS_Es_iAMPA_IAMPA(tix2,:))) sum(sum(data.FS_Es_iNMDA_INMDA(tix2,:))) sum(sum(data.FS_FS_iGABA_IGABA(tix2,:)))]' ...
[sum(sum(data.Ed_Es_iAMPA_IAMPA(tix3,:))) sum(sum(data.Ed_Es_iNMDA_INMDA(tix3,:))) sum(sum(data.Es_FS_iGABA_IGABA(tix3,:))) sum(sum(data.FS_Es_iAMPA_IAMPA(tix3,:))) sum(sum(data.FS_Es_iNMDA_INMDA(tix3,:))) sum(sum(data.FS_FS_iGABA_IGABA(tix3,:)))]']
iext=-[[sum(sum(data.Ed_iAMPA(tix1,:))) sum(sum(data.Ed_iNMDA(tix1,:))) sum(sum(data.Ed_iGABA(tix1,:))) sum(sum(data.FS_iAMPA(tix1,:))) sum(sum(data.FS_iNMDA(tix1,:))) sum(sum(data.FS_iGABA(tix1,:)))]' ...
 [sum(sum(data.Ed_iAMPA(tix2,:))) sum(sum(data.Ed_iNMDA(tix2,:))) sum(sum(data.Ed_iGABA(tix2,:))) sum(sum(data.FS_iAMPA(tix2,:))) sum(sum(data.FS_iNMDA(tix2,:))) sum(sum(data.FS_iGABA(tix2,:)))]' ...
 [sum(sum(data.Ed_iAMPA(tix3,:))) sum(sum(data.Ed_iNMDA(tix3,:))) sum(sum(data.Ed_iGABA(tix3,:))) sum(sum(data.FS_iAMPA(tix3,:))) sum(sum(data.FS_iNMDA(tix3,:))) sum(sum(data.FS_iGABA(tix3,:)))]']
isyn+iext

%PlotData(data,'variable',{'Ed_Es_iAMPA_IAMPA','Ed_Es_iNMDA_INMDA','Es_FS_iGABA_IGABA'},'xlim',[0 100],'ylim',[-1 1])
figure('position',[50 80 1800 850]); xlims=[0 data(1).time(end)];
subplot(5,1,1); plot(t,data(1).Ed_Es_iNMDA_INMDA(:,1)); ylabel('Ed.NMDA'); line(xlim,[0 0]); line(xlim,[-.2 -.2],'color','r'); axis tight
subplot(5,1,2); plot(t,data(1).Ed_Es_iAMPA_IAMPA(:,1)); ylabel('Ed.AMPA'); ylim([-1 0]); xlim(xlims); line(xlim,[0 0]); line(xlim,[-.2 -.2],'color','r');
subplot(5,1,3); plot(t,data(1).Es_FS_iGABA_IGABA(:,1)); ylabel('Es.GABA'); ylim([0 1]); xlim(xlims); line(xlim,[.2 .2],'color','r');
subplot(5,1,4); plot(t,data(1).Es_V(:,1)); ylabel('Es\_V'); xlim(xlims); ylim([-100 50]);
subplot(5,1,5); plot(t,data(1).Ed_input(:,1)); ylabel('input'); xlim(xlims); line(xlim,[0 0]); ylim([0 2]); line(xlim,[.2 .2],'color','r');

% sims w/ persistence:
% using fixed input tauNMDA=150 and variable synaptic tauNMDA: 
% vary={EL,'dcAMPA',dcAMPAe;EL,'dcNMDA',dcNMDAe; IL,'gAMPA',gAMPAi;IL,'gNMDA',gNMDAi;IL,'gGABA',gGABAi;IL,'bgAMPA',bgAMPAi;IL,'bgNMDA',bgNMDAi;IL,'bgGABA',bgGABAi;EL,'gAMPA',gAMPAe;EL,'gNMDA',gNMDAe;EL,'gGABA',gGABAe;EL,'bgAMPA',bgAMPAe;EL,'bgNMDA',bgNMDAe;EL,'bgGABA',bgGABAe;'FS->Es','gGABA',gGABAie;'Es->Ed','gNMDA',gNMDAee;'(Es->Ed,Es->FS)','tauNMDA',tauNMDA;'(Es->Ed,Es->FS)','tauNMDAr',tauNMDAr;'Es->Ed','gAMPA',gAMPAee;'(Es,Ed,FS)','onset',onset;'(Es,Ed,FS)','offset',offset;'(Es,Ed,FS)','Iapp',0;'Es->FS','gAMPA',gAMPAei;'Es->FS','gNMDA',gNMDAei;'FS->FS','gGABA',gGABAii};%;
  % w/o inhib: gAMPAee=.003/4; gNMDAee=.0001*150; gGABAie=.00;   gAMPAei=gAMPAee/2;  bgAMPAe=15e3; bgNMDAe=bgAMPAe; gGABAii=.0044; tauNMDA=95; tauNMDAr=10.6; gNMDAei=gAMPAei/50; onset=4000; offset=6000;
  % w/ inhib:  gAMPAee=.003/4; gNMDAee=.0001*150; gGABAie=.0010; gAMPAei=.00037*1.5; bgAMPAe=16e3; bgNMDAe=bgAMPAe; gGABAii=.0044; tauNMDA=95; tauNMDAr=10.6; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
  % w/ inhib:  gAMPAee=.00075; gNMDAee=.0175;     gGABAie=.0010; gAMPAei=.001;       bgAMPAe=16e3; dcAMPAe=20e3; dcNMDAe=dcAMPAe;bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;tauNMDA=95; tauNMDAr=10.6; 
  % dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.00075; gNMDAee=.0175; gGABAie=.0010; gAMPAei=.001; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=16e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
  % dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.0009; gNMDAee=.0175; gGABAie=.0010; gAMPAei=.001; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=16e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
  % dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.001; gNMDAee=.0175; gGABAie=.002; gAMPAei=.001; tauNMDA=95; tauNMDAr=10.6; bgAMPAe=16e3; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
  % dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.001; gNMDAee=.0175; gGABAie=.0015; gAMPAei=.0037; bgAMPAe=16.3e3; tauNMDA=95; tauNMDAr=10.6; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;
  % gAMPAi=.2e-3; gNMDAi=gAMPAi/50; gGABAi=gAMPAi*(.6/.74); dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.001; gNMDAee=.0175; gGABAie=.0015; gAMPAei=.0037; bgAMPAe=16.3e3; tauNMDA=95; tauNMDAr=10.6; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;

% set tauNMDAr=13.89 and tauNMDA=150 in Es->Ed and Es->FS (Methods in Neuronal Modeling)
  
  
%{
1. tune gAMPAe (and maybe baseline) s.t. <FR>=.5-3Hz (as in [DS02])
    gAMPAe=[.2 .4 .6 .8 1]*1e-4; dcAMPA=0; baseline=20e3; gAMPAee=0; gNMDAee=0; onset=500; offset=1000;
    gAMPAe=[.15 .2 .25 .3 .35]*1e-4; dcAMPA=0; baseline=20e3; gAMPAee=0; gNMDAee=0; onset=500; offset=1000;

2. tune dcAMPA s.t. <FR> increases noticeably
3. tune gNMDAee s.t. activity persists after offset
%}


%% PFC network with getBurstGating (inputs varying sync and freq; nets varying resonance and integration time (taum, tauSYN), and lateral inhibition (gIE,gEI))



