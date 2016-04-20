% PFC_competition.m
% model assumption: single-layer E/I network

% Path to mechanism files, get_PFC_1layer.m, and get_PFC_cell.m
model_dir='/home/jason/models/dynasim/mechanisms/JSS_PFC';
cd(model_dir);

TargetPop='Ed'; % name of population(s) to stimulate (for inputs and recurrent connections Es->Ed)
OutputPop='Es'; % name of population to analyze (for post-simulation analysis)
EiPop='Es';     % excitatory population input compartment for E/I connectivity
EoPop='Es';     % excitatory population output compartment for E/I connectivity
IPop='FS';      % inhibitory population input=output compartment for E/I connectivity

% set population sizes
assembly_size=8;                  % # E-cells per assembly {4,8,12,16,20}
num_assemblies=1;                 % # assemblies
Ne=assembly_size*num_assemblies;  % total # E-cells
Ni=.25*Ne;                        % total # I-cells

% load baseline model
base=get_PFC_1layer('DS02PYjs',Ne,'DS02FSjs',Ni,[],0);
TargetPopIndex=find(strcmp(TargetPop,{base.populations.name}));
OutputPopIndex=find(strcmp(OutputPop,{base.populations.name}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TUNING PROFILES
% Description: one assembly receiving poisson and burst inputs. assess
% spiking for varying input frequency and synchrony. compare for networks
% of cells with fast or slow input integration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec=base;

% generic inputs and dependent state equations
input_def={'input(V)=iPoisson(V)+iBurst(V); monitor input,iBurst,iPoisson; onset=50; offset=inf;';
           'iPoisson(V)=gAMPA.*s1a(k,:).*(V-EAMPA); EAMPA=0; gAMPA=0;';
           'iBurst(V)=gAMPA.*s1b(k,:).*(V-EAMPA);';
           's1a=getPoissonGating(baseline,dcAMPA1,acAMPA1,fAMPA1,phiAMPA1,onset,offset,tauAMPA,T,Npop,ones(1,Npop),kick,ramp_dc_flag1,ramp_ac_flag1);';
           's1b=getBurstGating(T,fBURST1,widthBURST1,nspksBURST1,dcBURST1,Npop,minIBI,meanIBI,tauBURST,kick,ones(1,Npop),shared_sources_flag,onset,offset,ramp_dc_flag1,ramp_ac_flag1,num_sources);';
           'baseline=0; dcAMPA1=0; acAMPA1=0; fAMPA1=0; phiAMPA1=0; tauAMPA=2; kick=1; ramp_dc_flag1=0; ramp_ac_flag1=0;';
           'fBURST1=0; widthBURST1=10; nspksBURST1=100; dcBURST1=0; minIBI=10; meanIBI=100; tauBURST=2; shared_sources_flag=0; num_sources=1;';
           };
state_equations=['dV/dt=(@current-input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;' input_def{:}];
input_parameters={}; % input and biophysical parameters (key/value pairs)

% connectivity (all one assembly for tuning)
% kernels
Kee=ones(Ne,Ne); % E->E, [N_pre x N_post]
Kei=ones(Ne,Ni); % E->I, [N_pre x N_post]
Kie=ones(Ni,Ne); % I->E, [N_pre x N_post]
Kii=ones(Ni,Ni); % I->I, [N_pre x N_post]
% max conductance
gie=.25;      % .5,  .25  E->I
gii=0;        % 0,   0    I->I
geiAMPA=.75;  % .25, .75  E->I
geiNMDA=.75;  % 1,   .75  E->I
geeAMPA=.1;   % 0,   .1   E->E
geeNMDA=.1;   % 1,   .1   E->E
% time constants
tauAMPA=2;
tauNMDA=95;
tauGABA=5;
% normalize kernels by number of presynaptic connections
Kee=Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1]);
Kei=Kei./repmat(max(1,sum(Kei,1)),[size(Kei,1) 1]);
Kie=Kie./repmat(max(1,sum(Kie,1)),[size(Kie,1) 1]);
Kii=Kii./repmat(max(1,sum(Kii,1)),[size(Kii,1) 1]);

% update model specification
% populations (add inputs and input parameters)
spec=ApplyModifications(spec,{TargetPop,'equations',state_equations}); % Ed
spec.populations(TargetPopIndex).parameters=cat(2,input_parameters,RemoveKeyval(spec.populations(TargetPopIndex).parameters,input_parameters(1:2:end)));
% connections (adjust connectivity kernels/parameters)
L=[EoPop '->' TargetPop]; % Es->Ed
spec=ApplyModifications(spec,{L,'netcon',Kee;L,'gAMPA',geeAMPA;L,'gNMDA',geeNMDA;L,'tauAMPA',tauAMPA;L,'tauNMDA',tauNMDA});
L=[EoPop '->' IPop];      % Es->FS
spec=ApplyModifications(spec,{L,'netcon',Kei;L,'gAMPA',geiAMPA;L,'gNMDA',geiNMDA;L,'tauAMPA',tauAMPA;L,'tauNMDA',tauNMDA});
L=[IPop '->' EiPop];      % FS->Es
spec=ApplyModifications(spec,{L,'netcon',Kie;L,'gGABA',gie;L,'tauGABA',tauGABA});
L=[IPop '->' IPop];       % FS->FS
spec=ApplyModifications(spec,{L,'netcon',Kii;L,'gGABA',gii;L,'tauGABA',tauGABA});

% -------------------------------------------------------------------------
% simulation controls
tspan=[0 2000];  % [beg end], ms
dt=.01;         % fixed time step, ms
solver='rk1';   % numerical integration method {'rk1','rk2','rk4'}
compile_flag=1; % whether to compile simulation
simulator_options={'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1};
mods=[];
vary=[];
% -------------------------------------------------------------------------
% data=SimulateModel(ApplyModifications(spec,mods),'vary',vary,simulator_options{:});  
% PlotData(data);


f=0; w=[10 10 20 20 30 30];
f=0; w=[10:40:130]; minIBI=max(w); meanIBI=2*minIBI; nspks=100;
f=0:10:60; w=10; nspks=100;
f=10; w=[5 10 20]; nspks=40;
vary={'Ed','gAMPA',1e-3;'Ed','fBURST1',f;'Ed','meanIBI',meanIBI;'Ed','minIBI',minIBI;'Ed','widthBURST1',w;'Ed','nspksBURST1',nspks;'Es->Ed','(gAMPA,gNMDA)',0;'FS->Es','gGABA',0};
data=SimulateModel(spec,'vary',vary,simulator_options{:});
PlotData(data,'plot_type','rastergram');
PlotData(data,'variable','Es_V');
PlotData(data,'variable','Ed_input');

data=CalcFR(data,'bin_size',30,'bin_shift',10);
gating=arrayfun(@(x)mean(x.model.fixed_variables.Ed_s1b(:)),data);
input=cellfun(@(x)mean(-x(:)),{data.Ed_input});
output=cellfun(@(x)mean(x(:)),{data.Es_V_FR});
figure('position',[190 420 1430 420]); 
plot(w,output./input); xlabel('1/sync'); ylabel('output/input');

figure('position',[190 420 1430 420]); 
subplot(1,4,1); plot(input,output,'o--'); xlabel('input'); ylabel('output');
subplot(1,4,2); plot(gating,input,'o--'); xlabel('gating'); ylabel('input');
subplot(1,4,3); plot(w,input,'o--'); xlabel('1/sync'); ylabel('input');
subplot(1,4,4); plot(w,output,'o--'); xlabel('1/sync'); ylabel('output');

PlotData(data,'variable','Ed_input');
PlotData(data,'variable','Es_V');

PlotData(data,'plot_type','rastergram');

% need to adjust:
% E->I is too strong (I-cells saturate with any E spikes)
% other connections are probably too strong as well
% inputs are too strong (all cells spike together)

% todo: check for relationship b/w (# generated spikes) and width using
% getBurstGating() directly.


% Experiments and analysis
% Experiment: tuning profile: (frequency & spike synchrony) -> (resonance & coincidence detection)
%  1. FR=f(sync), IBI~Exp (freq=0), fast (50%) and slow (200%) taum
%  2. FR=f(freq), IBI=1/freq, low and high sync
%  3. FR=f(sync,freq)

% FR=f(sync), IBI~Exp (freq=0), fast (50%) and slow (200%) taum
vary={'Ed','fBURST1',0;'Ed','widthBURST1',[.01:.01:.1];'(Es,Ed)','Cm',[1 1.5];'FS->Ed','gGABA',0};
data=SimulateModel(spec,'vary',vary);
PlotFR(data);

% FR=f(freq), IBI=1/freq, low and high sync
vary={'Ed','fBURST1',[0:10:60];'Ed','widthBURST1',[.01 .1];'FS->Ed','gGABA',0};
data=SimulateModel(spec,'vary',vary);
PlotFR(data);

% input/ouput
data=CalcFR(data); 
[mean(data(1).Ed_input(:)) mean(cellfun(@length,data(1).Es_V_spike_times))]

% add inputs (getBurstGating)
% INPUTS: poisson vs sync(exp_burstISI) vs sync(rhythmic_burstISI)
%   controlled sync: use getBurstGating w/ width and num_spikes
%   rhythmic controlled sync: use getBurstGating w/ width, num_spikes, and f>0
%   poisson: use getPoissonGating w/ baseline and DC
%   rhythmic: use getPoissonGating w/ freq and AC
input_def={'input(V)=I1(V)+I2(V); monitor input,I1,I2; onset=50; offset=inf;';
           'I1(V)=gAMPA.*(s1a(k,:)+s1b(k,:)).*(V-EAMPA); EAMPA=0; gAMPA=0;';
           'I2(V)=gAMPA.*(s2a(k,:)+s2b(k,:)).*(V-EAMPA); EAMPA=0; gAMPA=0;';
           's1a=getPoissonGating(bgAMPA,dcAMPA,acAMPA,freqAMPA,phiAMPA,onset,offset,tauAMPA,T,Npop); dcAMPA=0; acAMPA=0; freqAMPA=0; bgAMPA=0; phiAMPA=0; tauAMPA=2;';
           's1b=getBurstGating(bgAMPA,dcAMPA,acAMPA,freqAMPA,phiAMPA,onset,offset,tauAMPA,T,Npop); dcAMPA=0; acAMPA=0; freqAMPA=0; bgAMPA=0; phiAMPA=0; tauAMPA=2;';
           's2a=getPoissonGating(bgAMPA,dcAMPA,acAMPA,freqAMPA,phiAMPA,onset,offset,tauAMPA,T,Npop); dcAMPA=0; acAMPA=0; freqAMPA=0; bgAMPA=0; phiAMPA=0; tauAMPA=2;';
           's2b=getBurstGating(bgAMPA,dcAMPA,acAMPA,freqAMPA,phiAMPA,onset,offset,tauAMPA,T,Npop); dcAMPA=0; acAMPA=0; freqAMPA=0; bgAMPA=0; phiAMPA=0; tauAMPA=2;';
           };
state_equations=['dV/dt=(@current-input(V)+Iapp*(t>onset&t<offset))./Cm; Cm=1; Iapp=0; V(0)=-65;' input_def{:}];
s=ApplyModifications(spec,{'(Es,Ed)','equations',state_equations})

se1a=sprintf('s1a=getPoissonGating(baseline,dcAMPA1,acAMPA1,fAMPA1,phiAMPA1,onset,offset,tauAMPA,T,Npop,%s,kick,ramp_dc_flag1,ramp_ac_flag1)',toString(Ke1));
se1b=sprintf('s1b=getBurstGating(T,fBURST1,widthBURST1,nspksBURST1,dcBURST1,Npop,minIBI,meanIBI,tauBURST,kick,%s,shared_sources_flag,onset,offset,ramp_dc_flag1,ramp_ac_flag1,num_sources)',toString(Ke1));
se2a=sprintf('s2a=getPoissonGating(baseline,dcAMPA2,acAMPA2,fAMPA2,phiAMPA2,onset,offset,tauAMPA,T,Npop,%s,kick,ramp_dc_flag2,ramp_ac_flag2)',toString(Ke2));
se2b=sprintf('s2b=getBurstGating(T,fBURST2,widthBURST2,nspksBURST2,dcBURST2,Npop,minIBI,meanIBI,tauBURST,kick,%s,shared_sources_flag,onset,offset,ramp_dc_flag2,ramp_ac_flag2,num_sources)',toString(Ke2));
si1a=sprintf('s1a=getPoissonGating(baseline,dcAMPA1,acAMPA1,fAMPA1,phiAMPA1,onset,offset,tauAMPA,T,Npop,%s,kick,ramp_dc_flag1,ramp_ac_flag1)',toString(ones(Ni,1)));
si1b=sprintf('s1b=getBurstGating(T,fBURST1,widthBURST1,nspksBURST1,dcBURST1,Npop,minIBI,meanIBI,tauBURST,kick,%s,shared_sources_flag,onset,offset,ramp_dc_flag1,ramp_ac_flag1,num_sources)',toString(ones(Ni,1)));
E_input=sprintf('input(t)=-gAMPA.*(s1a(k,:)+s1b(k,:)+s2a(k,:)+s2b(k,:)).*(X-EAMPA); monitor input');
I_input=sprintf('input(t)=-gAMPA.*(s1a(k,:)+s1b(k,:)).*(X-EAMPA); monitor input');

% default input parameters (background + nonrhythmic sync bursts)
% E-cell inputs:
dcAMPA1=0; acAMPA1=0; fAMPA1=0; phiAMPA1=0; fBURST1=0; widthBURST1=10; nspksBURST1=100; dcBURST1=0; ramp_dc_flag1=0; ramp_ac_flag1=0; 
dcAMPA2=0; acAMPA2=0; fAMPA2=0; phiAMPA2=0; fBURST2=0; widthBURST2=10; nspksBURST2=100; dcBURST2=0; ramp_dc_flag2=0; ramp_ac_flag2=0;
gAMPA=1e-4; baseline=1; onset=50; offset=500; tauAMPA=2; tauBURST=2; minIBI=widthBURST1; meanIBI=10*minIBI; kick=1; EAMPA=0; shared_sources_flag=0;

'dcBURST1',dcBURST1,'fBURST1',fBURST1,'widthBURST1',widthBURST1,'nspksBURST1',nspksBURST1,'dcAMPA1',dcAMPA1,'acAMPA1',acAMPA1,'fAMPA1',fAMPA1,'phiAMPA1',phiAMPA1,'ramp_dc_flag1',ramp_dc_flag1,'ramp_ac_flag1',ramp_ac_flag1
'dcBURST2',dcBURST2,'fBURST2',fBURST2,'widthBURST2',widthBURST2,'nspksBURST2',nspksBURST2,'dcAMPA2',dcAMPA2,'acAMPA2',acAMPA2,'fAMPA2',fAMPA2,'phiAMPA2',phiAMPA2,'ramp_dc_flag2',ramp_dc_flag2,'ramp_ac_flag2',ramp_ac_flag2
'baseline',baseline,'onset',onset,'offset',offset,'tauAMPA',tauAMPA,'tauBURST',tauBURST,'minIBI',minIBI,'meanIBI',meanIBI,'kick',kick,'gAMPA',gAMPA,'EAMPA',EAMPA,'shared_sources_flag',shared_sources_flag

% null values:
'dcBURST1',0,'fBURST1',0,'widthBURST1',.1,'nspksBURST1',1,'dcAMPA1',0,'acAMPA1',0,'fAMPA1',0,'phiAMPA1',0,'ramp_dc_flag1',0,'ramp_ac_flag1',0
'dcBURST2',0,'fBURST2',0,'widthBURST2',.1,'nspksBURST2',1,'dcAMPA2',0,'acAMPA2',0,'fAMPA2',0,'phiAMPA2',0,'ramp_dc_flag2',0,'ramp_ac_flag2',0
'baseline',0,'onset',0,'offset',inf,'tauAMPA',2,'tauBURST',2,'minIBI',.1,'meanIBI',.1,'kick',1,'gAMPA',gAMPA,'EAMPA',EAMPA,'shared_sources_flag',shared_sources_flag

T=0:.01:1000; f=10; w=10; nspk=10; dc=0; Npop=2; minIBI=0; meanIBI=0; [sext,ts,Ptot,Peg]=getBurstGating(T,f,w,nspk,dc,Npop,minIBI,meanIBI,2,1,ones(1,Npop),0,0,inf,0,0,1);
PlotData(sext);

onset=50; offset=inf; Npop=Ne; T=0:.01:1000;
baseline=0; dcAMPA1=0; acAMPA1=0; fAMPA1=0; phiAMPA1=0; tauAMPA=2; kick=1; ramp_dc_flag1=0; ramp_ac_flag1=0;
fBURST1=0; widthBURST1=10; nspksBURST1=100; dcBURST1=0; minIBI=25; meanIBI=100; tauBURST=2; shared_sources_flag=0; num_sources=1;
s1b=getBurstGating(T,fBURST1,widthBURST1,nspksBURST1,dcBURST1,Npop,minIBI,meanIBI,tauBURST,kick,ones(1,Npop),shared_sources_flag,onset,offset,ramp_dc_flag1,ramp_ac_flag1,num_sources);
PlotData(s1b);

% I-cell inputs:
dcAMPAi=0; acAMPAi=0; fAMPAi=0; phiAMPAi=0; fBURSTi=0; widthi=10; num_spikesi=100; dcBURSTi=0; ramp_dc_flagi=0; ramp_ac_flagi=0; 
gAMPAi=1e-4; baselinei=1; onseti=50; offseti=500; tauAMPAi=2; tauBURSTi=2; minIBIi=widthBURST1; meanIBIi=10*minIBI; shared_sources_flagi=0;

spec.populations(TargetPopIndex).parameters=cat(2,E_input_parameters,RemoveKeyval(spec.populations(TargetPopIndex).parameters,E_input_parameters(1:2:end)));
% spec.populations(TargetPopIndex).parameters=cat(2,spec.populations(TargetPopIndex).parameters,E_input_parameters);
% p=cat(2,p,RemoveKeyval(s.populations(end).parameters,p(1:2:end)));

T=0:.01:1000; f=10; w=10; nspk=10; dc=0; Npop=2; minIBI=0; meanIBI=0; 
[sext,ts,Ptot,Peg]=getBurstGating(T,f,w,nspk,dc,1,Npop,minIBI,meanIBI,2,1,ones(Npop,1),0,0,inf,0,0);
PlotData(sext)
num_input_spikes=sum(Ptot,1);
num_output_spikes=cellfun(@length,data(1).spike_times_Es_V);
[num_input_spikes(roi1) num_input_spikes(roi2); num_output_spikes(roi1) num_output_spikes(roi2)]

% define assembly connectivity kernels
% 1. K1ee: one PY assembly (assembly_size=Ne/2), for tuning curves (set size=Ne/2 before sim)
% 2. K2ee: two PY assemblies (assembly_size=Ne/2), for competition

% Kee (E->E): default assemblies
Kee=zeros(Ne,Ne); % E->E, [N_pre x N_post]
  % connect E-cells within block
  block=ones(assembly_size)-eye(assembly_size);
  for i=1:num_assemblies
    ind=(i-1)*assembly_size+(1:assembly_size);
    Kee(ind,ind)=block;
  end
  if num_assemblies<=2
    % connect E-cells between blocks
    K12=(1-Kee-eye(Ne)); % connections b/w blocks 1 and 2
    Kee=Kee+bwblockee*K12;
  elseif bwblockee>0
    error('between-assembly recurrent excitation not supported for E-populations with more than 2 assemblies.');
  end
% default interneuron-related kernels:
Kei=zeros(Ne,Ni); % E->I, [N_pre x N_post] (none)
Kie=zeros(Ni,Ne); % I->E, [N_pre x N_post] (none)
Kii=ones(Ni,Ni); % I->I,
s=ApplyModifications(spec,{'Es->FS','netcon',Kei})
% ...



% experiments and analysis
% Experiment: tuning profile: (frequency & spike synchrony) -> (resonance & coincidence detection)
% Experiment: assembly competition (lateral inhibition)

% 1. compute tuning curves (1 assembly): 
%   - FR=f(sync), IBI~Exp (freq=0), fast (50%) and slow (200%) taum
%   - FR=f(freq), IBI=1/freq, low and high sync
%   - FR=f(sync,freq)
% 2. assembly competition (factors shaping relative activity b/w 2 assemblies)
%   - collect FRbias=f(sync1,freq1,sync2,freq2)
%   - compute FRi=dFR-FRbias
%   - split (sync,freq,inhib) contributions to dFR
%   - vary gEI, gIE, I0, E0

% FR=f(sync), IBI~Exp (freq=0), fast (50%) and slow (200%) taum
vary={'Ed','freq',0;'Ed','width',[.01:.01:.1];'(Es,Ed)','Cm',[1 1.5];'FS->Ed','gGABA',0};
data=SimulateModel(spec,'vary',vary);
PlotFR(data);

% FR=f(freq), IBI=1/freq, low and high sync
vary={'Ed','freq',[0:10:60];'Ed','width',[.01 .1];'FS->Ed','gGABA',0};
data=SimulateModel(spec,'vary',vary);
PlotFR(data);

% assembly competition
vary={'Ed','(f1,f2)',0;'Ed','[widthBURST1,widthBURST2]',.01;'FS->Ed','gGABA',[0:.1:1]};
data=SimulateModel(spec,'vary',vary);
PlotFR(data);

% run simulations
mods=[]; vary=[];
solver_options={'tspan',tspan,'solver','rk2','dt',dt,'compile_flag',compile_flag,'verbose_flag',1};
data=SimulateModel(ApplyModifications(s,mods),'vary',vary,solver_options{:});


% calc tuning curves
stats=CalcResonanceStats(data,'sweep_parameter','E_f1','repetition_parameter','E_rep');
fMUAmu=stats.E_v_Power_MUA.repetition_sets.PeakFreq_mu;
fMUAsd=stats.E_v_Power_MUA.repetition_sets.PeakFreq_sd;
FRmu=zeros(size(f)); FRse=zeros(size(f));
xc1mu=zeros(size(f)); xc1se=zeros(size(f));
for i=1:length(f)
  sel=([data.E_f1]==f(i));
  FRmu(i)=mean(stats.E_v_FR.pop_FR_mean(sel));
  FRse(i)=std(stats.E_v_FR.pop_FR_mean(sel));%/sqrt(nrep);
  xc1mu(i)=mean(xc1(sel));
  xc1se(i)=std(xc1(sel));%/sqrt(nrep)
end
FR=smooth(FRmu,2);
fc=f(find(FR==max(FR),1,'first'));

% spike sync b/w two ROIs
for i=1:length(data)
  stats(i)=CalcSpikeSync(data(i),'ROI_pairs',{'E_v',[0 .5],'E_v',[.5 1]});
end
xc1=arrayfun(@(x)x.pairs.xcmax_pops,stats);


% separate the populations being manipulated (i.e., those to which the
% input is added) and those that are being analyzed (i.e., those whose
% rastergrams and LFP/MUA are being compared). will be useful for extending
% beyond the single-layer PFC network to multi-layer networks and models of
% other regions (e.g., Ctx/TRN/Thalamus or Ctx/BG/Thalamus loops).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ASSEMBLY COMPETITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gAMPAee=3e-3*(100/Ne);   % 0.0150
    gNMDAee=gAMPAee/50;      % 0.0003
    gGABAie=.2e-3*(37/Ni);   % 0.0015
    gAMPAei=.74e-3*(100/Ne); % 0.0037
    gNMDAei=gAMPAei/50;      % .000074
    gGABAii=.6e-3*(37/Ni);   % 0.0044
    tauNMDAr=2.3; tauNMDA=95;
tauNMDAinp=150; 

dcAMPAe=15e3; dcNMDAe=dcAMPAe; gAMPAee=.001; gNMDAee=.0175; gGABAie=.0015; gAMPAei=.0037; bgAMPAe=16.3e3; tauNMDA=95; tauNMDAr=10.6; bgNMDAe=bgAMPAe; gGABAii=.0044; gNMDAei=gAMPAei/50; onset=4000; offset=5000;



% No competition (gIE=0):
% FR(w)
% FR(f)
% FR(w,f)
% repeat: taum=50%, taum=200%

% With competition (gIE>0):
% FRIbias(w1,w2,f1,f2)
  % f=0, IBI~exp; f!=0: f
  % w=MIN (max sync)

% where:
% w = width of burst ~ (1/sync)
% f = burst frequency = (1/IBI) where IBI = "inter-burst interval"

% input_def={'input(V)=iAMPA(V)+iNMDA(V)+iGABA(V); monitor input,iAMPA,iNMDA,iGABA; onset=50; offset=inf;';
%            'iAMPA(V)=-gAMPA.*sAMPA(k,:).*(V-EAMPA); EAMPA=0; gAMPA=0;';
%            'iNMDA(V)=-gNMDA.*sNMDA(k,:).*(1.50265./(1+0.33*exp(V./(-16)))).*(V-ENMDA); ENMDA=0; gNMDA=0;';
%            'iGABA(V)=-gGABA.*sGABA(k,:).*(V-EGABA); EGABA=-75; gGABA=0;';
%            'sAMPA=getPoissonGating(bgAMPA,dcAMPA,acAMPA,freqAMPA,phiAMPA,onset,offset,tauAMPA,T,Npop); dcAMPA=0; acAMPA=0; freqAMPA=0; bgAMPA=0; phiAMPA=0; tauAMPA=2;';
%            'sNMDA=getPoissonGating(bgNMDA,dcNMDA,acNMDA,freqNMDA,phiNMDA,onset,offset,tauNMDA,T,Npop); dcNMDA=0; acNMDA=0; freqNMDA=0; bgNMDA=0; phiNMDA=0; tauNMDA=150;';
%            'sGABA=getPoissonGating(bgGABA,dcGABA,acGABA,freqGABA,phiGABA,onset,offset,tauGABA,T,Npop); dcGABA=0; acGABA=0; freqGABA=0; bgGABA=0; phiGABA=0; tauGABA=5';
%            };
