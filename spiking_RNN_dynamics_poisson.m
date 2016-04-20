% Explore soft WTA and persistence in spiking network with RNN-based
% assemblies and Poisson-based inputs.
% wd: ~/models/dynasim/sandbox
% gist: https://gist.github.com/jsherfey/987453dccf5135631456

% Purpose: simulate Poisson-based inputs (steady or ramping) to two competing 
% assemblies E1 & E2 in E/I network with goal of determining factors controlling 
% persistence and competition.

% move to directory storing mechanism files
mechanism_dir='/home/jason/models/dynasim/mechanisms';
cd(mechanism_dir);

% Model:
%   ramping I1(t) ->  E1(RNN) <-> I <-> E2(RNN)  <- steady I2(t)=I0
%   vary tauI, tauNMDA, gNMDA, gAMPA/gNMDA.

% Observations:
% ** transition to competition: I->E.gGABA=(.5-1) and E->I.gAMPA=(.5-1) ** | E->I.gNMDA=1

cell_model_id=2; % (same model used for all E and I cells)
% 1. Hodgkin-Huxley model
% 2. Minimal PFC mechanisms (Durstewitz, Seamans, Sejnowski 2000; 2002; Durstewitz, Gabriel 2007)
Cm=1.2; gleak=.04; Eleak=-70;
gKDR=33.8; gNaF=86; % Durstewitz 2000
gKDR=50; gNaF=117;  % Durstewitz 2002

E_input_id=1; % (filtered Poisson inputs with different patterns of rate modulation, lambda(t))
% 0. no input (except for noise)
% 1. ramp to assembly 1, constant to all others. (lambda: mean(ramp)=constant=DC)
% 2. ramp to all assemblies (lambda: mean(ramp)=DC, max(ramp)=2*DC=ramp(end))
% 3. constant to all assemblies (lambda: constant=DC)
onset=0; offset=inf; % constrain input to interval [onset,offset]

I_input_id=0; % no inputs to I-cells have been implemented
% 0. no input (except for noise)

network_condition_id=4; % (assembly-interaction conditions implemented by varying connectivity kernels)
% 0. no connections
% 1. assemblies (RNN): without inhibition (explore max assembly activation)
% 2. assemblies (RNN): feedforward inhibition without feedback inhibition
% 3. assemblies (RNN): feedback inhibition without feedforward inhibition
% 4. assemblies (RNN): feedback + lateral inhibition

% solver controls (using rk2)
tspan=[0 500];  % [beg end], ms
dt=.01;         % fixed time step, ms
compile_flag=1; 

% set population sizes
assembly_size=8;                  % # E-cells per assembly {4,8,12,16,20}
num_assemblies=2;                 % # assemblies
Ne=assembly_size*num_assemblies;  % total # E-cells
Ni=.25*Ne;                        % total # I-cells

% define inputs (note: PFC cells get more inputs than other regions: see http://cercor.oxfordjournals.org/content/13/11/1124/F6.large.jpg)
% salva uses for BG MSNs: g=.00025, blext=25kHz, blpfc=2kHz, dcpfc=23kHz, acpfc=.15*dcpfc
% gexte=.0005; enoise=5; DCe=25; ACe=10;
gexte=.0075;      % .1,  .0075, max synaptic conductance for external input to E-cells [uS/cm2]
gexti=.0075;      % .1,  .0075, max synaptic conductance for external input to I-cells [uS/cm2]
enoise=.25;       % .01, .25,   baseline Poisson rate for E-cells [kHz] (1000 inputs at 5Hz)
inoise=.25;       % .01, .25,   baseline Poisson rate for I-cells [kHz]
DCe=2.25;         % .1,  2.25,  stimulus Poisson rate for E-cells [Khz] (eg, 25kHz = 5000 inputs at 5Hz)
ACe=2.25;         % .1,  2.25,  stimulus Poisson rate for I-cells [Khz]
            % g*DC= .011,.0187
kick=1; tau=2;
% kernels, conductance, current (NMDA time constant: http://www.ncbi.nlm.nih.gov/pubmed/7901401)
Ke1=zeros(1,Ne); Ke1(1:assembly_size)=1;    % input kernel for first E-cell assembly
Ke2=zeros(1,Ne); Ke2(assembly_size+1:Ne)=1; % input kernel for all other E-cell assemblies
se1=sprintf('s1=get_input(''poisson'',Npop,T,f1,DC1,AC1,tau,%s,baseline,phase1,kick,ramp_dc_flag1,ramp_ac_flag1,onset,offset); phase1=0',toString(Ke1));
se2=sprintf('s2=get_input(''poisson'',Npop,T,f2,DC2,AC2,tau,%s,baseline,phase2,kick,ramp_dc_flag2,ramp_ac_flag2,onset,offset); phase2=0',toString(Ke2));
si1=sprintf('s1=get_input(''poisson'',Npop,T,f1,DC1,AC1,tau,%s,baseline,phase1,kick,ramp_dc_flag1,ramp_ac_flag1,onset,offset); phase1=0',toString(ones(1,Ni)));
E_input=sprintf('input(t)=-gext.*(s1(k,:)+s2(k,:)).*(X-0); monitor input');
I_input=sprintf('input(t)=-gext.*s1(k,:).*(X-0); monitor input');
% adjustable parameters
switch E_input_id
  case 1 % ramp to assembly 1, constant to all others
    DCe1=2*DCe; ramp_dc_flag1=1; f1=0; ACe1=0; phase1=0; ramp_ac_flag1=0;
    DCe2=1*DCe; ramp_dc_flag2=0; f2=0; ACe2=0; phase2=0; ramp_ac_flag2=0;
  case 2 % ramp to all assemblies
    DCe1=2*DCe; ramp_dc_flag1=1; f1=0; ACe1=0; phase1=0; ramp_ac_flag1=0;
    DCe2=2*DCe; ramp_dc_flag2=1; f2=0; ACe2=0; phase2=0; ramp_ac_flag2=0;
  case 3 % constant to all assemblies
    DCe1=1*DCe; ramp_dc_flag1=0; f1=0; ACe1=0; phase1=0; ramp_ac_flag1=0;
    DCe2=1*DCe; ramp_dc_flag2=0; f2=0; ACe2=0; phase2=0; ramp_ac_flag2=0;
  otherwise % no input
    DCe1=0; ramp_dc_flag1=0; f1=0; ACe1=0; phase1=0; ramp_ac_flag1=0;
    DCe2=0; ramp_dc_flag2=0; f2=0; ACe2=0; phase2=0; ramp_ac_flag2=0;
end
E_input_parameters={'gext',gexte,'baseline',enoise,'DC1',DCe1,'DC2',DCe2,'AC1',ACe1,'AC2',ACe2,'f1',f1,'f2',f2,'ramp_dc_flag1',ramp_dc_flag1,'ramp_ac_flag1',ramp_ac_flag1,'ramp_dc_flag2',ramp_dc_flag2,'ramp_ac_flag2',ramp_ac_flag2,'phase1',phase1,'phase2',phase2,'kick',kick,'tau',tau,'onset',onset,'offset',offset};
switch I_input_id % input kernel for I-cells, [time x cells]    
  otherwise % no input
    DCi1=0;
    DCi2=0;
end
I_input_parameters={'gext',gexti,'baseline',inoise,'DC1',DCi1,'DC2',DCi2,'AC1',0,'AC2',0,'f1',0,'f2',0,'ramp_dc_flag1',0,'ramp_ac_flag1',0,'ramp_dc_flag2',0,'ramp_ac_flag2',0,'phase1',0,'phase2',0,'kick',kick,'tau',tau,'onset',onset,'offset',offset};

% add populations to model specification
equations='dv/dt=(@current+input(t))/Cm; Cm=1';
switch cell_model_id
  case 1 % Hodgkin-Huxley model
    mechanism_list={'iNa','iK','ileak'};
  case 2 % PFC mechanisms (Durstewitz, Seamans, Sejnowski 2000; Durstewitz, Gabriel 2007)
    mechanism_list={'iNaF','iKDR','ileak'};
end
s=[];
s.pops(1).name='E';
s.pops(1).size=Ne;
s.pops(1).equations=sprintf('%s; %s; %s; %s',equations,se1,se2,E_input);
s.pops(1).mechanism_list=mechanism_list;
s.pops(1).parameters={E_input_parameters{:},'gKDR',gKDR,'gNaF',gNaF,'Cm',Cm,'gleak',gleak,'Eleak',Eleak};
s.pops(2).name='I';
s.pops(2).size=Ni;
s.pops(2).equations=sprintf('%s; %s; %s',equations,si1,I_input);
s.pops(2).mechanism_list=mechanism_list;
s.pops(2).parameters={I_input_parameters{:},'gKDR',gKDR,'gNaF',gNaF,'Cm',Cm,'gleak',gleak,'Eleak',Eleak};

% define connectivity
         % gext=.1   .0075
gie=.25;      % .5,  .25
gii=0;        % 0,   0
geiAMPA=.75;  % .25, .75
geiNMDA=.75;  % 1,   .75
geeAMPA=.1;   % 0,   .1
geeNMDA=.1;   % 1,   .1
bwblockee=0; % kernal scaling for connection b/w E-assemblies (i.e., E1->E2 connected with gsyn=gee*bwblockee)
% connectivity kernels:
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
Kii=ones(Ni,Ni); % I->I, [N_pre x N_post] (all-to-all)
% override defaults based on network condition
switch network_condition_id
  case 0 % no connections
    % turn off recurrent connections that define assemblies
    Kee=zeros(Ne,Ne); % E->E, [N_pre x N_post]    
    Kii=zeros(Ni,Ni); % I->I, [N_pre x N_post]
  case 1 % assemblies (RNN): without inhibition (explore max assembly activation)
    % default
  case 2 % assemblies (RNN): feedforward inhibition without feedback inhibition
    % E(1) -> I(1) -> E(2) -> I(2) -> E(3) ...
    % split I-cells into (num_assemblies-1) subsets
    isz=Ni/(num_assemblies-1);
    esz=assembly_size;
    for i=1:num_assemblies-1
      % connect E(i) -> I(i)
      inde=(i-1)*esz+(1:esz); % source E(i)
      indi=(i-1)*isz+(1:isz); % target I(i)
      Kei(inde,indi)=1;
      % connect I(i) -> E(i+1)
      inde=i*esz+(1:esz);     % target E(i+1)
      indi=(i-1)*isz+(1:isz); % source I(i)
      Kie(indi,inde)=1;
    end
  case 3 % assemblies (RNN): feedback inhibition without feedforward inhibition
    % E(1) <-> I(1), E(2) <-> I(2), ...
    % split I-cells into (num_assemblies) subsets
    isz=Ni/num_assemblies;
    esz=assembly_size;
    for i=1:num_assemblies
      inde=(i-1)*esz+(1:esz); % E(i)
      indi=(i-1)*isz+(1:isz); % I(i)
      % connect E(i) -> I(i)
      Kei(inde,indi)=1;
      % connect I(i) -> E(i)
      Kie(indi,inde)=1;
    end
  case 4 % assemblies (RNN): feedback + lateral inhibition
    % {E(1),E(2),...} <-> I.  e.g.) E(1) <-> I <-> E(2)
    Kie=ones(Ni,Ne); % I->E, [N_pre x N_post]
    Kei=ones(Ne,Ni); % E->I, [N_pre x N_post]
end
% normalize kernels by number of presynaptic connections
Kee=Kee./repmat(max(1,sum(Kee,1)),[size(Kee,1) 1]);
Kei=Kei./repmat(max(1,sum(Kei,1)),[size(Kei,1) 1]);
Kie=Kie./repmat(max(1,sum(Kie,1)),[size(Kie,1) 1]);
Kii=Kii./repmat(max(1,sum(Kii,1)),[size(Kii,1) 1]);

% synaptic time constants
% HPC NMDA time constant: http://www.ncbi.nlm.nih.gov/pubmed/7901401
% "The average fast (tau f) and slow (tau s) time constants of decay were, respectively, 66.5 and 353.9 ms in PCs, and 34.4 and 212.5 ms in M-INs."
% PFC NMDA time constant: http://www.pnas.org/content/105/43/16791.full
tauNMDAe=150;
tauNMDAi=150;
tauGABA=5;

% add connections to model specification
s.cons(1).direction='E->I';
s.cons(1).mechanism_list={'iAMPA','iNMDA'};
s.cons(1).parameters={'netcon',Kei,'gAMPA',geiAMPA,'gNMDA',geiNMDA,'tauNMDA',tauNMDAi};
s.cons(2).direction='E->E';
s.cons(2).mechanism_list={'iAMPA','iNMDA'};
s.cons(2).parameters={'netcon',Kee,'gAMPA',geeAMPA,'gNMDA',geeNMDA,'tauNMDA',tauNMDAe};
s.cons(3).direction='I->E';
s.cons(3).mechanism_list='iGABA';
s.cons(3).parameters={'netcon',Kie,'gGABA',gie,'tauGABA',tauGABA};
s.cons(4).direction='I->I';
s.cons(4).mechanism_list='iGABA';
s.cons(4).parameters={'netcon',Kii,'gGABA',gii,'tauGABA',tauGABA};

% simulate model
mods=[]; vary=[];
% mods={'E->E','gAMPA',.1;   'E->E','gNMDA',.1; 
%       'E->I','gAMPA',.75; 'E->I','gNMDA',.75;
%       'I->E','gGABA',.25;  'I->I','gGABA',0;
%       '(E->E,E->I)','tauNMDA',150;'(I->I,I->E)','tauGABA',5'};

% vary={'E->I','gAMPA',[.25 .5 .75];'I->E','gGABA',[.25 .5 .75]};
% vary={'E->E','gAMPA',0:.05:.2;'E->E','gNMDA',0:.05:.2};
% vary={'I->I','gGABA',[0 .25 .5 .75]};
% vary={'E','rep',1:5};
% vary={'E->E','gAMPA',0:.05:.2;'E->E','gNMDA',0:.05:.2};
% vary={'E->E','(gAMPA,gNMDA)',[0];'E','(DC1,DC2)',[3:6]};
  
vary={'E->I','gNMDA',0:.25:1;'E->I','gAMPA',0:.25:1;'I->E','gGABA',.25};
% vary={'(E->E,E->I)','tauNMDA',[150 250 350 450];'(I->I,I->E)','tauGABA',[5:5:30]};
% vary={'E','gleak',[.01 .04 .1];'I','gleak',[.01 .04 .1]};
% vary={'E','Cm',[.75 1 1.25];'I','Cm',[.75 1 1.25]};

% vary={'(E,I)','gext',.001:.002:.01;'E','(DC1,DC2)',1:2:10;'(E,I)','baseline',0};
% vary={'(E,I)','gext',.006:.0005:.008;'E','(DC1,DC2)',2:.5:4;'(E,I)','baseline',0};
% vary={'(E,I)','gext',.008;'E','(DC1,DC2)',[2 2.5];'(E,I)','baseline',[0 .25 .5]};
% vary={'(E,I)','gext',.0075;'E','(DC1,DC2)',[1.75 2 2.25];'(E,I)','baseline',[.25 .5]};
% vary={'I','baseline',[0:.05:.4]};

solver_options={'tspan',tspan,'solver','rk2','dt',dt,'compile_flag',compile_flag,'verbose_flag',1};
data=SimulateModel(ApplyModifications(s,mods),'vary',vary,solver_options{:});

% plot results
%PlotData(data,'plot_type','waveform','variable','E_v');
%PlotData(data,'plot_type','power','variable','E_v');
% PlotData(data,'plot_type','waveform','variable','input');
PlotData(data,'plot_type','rastergram','threshold',-20);
PlotFR(data,'threshold',-20,'bin_size',50)

%% analysis of single simulation
k=length(data);         % simulation index
threshold=-20;          % mV, voltage threshold for spike detection
kwidth=20;              % ms, width of gaussian for kernel regression
Ts=1;                   % ms, set to this effective time step for rate process before regression
toilim=[20 tspan(2)];   % [beg end], limits on times to analyze
roi1=1:assembly_size;   % assembly roi 1
roi2=assembly_size+1:Ne;% assembly roi 2
toi=data(k).time>=toilim(1)&data(k).time<=toilim(2);
t=data(k).time(toi);
V=data(k).E_v(toi,:);
I=data(k).E_input(toi,:);
raster=computeRaster(t,V,threshold);
  % raster(:,1) -> spike times
  % raster(:,2) -> cell index for each spike
% assembly outputs
[iFR1,itime]=NWgaussKernelRegr(t,raster,roi1,kwidth,Ts); iFR1=1000*iFR1;
[iFR2,itime]=NWgaussKernelRegr(t,raster,roi2,kwidth,Ts); iFR2=1000*iFR2;
% assembly inputs
I1=mean(I(:,roi1),2);
I2=mean(I(:,roi2),2);
y=ksr(itime,resample(I1,1,Ts/dt,2),kwidth);
iI1=y.f; % guassian-regressed input to assembly 1 resampled for plotting vs. inst. rate
y=ksr(itime,resample(I2,1,Ts/dt,2),kwidth);
iI2=y.f; % guassian-regressed input to assembly 2 resampled for plotting vs. inst. rate
% analysis plots
figure('position',[285 370 1250 420]); 
ylims=[min([iFR1(:);iFR2(:)]) max([iFR1(:);iFR2(:)])];
subplot(2,3,1); plot(itime,iFR1,'b'); ylabel('iFR1(t) [Hz]'); title('assembly 1'); xlim(toilim); ylim(ylims);
subplot(2,3,2); plot(itime,iFR2,'r'); ylabel('iFR2(t) [Hz]'); title('assembly 2'); xlim(toilim); ylim(ylims);
subplot(2,3,3); plot(itime,iFR1,'b',itime,iFR2,'r'); title('outputs'); xlim(toilim); ylim(ylims);
legend('roi1','roi2','Location','NorthWest')
ylims=[min([iI1(:);iI2(:)]) max([iI1(:);iI2(:)])];
subplot(2,3,4); plot(itime,iI1,'b'); ylabel('<I1(t)>'); xlim(toilim); ylim(ylims);
subplot(2,3,5); plot(itime,iI2,'r'); ylabel('<I2(t)>'); xlabel('time (ms)'); xlim(toilim); ylim(ylims);
subplot(2,3,6); plot(itime,iI1,'b',itime,iI2,'r'); title('inputs'); xlim(toilim); ylim(ylims);
legend('I1','I2','Location','NorthWest');
% Relative activity and Competition
figure('position',[380 12 1060 275]); 
dr=(iFR1-iFR2)./(iFR1+iFR2);
subplot(1,2,1); plot(itime,abs(dr),'r--',itime,dr,'b'); line(xlim,[0 0]); ylim([-1 1]); xlim(toilim);
xlabel('time (ms)'); ylabel('rel activity'); legend('abs(dr)','dr','Location','SouthEast')
subplot(1,2,2); plot(iI1-iI2,dr); line(xlim,[0 0]); ylim([-1 1]); 
xlabel('rel input'); ylabel('rel activity');

% figure
% subplot(2,1,1); plot(t,I1,'b'); ylim([0 5]); ylabel('<I1(t)>'); xlim(toilim);
% subplot(2,1,2); plot(t,I2,'r'); ylim([0 5]); ylabel('<I2(t)>'); xlabel('time (ms)'); xlim(toilim);

% Spectral analysis of instantaneous firing rates
NFFT=2^(nextpow2(length(itime)-1)-2);
WINDOW=2^(nextpow2(NFFT-1)-3);
[Y,F,T]=spectrogram(iFR1,WINDOW,[],NFFT,1/(Ts/1000));
P1=(abs(Y)+eps).^2;
[Y,F,T]=spectrogram(iFR2,WINDOW,[],NFFT,1/(Ts/1000));
P2=(abs(Y)+eps).^2;
mu=repmat(mean([P1 P2],2),[1 size(P1,2)]);
sd=repmat(std([P1 P2],[],2),[1 size(P1,2)]);
P1z=(P1-mu)./sd;
P2z=(P2-mu)./sd;
% plot spectra
figure('position',[1450 370 485 400]); 
fsel=F>(0)&F<(200); tmp=[P1(fsel,:);P2(fsel,:)]; 
clims=prctile(tmp(:),[0 97]);
subplot(2,2,1); imagesc(T,F(fsel),P1(fsel,:)); title('TFR(iFR1) assembly 1'); ylabel('freq (Hz)'); xlabel('time (sec)'); axis xy; axis tight; caxis(clims); colorbar
subplot(2,2,3); imagesc(T,F(fsel),P2(fsel,:)); title('TFR(iFR2) assembly 2'); ylabel('freq (Hz)'); xlabel('time (sec)'); axis xy; axis tight; caxis(clims); colorbar
clims=[-3 3];
subplot(2,2,2); imagesc(T,F(fsel),P1z(fsel,:)); title('TFRz(iFR1) assembly 1'); ylabel('freq (Hz)'); xlabel('time (sec)'); axis xy; axis tight; caxis(clims); colorbar
subplot(2,2,4); imagesc(T,F(fsel),P2z(fsel,:)); title('TFRz(iFR2) assembly 2'); ylabel('freq (Hz)'); xlabel('time (sec)'); axis xy; axis tight; caxis(clims); colorbar

%%

run_experiment_flag=0;
if run_experiment_flag
  % ??? (gleak,Cm),(tauGABA,tauNMDA),(gNMDA,gNMDA/gAMPA) ???
  vars={'E_v','E_input','I_input','time'};
  %               **             **    **           **
  gie_=.25;    % .25       0     0    .25     .25   .25
  gee_=.1;     % .1        .1    .1    .1      0    .1
  gei_=1.25;   % .75       .75   .75   .75   .75    1
  DC_=2.75;    % 2.25      2.25  1    2.25    3     2.25
  AC_=2.75;    % 2.25      2.25  1    1.125   3     2.25
  f_=0:10:90; % 0:10:90
  % 1) tuning curve
  f=f_; DC=DC_; AC=AC_; nrep=3; % f=0:5:45;
  mods={'E','(ramp_dc_flag1,ramp_dc_flag2,ramp_ac_flag1,ramp_ac_flag2)',0;
        'E','(DC1,DC2)',DC;'E','(AC1,AC2)',AC;'I->E','gGABA',gie_;'E->E','(gAMPA,gNMDA)',gee_;'E->I','(gAMPA,gNMDA)',gei_};
  vary={'E','(f1,f2)',f;'E','rep',1:nrep};
  data=SimulateModel(ApplyModifications(s,mods),'vary',vary,solver_options{:});
  data=rmfield(data,setdiff(data(1).labels,vars)); [data.labels]=deal(vars);
  sims=1:(2*nrep):length(data);
  PlotData(data(sims),'plot_type','rastergram','threshold',-20);
  PlotData(data(sims),'plot_type','waveform','variable','input');
  PlotData(data(sims),'plot_type','power','variable','E_v')
  PlotData(data(sims),'plot_type','power','variable','input')  
  clear stats1
  for i=1:length(data)
    stats1(i)=CalcSpikeSync(data(i),'ROI_pairs',{'E_v',[0 .5],'E_v',[.5 1]});
  end
  xc1=arrayfun(@(x)x.pairs.xcmax_pops,stats1);
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
  if 0
    figure('position',[750 135 460 740]); 
    subplot(4,1,1); plot(f,FR,'ko-'); xlabel('fext [Hz]'); ylabel('<FR>'); line([fc fc],ylim);
    subplot(4,1,2); plot(f,fMUAmu,'ko-'); xlabel('fext [Hz]'); ylabel('<fMUA>'); ylim([min(f) max(f)]); line([fc fc],ylim);
    subplot(4,1,3); plot(FR,fMUAmu,'go-'); xlabel('<FR>'); ylabel('<fMUA>');
    subplot(4,1,4); plot(f,xc1mu,'ko-'); xlabel('fext [Hz]'); ylabel('spkcoh(E1,E2)');
  end
  
  % 2) rhythmic vs equal-strength nonrhythmic
  f1=f; f2=0; nrep=3;
  DConly=DC; t=0:1e-5:1; % one cycle
  for i=1:length(DC)
    DConly(i)=sum(max(0,AC*sin(2*pi*t)+DC(i)))/length(t);
  end
  mods={'E','(ramp_dc_flag1,ramp_dc_flag2,ramp_ac_flag1,ramp_ac_flag2)',0;
        'E','DC1',DC;'E','DC2',DConly;'E','(AC1,AC2)',AC;'I->E','gGABA',gie_;'E->E','(gAMPA,gNMDA)',gee_;'E->I','(gAMPA,gNMDA)',gei_};
  vary={'E','f1',f1;'E','f2',f2;'E','rep',1:nrep};
  data2=SimulateModel(ApplyModifications(s,mods),'vary',vary,solver_options{:});  
  data2=rmfield(data2,setdiff(data2(1).labels,vars)); [data2.labels]=deal(vars);
  clear stats2
  for i=1:length(data2)
    stats2(i)=CalcSpikeSync(data2(i),'ROI_pairs',{'E_v',[0 .5],'E_v',[.5 1]});
  end
  xc2=arrayfun(@(x)x.pairs.xcmax_pops,stats2);
  dr=arrayfun(@(x)x.pairs.dNsumN,stats2);
  drmu2=zeros(size(f)); drse2=zeros(size(f));
  xc2mu=zeros(size(f)); xc2se=zeros(size(f));
  for i=1:length(f)
    sel=[data2.E_f1]==f(i);
    drmu2(i)=mean(dr(sel));
    drse2(i)=std(dr(sel));%/sqrt(nrep);
    xc2mu(i)=mean(xc2(sel));
    xc2se(i)=std(xc2(sel));%/sqrt(nrep)
  end
  if 0
    figure('position',[690 430 290 420]); 
    plot_CI(f,drmu2,drse2,'b','o'); xlabel('fext [Hz]'); ylabel('rel activity: E1(f1)-E2(0)'); ylim([-1 1]); line(xlim,[0 0]); line([fc fc],ylim);
  end
  
  % 3) rhythmic vs rhythmic  
  f1=f; f2=fc; nrep=3;
  mods={'E','(ramp_dc_flag1,ramp_dc_flag2,ramp_ac_flag1,ramp_ac_flag2)',0;
        'E','(DC1,DC2)',DC;'E','(AC1,AC2)',AC;'I->E','gGABA',gie_;'E->E','(gAMPA,gNMDA)',gee_;'E->I','(gAMPA,gNMDA)',gei_};
  vary={'E','f1',f1;'E','f2',f2;'E','rep',1:nrep};
  data3=SimulateModel(ApplyModifications(s,mods),'vary',vary,solver_options{:});  
  data3=rmfield(data3,setdiff(data3(1).labels,vars)); [data3.labels]=deal(vars);
  clear stats3
  for i=1:length(data3)
    stats3(i)=CalcSpikeSync(data3(i),'ROI_pairs',{'E_v',[0 .5],'E_v',[.5 1]});
  end
  xc3=arrayfun(@(x)x.pairs.xcmax_pops,stats3);
  dr=arrayfun(@(x)x.pairs.dNsumN,stats3);
  drmu3=zeros(size(f)); drse3=zeros(size(f));
  xc3mu=zeros(size(f)); xc3se=zeros(size(f));
  for i=1:length(f)
    sel=[data3.E_f1]==f(i);
    drmu3(i)=mean(dr(sel));
    drse3(i)=std(dr(sel));%/sqrt(nrep);
    xc3mu(i)=mean(xc3(sel));
    xc3se(i)=std(xc3(sel));%/sqrt(nrep)
  end
  dfc=abs(f1-fc)-abs(f2-fc);
  if 0
    figure('position',[1300 430 560 420]); 
    subplot(1,2,1); plot_CI(f,drmu3,drse3,'r','o'); xlabel('f1 [Hz] | f2=fc'); ylabel('rel activity: E1(f1)-E2(fc)'); ylim([-1 1]); line(xlim,[0 0]); line([fc fc],ylim);
    subplot(1,2,2); plot_CI(dfc,drmu3,drse3,'r','o'); xlabel('|f1-fc|-|f2-fc|'); ylabel('rel activity'); ylim([-1 1]); line(xlim,[0 0]);
  end
  
  % calculate biases
  resbias2=(FR-FR(f==0))./(FR+FR(f==0));   % resonance-mediated bias (rhythmic vs nonrhythmic)
  inhbias2=drmu2'-resbias2;                 % inhibition-mediated bias (rhythmic vs nonrhythmic)
  resbias3=(FR-FR(f==fc))./(FR+FR(f==fc)); % resonance-mediated bias (rhythmic vs rhythmic)
  inhbias3=drmu3'-resbias3;                 % inhibition-mediated bias (rhythmic vs rhythmic)
  
  % ----- summary plots ----- %
  figure('position',[680 25 705 950]);
  gi=data(1).model.parameters.E_I_iGABA_gGABA;
  gee=data(1).model.parameters.E_E_iAMPA_gAMPA;
  gei=data(1).model.parameters.I_E_iAMPA_gAMPA;
  ac=data(1).model.parameters.E_AC1;
  dc=data(1).model.parameters.E_DC1;
  subplot(2,2,1); plot_CI(f,FRmu,FRse,'g','o'); hold on
  plot(f,FR,'ko-','linewidth',3); xlabel('fext [Hz]'); ylabel('<FR> [spk/s]'); line([fc fc],ylim,'color','k'); axis tight; 
  str=sprintf('network tuning profile\n(g(i->e)=%g,gee=%g,gei=%g,ac=%g,dc=%g,ac/dc=%g)',gi,gee,gei,ac,dc,ac/dc); title(str);
  subplot(2,2,3); plot_CI(f,fMUAmu,fMUAsd,'g','o'); xlabel('fext [Hz]'); ylabel('<fMUA> [Hz]'); ylim([min(f) max(f)]); line([fc fc],ylim); line(xlim,[fc fc]); title('network rhythmicity');
  subplot(2,2,2); plot_CI(f,drmu2,drse2,'b','o'); xlabel('f1 [Hz] | f2=0'); ylabel('rel activity: E1(f1)-E2(0)'); ylim([-1 1]); line(xlim,[0 0]); line([fc fc],ylim,'color','k'); title('rhythmic vs. nonrhythmic');
  hold on; plot(f,resbias2,'go-','linewidth',3); %plot(f,inhbias2,'ro-');
  subplot(2,2,4); plot_CI(f,drmu3,drse3,'b','o'); xlabel('f1 [Hz] | f2=fc'); ylabel('rel activity: E1(f1)-E2(fc)'); ylim([-1 1]); line(xlim,[0 0]); line([fc fc],ylim,'color','k'); title('rhythmic vs. rhythmic');
  hold on; plot(f,resbias3,'go-','linewidth',3); %plot(f,inhbias3,'ro-');
  %file='fig1_gie.25_gee.1_gei.75_ac1.25_dc2.25_f0-90Hz_t0-500ms_summary'; set(gcf,'PaperPositionMode','auto'); print(gcf,[file '.jpg'],'-djpeg'); print(gcf,[file '.eps'],'-depsc'); saveas(gcf,[file '.fig']); plot2svg([file '.svg'],gcf);
  
  figure('position',[1400 50 333 913]); tmp=[xc1mu xc2mu xc3mu]; ylims=[min(tmp)-.1 max(tmp)+.1];%ylims=[0 .5];
  subplot(3,1,1); plot_CI(f,xc1mu,xc1se,'g','o'); xlabel('fext [Hz]'); ylabel('spkcoh(E1,E2)'); ylim(ylims); line([fc fc],ylim,'color','k'); title(str)
  subplot(3,1,2); plot_CI(f,xc2mu,xc2se,'b','o'); xlabel('f1 [Hz] | f2=0'); ylabel('spkcoh(E1,E2)'); title('rhythmic vs. nonrhythmic'); ylim(ylims); line([fc fc],ylim,'color','k');
  subplot(3,1,3); plot_CI(f,xc3mu,xc3se,'b','o'); xlabel('f1 [Hz] | f2=fc'); ylabel('spkcoh(E1,E2)'); title('rhythmic vs. rhythmic'); ylim(ylims); line([fc fc],ylim,'color','k');
  
%   figure('position',[680 25 330 950]);
%   gi=data(1).model.parameters.E_I_iGABA_gGABA;
%   ac=data(1).model.parameters.E_AC1;
%   dc=data(1).model.parameters.E_DC1;
%   subplot(3,1,1); plot_CI(f,FRmu,FRse,'g','o'); hold on
%   plot(f,FR,'ko-','linewidth',3); xlabel('fext [Hz]'); ylabel('<FR>'); line([fc fc],ylim,'color','k'); axis tight; 
%   str=sprintf('network tuning profile\n(g(i->e)=%g,dc=%g,ac=%g,ac/dc=%g)',gi,ac,dc,ac/dc); title(str);
%   subplot(3,1,2); plot_CI(f,drmu2,drse2,'b','o'); xlabel('f1 [Hz] | f2=0'); ylabel('rel activity: E1(f1)-E2(0)'); ylim([-1 1]); line(xlim,[0 0]); line([fc fc],ylim,'color','k'); title('rhythmic vs. nonrhythmic');
%   subplot(3,1,3); plot_CI(f,drmu3,drse3,'r','o'); xlabel('f1 [Hz] | f2=fc'); ylabel('rel activity: E1(f1)-E2(fc)'); ylim([-1 1]); line(xlim,[0 0]); line([fc fc],ylim,'color','k'); title('rhythmic vs. rhythmic');
end


%% Useful parameter spaces
% GIVEN: cell_model_id=2; E_input_id=1;

% noisy f/I response: Probe FR = f(DC, baseline) (network_condition_id=0)
% vary={'(E,I)','gext',[.0005 .001 .0015 .002];'E','(DC1,DC2)',10:10:50;'(E,I)','baseline',5};
% vary={'(E,I)','gext',.0015;'E','(DC1,DC2)',5:12;'(E,I)','baseline',0:5};
% vary={'(E,I)','gext',.001:.002:.01;'E','(DC1,DC2)',1:2:10;'(E,I)','baseline',0};
% vary={'(E,I)','gext',.006:.0005:.008;'E','(DC1,DC2)',2:.5:4;'(E,I)','baseline',0};
% vary={'(E,I)','gext',.0075;'E','(DC1,DC2)',[0 .5 1 1.5];'(E,I)','baseline',[0 .5 1 1.5]};
% vary={'(E,I)','gext',[.01 .1 1];'E','(DC1,DC2)',[0 .01 .1 1];'(E,I)','baseline',0};
% vary={'(E,I)','gext',.1;'E','(DC1,DC2)',0:.05:.15;'(E,I)','baseline',0:.01:.04};
% vary={'(E,I)','gext',.1;'E','(DC1,DC2)',.1;'(E,I)','baseline',.01};

% examine impact of recurrent excitation (network_condition_id=1)
% vary={'E->E','gAMPA',0:.05:.2;'E->E','gNMDA',0:.5:2};

% examine impact of feedforward inhibition (network_condition_id=2,3,4)
% vary={'I->E','gGABA',[.25 .5 .75];'E->I','gAMPA',[.25 .5 .75];'E->I','gNMDA',[0 1]};
% vary={'I->I','gGABA',[0 .25 .5 .75]};

% examine effects of biophysical parameters on spiking: (network_condition_id=0)
% vary={'E','gKDR',[10 20 30 40 50];'E','gNaF',[75 85 95 105 115];'E','Iamp',[10];'E','noise',[4];'E->E','(gNMDA,gAMPA)',0};

if 0
network_condition_id=4; cell_model_id=2; assembly_size=4; num_assemblies=2;
E_input_id=1; onset=0; offset=inf; I_input_id=0;  
vary={'E','Iamp',0:3;'E','noise',0:.5:2};
print(gcf,fullfile(outdir,'competition_Eid2_inpid1_netid4_Enoise1_Eamp20_Inoise2_vary-Eamp-Enoise.jpg'),'-djpeg')

  
network_condition_id=4; cell_model_id=2; assembly_size=4; num_assemblies=2;
E_input_id=1; onset=0; offset=inf; I_input_id=0;  
vary={'E->I','gAMPA',0:.5:2;'I->E','gGABA',0:.5:2};
print(gcf,fullfile(outdir,'competition_Eid2_inpid1_netid4_Enoise1_Eamp20_Inoise2_vary-gAMPA-gGABA.jpg'),'-djpeg')

end


% ramp_dc_flag1=0; ramp_ac_flag1=1; f1=10; DC1=1; AC1=1; tau=2; baseline=0; phase1=0; kick=1; 
% s1=get_input('poisson',5,(tspan(1):dt:tspan(2))',f1,DC1,AC1,tau,ones(1,5),baseline,phase1,kick,ramp_dc_flag1,ramp_ac_flag1);
% PlotData(s1);
% ramp_dc_flag1=0; ramp_ac_flag1=0; f1=10; DC1=1; AC1=1; tau=2; baseline=0; phase1=0; kick=1; 
% ramp_dc_flag1=0; ramp_ac_flag1=0; f1=10; DC1=.1; AC1=.1; tau=2; baseline=.01; phase1=0; kick=1; 
% s1=get_input('poisson',5,(tspan(1):dt:tspan(2))',f1,DC1,AC1,tau,ones(1,5),baseline,phase1,kick,ramp_dc_flag1,ramp_ac_flag1);
% PlotData(s1);

