% N-layer PFC model with variable populations and Poisson-based inputs. (script created by JSS on 11-Apr-2016, contact: sherfey@bu.edu)
% wd: ~/models/dynasim/mechanisms
% gist: https://gist.github.com/jsherfey/8010719b0461ddc961829251389a4f0b
% github: ...

% Purpose: shareable basis for extended PFC-related models
% Default model: "L2/3" (E-Ed,FS,RSNP) -> "L5/6" (E-Ed,FS), with Poisson-based inputs.
% Easily extends to more layers with different populations per layer.

% todo: prepare separate function form s.t. (tip: rename this m-file to LaminarPFCdemo.m)
% s=LaminarPFC('size',{[],[8 8],[2 2],[2 0]});
% s=PrepCompetition(s,'inhib_mode','lateral','input_mode','ramp1'); % optional
%   'inhib_mode' {'lateral','feedforward','feedback'}
%   'input_mode' {'ramp1','ramp2','steady','sync1','rhythm1',...}
% d=SimulateModel(s,'vary',vary,solver_options{:});
% PlotData(d,'plot_type','rastergram');

% move to directory storing mechanism files and associated function PFC_cell_specification
mechanism_dir='~/models/dynasim/mechanisms';
cd(mechanism_dir);

% choose cell models for each compartment/population in each layer (arbitrary number of layers)
Esoma_type={'DS00Esoma' ,'DS00Esoma'};  % {layer1_type, layer2_type ...} for regular spiking E-cell soma compartments in each layer
Edend_type={'DS00Edend1','DS00Edend1'}; % {layer1_type, layer2_type ...} for E-cell dendritic compartments in each layer
FS_type   ={'WB96FS'    ,'WB96FS'};     % FS cells are fast spiking PV+ interneurons that target E-cell soma with fast GABA-A inhibition (tau=5ms)
LTS_type  ={'KK08LTSmod','KK08LTSmod'}; % LTS cells are low threshold spiking CB+ interneurons that target E-cell dendrites with longer-lasting GABA-A inhibition (tau=13ms)
% tip: for a complete list of available models and references, execute:
%  help PFC_cell_specification

% NMDA colocalizes w/ PV+ but not CB+:
% ref: DeFelipe, J. (1997). Types of neurons, synaptic connections and chemical characteristics of cells immunoreactive for calbindin-D28K, parvalbumin and calretinin in the neocortex. Journal of chemical neuroanatomy, 14(1), 1-19.

% PY taum is slower in superficial layers (L2/3 25+/-20ms vs L5/6 14+/-7ms)
% L2/3: http://www.neuroelectro.org/neuron/110/
% L5/6: http://www.neuroelectro.org/neuron/111/
% Tip: decrease Rm (or Cm) in deep layer wrt superficial layer

s=ApplyModifications(spec,{'Es->FS','netcon',ones(8,1)})
s=ApplyModifications(spec,{'(Es,Ed)','equations',state_equations})
s=ApplyModifications(spec,{'Es','name','L1Es'})
s=ApplyModifications(spec,{'(L1Es,L2Es)','Cm',1.2})
% directly override parameters set in specification
p=cat(2,p,RemoveKeyval(s.populations(end).parameters,p(1:2:end)));


% set number of cells per population in each layer 
E_size  =[8 8];
FS_size =[2 2]; % typically: .25*(# E-cells)
LTS_size=[2 0];
% tip: exclude a population by setting size to 0

num_layers=length(Esoma_type);

% generic input definition
input_def={'input(t)=-gext.*sext(k,:).*(X-0); monitor input;';
           'sext=getPoissonGating(baseline,dc,ac,freq,phase,onset,offset,2,T,Npop,kernel); kernel=ones(1,Npop)'};
% generic model equations
state_equations=['dV/dt=(@current+input(t))./Cm; Cm=1;' input_def{:}];

% set external input parameters for each population in each layer
% layer- and population-specific input parameters
null_input_parameters={'gext',0,'baseline',0,'dc',0,'ac',0,'freq',0,'onset',0,'offset',inf};
% layer 1
layer=1;
E_input_parameters{layer}  ={'gext',1,'baseline',1,'dc',0,'ac',0,'freq',0,'onset',0,'offset',inf};
FS_input_parameters{layer} =null_input_parameters;
LTS_input_parameters{layer}=null_input_parameters;
% note: inputs to E-cells target Edend if present, else Esoma
% other layers (default: no input)
for layer=2:num_layers
  E_input_parameters{layer}  =null_input_parameters;
  FS_input_parameters{layer} =null_input_parameters;
  LTS_input_parameters{layer}=null_input_parameters;
end

% connections between populations
tauAMPA=2; tauNMDA=150; tauIslow=13; tauIfast=5;
% within-layer max synaptic conductance and connectivity kernels
% E->*: E->E (all-to-all), E->FS (all-to-all), E->LTS (all-to-all)
geeAMPA  =1*ones(1,num_layers); % E->E
geeNMDA  =1*ones(1,num_layers); % E->E, typically: (1 to 1.25)*geeAMPA
gefsAMPA =1*ones(1,num_layers); % E->FS
gefsNMDA =1*ones(1,num_layers); % E->FS
geltsAMPA=1*ones(1,num_layers); % E->LTS
geltsNMDA=1*ones(1,num_layers); % E->LTS
% FS->*: FS->E (all-to-all), FS->LTS (none), FS->FS (none)
gfse   =1*ones(1,num_layers); % FS->E
gfsfs  =0*ones(1,num_layers); % FS->FS
gfslts =0*ones(1,num_layers); % FS->LTS
% LTS->*: LTS->E (all-to-all), LTS->LTS (none), LTS->FS (none)
gltse  =1*ones(1,num_layers); % LTS->E
gltsfs =0*ones(1,num_layers); % LTS->FS
gltslts=0*ones(1,num_layers); % LTS->LTS
for layer=1:num_layers
  Kee{layer}    =ones(E_size(layer)  ,E_size(layer)  ); % N_pre x N_post
  Kefs{layer}   =ones(E_size(layer)  ,FS_size(layer) ); % E->FS
  Kelts{layer}  =ones(E_size(layer)  ,LTS_size(layer)); % E->LTS
  Kfse{layer}   =ones(FS_size(layer) ,E_size(layer)  ); % FS->E
  Kfsfs{layer}  =ones(FS_size(layer) ,FS_size(layer) ); % FS->FS
  Kfslts{layer} =ones(FS_size(layer) ,LTS_size(layer)); % FS->LTS
  Kltse{layer}  =ones(LTS_size(layer),E_size(layer)  ); % LTS->E
  Kltsfs{layer} =ones(LTS_size(layer),FS_size(layer) ); % LTS->FS
  Kltslts{layer}=ones(LTS_size(layer),LTS_size(layer)); % LTS->LTS
end
% feedforward from layer i to layer i+1
% E->E:
geeAMPAff=1*ones(1,num_layers-1);
geeNMDAff=1*ones(1,num_layers-1);
for layer=1:num_layers-1
  Keeff{layer}=ones(E_size(layer),E_size(layer+1)); % E(i)->E(i+1)
end
% normalize kernels by presynaptic connections (s.t. total gating = 1)
for layer=1:num_layers
  K=Kee{layer};     x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kee{layer}=K;
  K=Kefs{layer};    x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kefs{layer}=K;
  K=Kelts{layer};   x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kelts{layer}=K;
  K=Kfse{layer};    x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kfse{layer}=K;
  K=Kfsfs{layer};   x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kfsfs{layer}=K;
  K=Kfslts{layer};  x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kfslts{layer}=K;
  K=Kltse{layer};   x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kltse{layer}=K;
  K=Kltsfs{layer};  x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kltsfs{layer}=K;
  K=Kltslts{layer}; x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Kltslts{layer}=K;
  if layer<num_layers
    K=Keeff{layer}; x=sum(K,1); x(x==0)=1; K=K./repmat(x,[size(K,1) 1]); Keeff{layer}=K;
  end
end

% add populations to model specification
spec=[];
for layer=1:num_layers % loop over layers
  % determine if E-cell inputs go to soma or dendrite (rule: to dendrite if present)
  if ~isempty(Edend_type{layer}) % input to dendrite
    Esoma_input_parameters=null_input_parameters;
    Edend_input_parameters=E_input_parameters{layer};
  else % input to soma
    Esoma_input_parameters=E_input_parameters{layer};
    Edend_input_parameters=null_input_parameters;
  end
  % collect population data for this layer
  names={'Es','Ed','FS','LTS'};
  types={Esoma_type{layer},Edend_type{layer},FS_type{layer},LTS_type{layer}};
  sizes=[E_size(layer),E_size(layer),FS_size(layer),LTS_size(layer)];
  input_parameters={Esoma_input_parameters,Edend_input_parameters,FS_input_parameters{layer},LTS_input_parameters{layer}};
  % add to specification each population in this layer
  for pop=1:length(types)
    if ~isempty(types{pop})
      cellspec                     = PFC_cell_specification(types{pop});
      spec.pops(end+1).name        = sprintf('L%g%s',layer,names{pop});
      spec.pops(end).size          = sizes(pop);
      spec.pops(end).equations     = state_equations;
      spec.pops(end).mechanism_list= cellspec.pops.mechanism_list;
      spec.pops(end).parameters    = cat(2,input_parameters{pop},cellspec.pops.parameters);
    end
  end    
end

% add connections to model specification
for layer=1:num_layers % loop over layers
  % within-layer connections
  if ~isempty(Edend_type{layer}) % connect E and LTS to dendrite
    Etarget=sprintf('L%gEd',layer);
  else % connect E and LTS to soma
    Etarget=sprintf('L%gEs',layer);
  end
  % E->* (iAMPA,iNMDA)
  targets={Etarget,sprintf('L%gFS',layer),sprintf('L%gLTS',layer)};
  Ke={Kee{layer},Kefs{layer},Kelts{layer}};
  geAMPA=[geeAMPA(layer),gefsAMPA(layer),geltsAMPA(layer)];
  geNMDA=[geeNMDA(layer),gefsNMDA(layer),geltsNMDA(layer)];
  for con=1:length(targets)
    spec.cons(end+1).direction = [sprintf('L%gEs',layer) '->' targets{con}];
    spec.cons(end).mechanism_list={'iAMPA','iNMDA'};
    spec.cons(end).parameters={'netcon',Ke{con},'gAMPA',geAMPA(con),'gNMDA',geNMDA(con),'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};
  end
  % FS->* (iGABA with tauIfast; always target Esoma)
  targets={sprintf('L%gEs',layer),sprintf('L%gFS',layer),sprintf('L%gLTS',layer)};
  Kfs={Kfse{layer},Kfsfs{layer},Kfslts{layer}};
  gfs=[gfse(layer),gfsfs(layer),gfslts(layer)];
  for con=1:length(targets)
    spec.cons(end+1).direction = [sprintf('L%gFS',layer) '->' targets{con}];
    spec.cons(end).mechanism_list={'iGABA'};
    spec.cons(end).parameters={'netcon',Kfs{con},'gGABA',gfs(con),'tauGABA',tauIfast};
  end
  % LTS->* (iGABA with tauIslow; target Edend if present else Esoma)
  targets={Etarget,sprintf('L%gFS',layer),sprintf('L%gLTS',layer)};
  Klts={Kltse{layer},Kltsfs{layer},Kltslts{layer}};
  glts=[gltse(layer),gltsfs(layer),gltslts(layer)];
  for con=1:length(targets)
    spec.cons(end+1).direction = [sprintf('L%gLTS',layer) '->' targets{con}];
    spec.cons(end).mechanism_list={'iGABA'};
    spec.cons(end).parameters={'netcon',Klts{con},'gGABA',glts(con),'tauGABA',tauIslow};
  end
  % intercompartmental connections: Ed(i)<->Es(i)
  % ...
  % feeforward connections: E(i)->E(i+1)
  if layer<num_layers
    if ~isempty(Edend_type{layer+1}) % connect to dendrite
      target=sprintf('L%gEd',layer+1);
    else % connect to soma
      target=sprintf('L%gEs',layer+1);
    end
    spec.cons(end+1).direction = [sprintf('L%gEs',layer) '->' target];
    spec.cons(end).mechanism_list={'iAMPA','iNMDA'};
    spec.cons(end).parameters={'netcon',Keeff{layer},'gAMPA',geeAMPAff(layer),'gNMDA',geeNMDAff(layer),'tauAMPA',tauAMPA,'tauNMDA',tauNMDA};
  end
end

% -------------------------------------------------------------------------
% simulation controls
tspan=[0 500];  % [beg end], ms
dt=.01;         % fixed time step, ms
solver='rk2';   % numerical integration method {'rk1','rk2','rk4'}
compile_flag=1; % whether to compile simulation
simulator_options={'tspan',tspan,'solver',solver,'dt',dt,'compile_flag',compile_flag,'verbose_flag',1};
run_sim_flag=0; % whether to run simulation
% (optional) set model modifications to override defaults (same values will be used for all simulations)
mods=[];
% (optional) set model parameters to vary across simulations
vary=[];
% -------------------------------------------------------------------------

if run_sim_flag
  % simulate model
  spec_=ApplyModifications(spec,mods);
  data=SimulateModel(spec_,'vary',vary,simulator_options{:});  
  % plot results
  vars={'L1Ed_V','L1Es_V','L2Ed_V','L2Es_V'};
  PlotData(data,'plot_type','waveform','variable',vars);
  PlotData(data,'plot_type','waveform','variable',{'L1Ed_input','L1Es_input'});
  PlotData(data,'plot_type','power','variable',vars);
  PlotData(data,'plot_type','rastergram','threshold',-20);
  PlotFR(data,'threshold',-20,'bin_size',50)
end
