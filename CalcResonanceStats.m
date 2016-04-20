function [stats,data]=CalcPopulationStats(data,varargin)
% Characterize population responses and effects on them of varying model
% components.
% 
% options=CheckOptions(varargin,{...
%   'bin_size',30,[],... % ms
%   'bin_shift',10,[],... % ms
%   'repetition_parameter','realization',[],...  % e.g., 'realization'
%   'sweep_parameter',[],[],...                  % e.g., 'f'
%   'plot_flag',1,[],...
%   },false);
% 
% Output:
% stats.(FR_var).
%     [time_FR]           = num_bins         x 1
%     [pop_FR]            = num_bins         x num_sims, average over cells
%     [pop_FR_mean]       = 1                x num_sims, average over cells and time
%     [parameters]        = num_varied       x num_sims
%     [varied]            = num_varied       x 1
% stats.(FR_var).repetition_sets.
%     [pop_FR_mu]         = num_bins         x num_rep_sets(=num_sims/num_repetitions)
%     [pop_FR_sd]         = num_bins         x num_rep_sets
%     [parameters]        = (num_varied-1)   x num_rep_sets, excludes repetition_parameter
%     [varied]            = (num_varied-1)   x 1
% stats.(FR_var).sweep_sets.
%     [sweep_pop_FR]      = num_sweep_values x num_sweeps(=num_sims/num_sweep_values), note: is actually pop_FR_mean reorganized by sweep parameter
%     [sweep_pop_FR_max]  = 1                x num_sweeps
%     [parameters]        = (num_varied-1)   x num_sweeps, excludes sweep_parameter
%     [varied]            = (num_varied-1)   x 1
% stats.(FR_var).sweep_sets.repetition_sets.
%     [sweep_pop_FR_mu]   = num_sweep_values x (num_sweeps/num_repetitions)
%     [sweep_pop_FR_sd]   = num_sweep_values x (num_sweeps/num_repetitions)
%     [sweep_pop_FR_max_mu] = 1              x (num_sweeps/num_repetitions)
%     [sweep_pop_FR_max_sd] = 1              x (num_sweeps/num_repetitions)
%     [parameters]        = (num_varied-2)   x (num_sweeps/num_repetitions), excludes sweep_parameter & repetition_parameter
%     [varied]            = (num_varied-2)   x 1
% 
% stats.(MUA_var).
%     [Pxx]               = num_freqs       x num_sims
%     [frequency]         = num_freqs       x 1
%     [PeakFreq]          = 1               x num_sims
%     [PeakArea]          = 1               x num_sims
%     [parameters]        = num_varied      x num_sims
%     [varied]            = num_varied      x 1
% stats.(MUA_var).repetition_sets.
%     [Pxx_mu]            = num_freqs       x num_rep_sets(=num_sims/num_repetitions)
%     [Pxx_sd]            = num_freqs       x num_rep_sets
%     [PeakFreq_mu]       = 1               x num_rep_sets
%     [PeakFreq_sd]       = 1               x num_rep_sets
%     [parameters]        = (num_varied-1)  x num_rep_sets, excludes repetition_parameter
%     [varied]            = (num_varied-1)  x 1
% stats.(MUA_var).sweep_sets.
% 
% stats.(SUA_var).
%     [Pxx]               = num_freqs       x num_sims, averaged over cells
%     [frequency]         = num_freqs       x 1
%     [PeakFreq]          = 1               x num_sims, averaged over cells
%     [PeakArea]          = 1               x num_sims, averaged over cells
%     [parameters]        = num_varied      x num_sims
%     [varied]            = num_varied      x 1
% stats.(SUA_var).repetition_sets.
%     (same as stats.(MUA_var).repetition_sets.)
% stats.(SUA_var).sweep_sets.
%     (same as stats.(MUA_var).sweep_sets.)

% Check options
options=CheckOptions(varargin,{...
  'bin_size',30,[],... % ms
  'bin_shift',10,[],... % ms
  'repetition_parameter','realization',[],...  % e.g., 'realization'
  'sweep_parameter',[],[],...                  % e.g., 'f'
  'plot_flag',1,[],...
  },false);

% Check data
% data=CheckData(data);

stats=[];

%% 1.0 Collect information on what was varied

% number of data sets (one per simulation)
num_sims=length(data); 

% what was varied
if isfield(data,'varied')
  varied=data(1).varied'; % name of parameters varied
  num_varied=length(varied); % number of model components varied across simulations
  params_all=zeros(num_varied,num_sims); % values for each simulation
  % loop over varied components and collect values
  for i=1:num_varied
    params_all(i,:)=[data.(varied{i})]; % values for each simulation
      % each row of params_all contains the varied values for a single simulation.
      % [num_sims x num_varied]
  end
  % check for presence of a parameter indexing multiple realizations of each simulation
  if ~ismember(options.repetition_parameter,varied)
    options.repetition_parameter=[];
    stats.num_repetitions=1;
  else
    repetition_varied_idx=ismember(varied,options.repetition_parameter);
    num_repetition_values=length(unique(params_all(repetition_varied_idx,:)));
    stats.num_repetitions=num_repetition_values;
  end
  % check for presence of a parameter over which others were swept
  if ~isempty(options.sweep_parameter) && ~ismember(options.sweep_parameter,varied)
    options.sweep_parameter=[];
  elseif ~isempty(options.sweep_parameter)
    sweep_varied_idx=ismember(varied,options.sweep_parameter);
    sweep_values=unique(params_all(sweep_varied_idx,:),'stable');
    num_sweep_values=length(sweep_values);
  end
else
  varied={};
  num_varied=0;
  params_all=nan;
  options.repetition_parameter=[];
  options.sweep_parameter=[];
end

%% 2.0 Firing rate analysis
data=CalcFR(data,varargin{:});

% get list of firing rate variables
FR_vars=data(1).results(~cellfun(@isempty,regexp(data(1).results,'.*_FR$')));
FR_vars=setdiff(FR_vars,'time_FR');
stats.FR_variables=FR_vars;
stats.time_FR=data(1).time_FR;
num_bins=length(stats.time_FR);

% stats.(FR_field).
%     [time_FR]           = num_bins         x 1
%     [pop_FR]            = num_bins         x num_sims, average over cells
%     [pop_FR_mean]       = 1                x num_sims, average over cells and time
%     [parameters]        = num_varied       x num_sims
%     [varied]            = num_varied       x 1
% stats.(FR_field).repetition_sets.
%     [pop_FR_mu]         = num_bins         x num_rep_sets(=num_sims/num_repetitions)
%     [pop_FR_sd]         = num_bins         x num_rep_sets
%     [parameters]        = (num_varied-1)   x num_rep_sets, excludes repetition_parameter
%     [varied]            = (num_varied-1)   x 1
% stats.(FR_field).sweep_sets.
%     [sweep_pop_FR]      = num_sweep_values x num_sweeps(=num_sims/num_sweep_values), note: is actually pop_FR_mean reorganized by sweep parameter
%     [sweep_pop_FR_max]  = 1                x num_sweeps
%     [parameters]        = (num_varied-1)   x num_sweeps, excludes sweep_parameter
%     [varied]            = (num_varied-1)   x 1
% stats.(FR_field).sweep_sets.repetition_sets.
%     [sweep_pop_FR_mu]   = num_sweep_values x (num_sweeps/num_repetitions)
%     [sweep_pop_FR_sd]   = num_sweep_values x (num_sweeps/num_repetitions)
%     [sweep_pop_FR_max_mu] = 1              x (num_sweeps/num_repetitions)
%     [sweep_pop_FR_max_sd] = 1              x (num_sweeps/num_repetitions)
%     [parameters]        = (num_varied-2)   x (num_sweeps/num_repetitions), excludes sweep_parameter & repetition_parameter
%     [varied]            = (num_varied-2)   x 1

% 2.1 population analysis
% calculate population means for each FR variable
for v=1:length(FR_vars)
  var=FR_vars{v};
  pop_FR=zeros(num_bins,num_sims);
  for i=1:num_sims
    pop_FR(:,i)=mean(data(i).(var),2); % average over cells
  end
  % store results for this population variable
  stats.(var).pop_FR=pop_FR;              % num_bins x num_sims
  stats.(var).pop_FR_mean=mean(pop_FR,1); % average over time
  stats.(var).parameters=params_all;
  stats.(var).varied=varied;
end

% 2.2 calc stats over multiple realizations
if ~isempty(options.repetition_parameter)
  % map repetitions of the same parameter set
  varied_repetitions=varied(~repetition_varied_idx);
  params_repetitions=unique(params_all(~repetition_varied_idx,:)','rows','stable')';
  num_repetition_sets=size(params_repetitions,2);
  [~,LOCB]=ismember(params_all(~repetition_varied_idx,:)',params_repetitions','rows');
  % define data set indices for each set of repetitions
  repetition_set_indices=nan(num_repetition_sets,num_repetition_values);
  for i=1:num_repetition_sets
    repetition_set_indices(i,:)=find(LOCB==i);
  end
  % calculate mean/std over repetitions
  for v=1:length(FR_vars)
    var=FR_vars{v};
    pop_FR_mu=zeros(num_bins,num_repetition_sets);
    pop_FR_sd=zeros(num_bins,num_repetition_sets);
    for i=1:num_repetition_sets
      sims=repetition_set_indices(i,:);
      pop_FR_mu(:,i)=mean(stats.(var).pop_FR(:,sims),2);
      pop_FR_sd(:,i)=std(stats.(var).pop_FR(:,sims),[],2);      
    end
    stats.(var).repetition_sets.pop_FR_mu=pop_FR_mu;
    stats.(var).repetition_sets.pop_FR_sd=pop_FR_sd;
    stats.(var).repetition_sets.parameters=params_repetitions;
    stats.(var).repetition_sets.varied=varied_repetitions;           
  end
end

% 2.3 calc stats over parameter sweeps
if ~isempty(options.sweep_parameter)
  % get varied params without sweep parameter
  varied_sweeps=varied(~sweep_varied_idx);
  params_sweeps=unique(params_all(~sweep_varied_idx,:)','rows','stable')';
    % each row of params_sweeps contains the varied values besides input
    % frequency (each row has num_freqs corresponding sims).
    % [params_sweeps] = (num_sims/num_freqs) x (num_varied-1)
  num_sweeps=size(params_sweeps,2); % = (num_sims/num_freqs)
  % get indices to sims for each sweep
  [~,LOCB]=ismember(params_all(~sweep_varied_idx,:)',params_sweeps','rows');
    % LOCB is an array indicating which sets of params varied match the
    % subsets in params_sweeps
  sweep_indices=nan(num_sweeps,num_sweep_values);
  for i=1:num_sweeps
    sweep_indices(i,:)=find(LOCB==i);
      % each row of sweep_indices contains indices into data for one
      % sweep across input frequencies.
  end
  if ~isempty(options.repetition_parameter)
    % collect params for sweep sets (multiple realizations of a sweep)
    other_idx_sweeps=(~ismember(varied_sweeps,options.repetition_parameter));
    varied_sweeps_2=varied_sweeps(other_idx_sweeps);
    params_sweeps_sets=unique(params_sweeps(other_idx_sweeps,:)','rows','stable')';
    num_sweeps_sets=size(params_sweeps_sets,2);
    [~,LOCB]=ismember(params_sweeps(other_idx_sweeps,:)',params_sweeps_sets','rows');
    sweep_set_indices=nan(num_sweeps_sets,num_repetition_values);
    for i=1:num_sweeps_sets
      sweep_set_indices(i,:)=find(LOCB==i);
    end
  end
  
  % 2.3.1 calc stats for each sweep
  for v=1:length(FR_vars)
    var=FR_vars{v};
    sweep_pop_FR=zeros(num_sweep_values,num_sweeps);
    sweep_pop_FR_max=zeros(1,num_sweeps);
    sweep_pop_FR_max_amp=zeros(1,num_sweeps);
    for i=1:num_sweeps
      dat=data(sweep_indices(i,:));
      sweep_pop_FR(:,i)=cellfun(@(x)mean(x(:)),{dat.(var)});
      maxval=max(sweep_pop_FR(:,i));
      maxind=find(sweep_pop_FR(:,i)==maxval,1,'first');
      if maxval==0 || length(maxind)>1
        sweep_pop_FR_max(1,i)=nan;
        sweep_pop_FR_max_amp(1,i)=nan;
      else
        sweep_pop_FR_max(1,i)=sweep_values(maxind);
        sweep_pop_FR_max_amp(1,i)=maxval;
      end
    end
    stats.(var).sweep_sets.sweep_pop_FR=sweep_pop_FR;     
    stats.(var).sweep_sets.sweep_pop_FR_max=sweep_pop_FR_max;
    stats.(var).sweep_sets.sweep_pop_FR_max_amp=sweep_pop_FR_max_amp;
    stats.(var).sweep_sets.parameters=params_sweeps; 
    stats.(var).sweep_sets.varied=varied_sweeps; 
    stats.(var).sweep_sets.sweep_parameter=options.sweep_parameter;
    stats.(var).sweep_sets.(options.sweep_parameter)=sweep_values;

    % 2.3.2 calc sweep stats over multiple realizations
    if ~isempty(options.repetition_parameter) 
      % calculate mean/std resonance frequency
      sweep_pop_FR_mu=zeros(num_sweep_values,num_sweeps_sets);
      sweep_pop_FR_sd=zeros(num_sweep_values,num_sweeps_sets);
      sweep_pop_FR_max_mu=zeros(1,num_sweeps_sets);
      sweep_pop_FR_max_sd=zeros(1,num_sweeps_sets);
      sweep_pop_FR_max_amp_mu=zeros(1,num_sweeps_sets);
      sweep_pop_FR_max_amp_sd=zeros(1,num_sweeps_sets);
      for i=1:num_sweeps_sets
        inds=sweep_set_indices(i,:);
        sweep_pop_FR_mu(:,i)=mean(sweep_pop_FR(:,inds),2);
        sweep_pop_FR_sd(:,i)=std(sweep_pop_FR(:,inds),[],2);
        sweep_pop_FR_max_mu(:,i)=mean(sweep_pop_FR_max(1,inds),2);
        sweep_pop_FR_max_sd(:,i)=std(sweep_pop_FR_max(1,inds),[],2);          
        sweep_pop_FR_max_amp_mu(:,i)=mean(sweep_pop_FR_max_amp(1,inds),2);
        sweep_pop_FR_max_amp_sd(:,i)=std(sweep_pop_FR_max_amp(1,inds),[],2);          
      end
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_mu=sweep_pop_FR_mu;               
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_sd=sweep_pop_FR_sd;               
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_mu=sweep_pop_FR_max_mu;        
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_sd=sweep_pop_FR_max_sd;        
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_amp_mu=sweep_pop_FR_max_amp_mu;        
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_amp_sd=sweep_pop_FR_max_amp_sd;        
      stats.(var).sweep_sets.repetition_sets.parameters=params_sweeps_sets;
      stats.(var).sweep_sets.repetition_sets.varied=varied_sweeps_2;    
    else
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_mu=sweep_pop_FR;               
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_sd=zeros(size(sweep_pop_FR));          
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_mu=sweep_pop_FR_max;        
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_sd=zeros(size(sweep_pop_FR_max));        
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_amp_mu=sweep_pop_FR_max_amp;        
      stats.(var).sweep_sets.repetition_sets.sweep_pop_FR_max_amp_sd=zeros(size(sweep_pop_FR_max_amp));        
      stats.(var).sweep_sets.repetition_sets.parameters=params_sweeps;
      stats.(var).sweep_sets.repetition_sets.varied=varied_sweeps;            
    end    
  end
end

%% 3.0 Power spectral analysis
data=CalcPower(data,varargin{:});

% get list of variables whose power should be averaged over realizations
MUA_vars=data(1).results(~cellfun(@isempty,regexp(data(1).results,'.*_MUA$')));
SUA_vars=data(1).results(~cellfun(@isempty,regexp(data(1).results,'.*_SUA$')));
stats.MUA_variables=MUA_vars;
stats.SUA_variables=SUA_vars;
frequency_MUA=data(1).(MUA_vars{1}).frequency;
frequency_SUA=data(1).(SUA_vars{1}).frequency;
num_freqs=length(frequency_MUA);

% stats.(MUA_var).
%     [Pxx]               = num_freqs       x num_sims
%     [frequency]         = num_freqs       x 1
%     [PeakFreq]          = 1               x num_sims
%     [PeakArea]          = 1               x num_sims
%     [parameters]        = num_varied      x num_sims
%     [varied]            = num_varied      x 1
% stats.(MUA_var).repetition_sets.
%     [Pxx_mu]            = num_freqs       x num_rep_sets(=num_sims/num_repetitions)
%     [Pxx_sd]            = num_freqs       x num_rep_sets
%     [PeakFreq_mu]       = 1               x num_rep_sets
%     [PeakFreq_sd]       = 1               x num_rep_sets
%     [parameters]        = (num_varied-1)  x num_rep_sets, excludes repetition_parameter
%     [varied]            = (num_varied-1)  x 1
% stats.(MUA_var).sweep_sets.

% stats.(SUA_var).
%     [Pxx]               = num_freqs       x num_sims, averaged over cells
%     [frequency]         = num_freqs       x 1
%     [PeakFreq]          = 1               x num_sims, averaged over cells
%     [PeakArea]          = 1               x num_sims, averaged over cells
%     [parameters]        = num_varied      x num_sims
%     [varied]            = num_varied      x 1
% stats.(SUA_var).repetition_sets.
%     (same as stats.(MUA_var).repetition_sets.)
% stats.(SUA_var).sweep_sets.
%     (same as stats.(MUA_var).sweep_sets.)

% 2.1 MUA

% 2.1.1 population analysis
% calculate population means for each Power_MUA variable
for v=1:length(MUA_vars)
  var=MUA_vars{v};
  stats.(var).frequency=frequency_MUA;
  Pxx=zeros(num_freqs,num_sims);
  PeakFreq=zeros(1,num_sims);
  PeakArea=zeros(1,num_sims);
  for i=1:num_sims
    Pxx(:,i)=mean(data(i).(var).Pxx,2); % average over cells
    PeakFreq(i)=mean(data(i).(var).PeakFreq,2); % average over cells
    PeakArea(i)=mean(data(i).(var).PeakArea,2); % average over cells
  end
  % store results for this population variable
  stats.(var).Pxx=Pxx;              % num_freqs x num_sims
  stats.(var).PeakFreq=PeakFreq;
  stats.(var).PeakArea=PeakArea;
  stats.(var).parameters=params_all;
  stats.(var).varied=varied;
  % 2.1.2 calc stats over multiple realizations
  if ~isempty(options.repetition_parameter)
    Pxx_mu=zeros(num_freqs,num_repetition_sets);
    Pxx_sd=zeros(num_freqs,num_repetition_sets);
    PeakFreq_mu=zeros(1,num_repetition_sets);
    PeakFreq_sd=zeros(1,num_repetition_sets);
    PeakArea_mu=zeros(1,num_repetition_sets);
    PeakArea_sd=zeros(1,num_repetition_sets);
    for i=1:num_repetition_sets
      sims=repetition_set_indices(i,:);
      Pxx_mu(:,i)=mean(Pxx(:,sims),2);
      Pxx_sd(:,i)=std(Pxx(:,sims),[],2);
      PeakFreq_mu(:,i)=mean(PeakFreq(sims));
      PeakFreq_sd(:,i)=std(PeakFreq(sims));      
      PeakArea_mu(:,i)=mean(PeakArea(sims));
      PeakArea_sd(:,i)=std(PeakArea(sims));      
    end
    Pxx_mu(isnan(Pxx_mu))=0; Pxx_sd(isnan(Pxx_sd))=0;
    PeakFreq_mu(isnan(PeakFreq_mu))=0; PeakFreq_sd(isnan(PeakFreq_sd))=0;
    PeakFreq_mu(isnan(PeakFreq_mu))=0; PeakArea_sd(isnan(PeakArea_sd))=0;
    stats.(var).repetition_sets.Pxx_mu=Pxx_mu;
    stats.(var).repetition_sets.Pxx_sd=Pxx_sd;
    stats.(var).repetition_sets.PeakFreq_mu=PeakFreq_mu;
    stats.(var).repetition_sets.PeakFreq_sd=PeakFreq_sd;
    stats.(var).repetition_sets.PeakArea_mu=PeakArea_mu;
    stats.(var).repetition_sets.PeakArea_sd=PeakArea_sd;
    stats.(var).repetition_sets.parameters=params_repetitions;
    stats.(var).repetition_sets.varied=varied_repetitions;           
  end
end

% 2.2 SUA

% 2.2.1 population analysis
% calculate population means for each Power_SUA variable
for v=1:length(SUA_vars)
  var=SUA_vars{v};
  stats.(var).frequency=frequency_SUA;
  Pxx=zeros(num_freqs,num_sims);
  PeakFreq=zeros(1,num_sims);
  PeakArea=zeros(1,num_sims);
  for i=1:num_sims
    Pxx(:,i)=mean(data(i).(var).Pxx,2); % average over cells
    PeakFreq(i)=mean(data(i).(var).PeakFreq,2); % average over cells
    PeakArea(i)=mean(data(i).(var).PeakArea,2); % average over cells
  end
  % store results for this population variable
  stats.(var).Pxx=Pxx;              % num_freqs x num_sims
  stats.(var).PeakFreq=PeakFreq;
  stats.(var).PeakArea=PeakArea;
  stats.(var).parameters=params_all;
  stats.(var).varied=varied;
  % 2.2.2 calc stats over multiple realizations
  if ~isempty(options.repetition_parameter)
    Pxx_mu=zeros(num_freqs,num_repetition_sets);
    Pxx_sd=zeros(num_freqs,num_repetition_sets);
    PeakFreq_mu=zeros(1,num_repetition_sets);
    PeakFreq_sd=zeros(1,num_repetition_sets);
    PeakArea_mu=zeros(1,num_repetition_sets);
    PeakArea_sd=zeros(1,num_repetition_sets);
    for i=1:num_repetition_sets
      sims=repetition_set_indices(i,:);
      Pxx_mu(:,i)=mean(Pxx(:,sims),2);
      Pxx_sd(:,i)=std(Pxx(:,sims),[],2);
      PeakFreq_mu(:,i)=mean(PeakFreq(sims));
      PeakFreq_sd(:,i)=std(PeakFreq(sims));      
      PeakArea_mu(:,i)=mean(PeakArea(sims));
      PeakArea_sd(:,i)=std(PeakArea(sims));      
    end
    Pxx_mu(isnan(Pxx_mu))=0; Pxx_sd(isnan(Pxx_sd))=0;
    PeakFreq_mu(isnan(PeakFreq_mu))=0; PeakFreq_mu(isnan(PeakFreq_mu))=0;
    PeakFreq_mu(isnan(PeakFreq_mu))=0; PeakArea_sd(isnan(PeakArea_sd))=0;
    stats.(var).repetition_sets.Pxx_mu=Pxx_mu;
    stats.(var).repetition_sets.Pxx_sd=Pxx_sd;
    stats.(var).repetition_sets.PeakFreq_mu=PeakFreq_mu;
    stats.(var).repetition_sets.PeakFreq_sd=PeakFreq_sd;
    stats.(var).repetition_sets.PeakArea_mu=PeakArea_mu;
    stats.(var).repetition_sets.PeakArea_sd=PeakArea_sd;
    stats.(var).repetition_sets.parameters=params_repetitions;
    stats.(var).repetition_sets.varied=varied_repetitions;           
  end  
end

