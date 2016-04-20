% inputs

% square wave
ac=2;       % positive peak amplitude (peak-to-peak amp = 2*ac)
duty=25;    % percent of cycle in which signal is positive
f=5;        % frequency, Hz
t=0:1e-4:1; % time, sec
y = ac*square(2*pi*f*t,duty); % range -ac to ac
figure; plot(t,y); ylim([-1.5 1.5]*ac);

% sin wave
ac=2;       % positive peak amplitude (peak-to-peak amp = 2*ac)
f=5;        % frequency, Hz
t=0:1e-4:1; % time, sec
y = ac*sin(2*pi*f*t); % range -ac to ac
figure; plot(t,y); ylim([-1.5 1.5]*ac);

% homogeneous poisson process
dt=.01/1000;  % sec
tdur=30;      % sec
t=0:dt:tdur;  % sec, time vector
DC=10; % Hz
P=poissrnd(DC*ones(size(t))*dt); % poisson process
sum(P)/tdur % spikes/sec

%% nonhomogeneous poisson with sinusoidal rate modulation
dt=.01/1000;  % sec
f=5;          % modulation frequency, Hz
t=0:dt:1;     % time, sec
y = sin(2*pi*f*t); % -1 to 1
DC=5000;  % Hz
AC=DC;    % Hz
lambda=max(DC+AC*y,0); % Hz
P=poissrnd(lambda*dt); % poisson process
% examine instantaneous rate of simulated input
kwidth=.02;              % s, width of gaussian for kernel regression
Ts=.01;                   % s, set to this effective time step for rate process before regression
s=find(P); % when spikes occur
raster=ones(length(s),2); % raster(:,1) -> spike times, raster(:,2) -> cell index
raster(:,1)=t(s);
[iFR,itime]=NWgaussKernelRegr(t,raster,1,kwidth,Ts);
figure;
subplot(2,1,1); plot(t,lambda/1000); xlabel('time'); ylabel('lambda [kHz]'); line(xlim,[DC DC]/1000);
subplot(2,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]'); line(xlim,[DC DC]/1000);
% ylim([0 DC+AC]);

%% nonhomogeneous poisson population with sinusoidal rate modulation
Ninputs=100;
dt=.01/1000;  % sec
tdur=1;       % sec
f=5;          % Hz, modulation frequency
t=0:dt:tdur;  % sec, time vector
y=sin(2*pi*f*t); % -1 to 1
DC=10000; % Hz, steady lambda component
AC=DC;   % Hz, modulation lambda component
lambda=max(DC+AC*y,0); % Hz, Poisson rate modulation
lambda=repmat(lambda/Ninputs,[Ninputs 1]);
P=poissrnd(lambda*dt); % poisson process
sum(P(:))/tdur % spikes/sec

% examine instantaneous rate of simulated input
kwidth=.02;              % s, width of gaussian for kernel regression
Ts=.01;                   % s, set to this effective time step for rate process before regression
% collect spikes in raster
raster=[];
for i=1:Ninputs
  spks=find(P(i,:));
  raster=cat(1,raster,[t(spks)' i*ones(length(spks),1)]);
end
[iFR,itime]=NWgaussKernelRegr(t,raster,1:Ninputs,kwidth,Ts);
figure('position',[200 80 1350 900]);
subplot(3,1,1); plot(t,lambda); xlabel('time'); ylabel('lambda [kHz]');
subplot(3,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]');
subplot(3,1,3); imagesc(t,1:Ninputs,P); colormap(1-gray); axis xy

%% sinusoidal poisson with controlled spike bursts
% allowing for AC>DC, set DC and AC s.t. Nspikes occur over Tburst/2
dt=.01/1000;  % sec
tdur=1;       % sec
f=5;          % Hz, modulation frequency
Tburst=1/f;   % sec, width of burst
Nspikes=250;  % # spikes per burst
DC=2000;
AC=(Nspikes-DC*(Tburst/2))*(pi/(Tburst*dt))/100000;
t=0:dt:tdur;            % sec, time vector
y=sin(2*pi*f*t);        % -1 to 1
lambda=max(DC+AC*y,0);  % kHz, nonhomogeneous poisson rate, 0 to (DC+AC)
P=poissrnd(lambda*dt);  % poisson process
[sum(P(1:nearest(t,Tburst/2))) round(sum(P)/(tdur*f))]

kwidth=.02;              % s, width of gaussian for kernel regression
Ts=.01;                   % s, set to this effective time step for rate process before regression
s=find(P); % when spikes occur
raster=ones(length(s),2); % raster(:,1) -> spike times, raster(:,2) -> cell index
raster(:,1)=t(s);
[iFR,itime]=NWgaussKernelRegr(t,raster,1,kwidth,Ts);
figure;
subplot(2,1,1); plot(t,lambda/1000); xlabel('time'); ylabel('lambda [kHz]'); line(xlim,[DC DC]/1000);
subplot(2,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]'); line(xlim,[DC DC]/1000);

% POPULATION
Ninputs=100;
lambda=max(DC+AC*y,0); % Hz, Poisson rate modulation
lambda=repmat(lambda/Ninputs,[Ninputs 1]);
P=poissrnd(lambda*dt); % poisson process
[sum(sum(P(:,1:nearest(t,Tburst/2)))) round(sum(P(:))/(tdur*f))]
% collect spikes in raster
raster=[];
for i=1:Ninputs
  spks=find(P(i,:));
  raster=cat(1,raster,[t(spks)' i*ones(length(spks),1)]);
end
[iFR,itime]=NWgaussKernelRegr(t,raster,1:Ninputs,kwidth,Ts);
figure('position',[200 80 1350 900]);
subplot(3,1,1); plot(t,lambda); xlabel('time'); ylabel('lambda [kHz]');
subplot(3,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]');
subplot(3,1,3); imagesc(t,1:Ninputs,P); colormap(1-gray); axis xy

%% sinusoidal poisson spike bursts with fixed interburst-interval

% todo: Note, adjusting AC wrt DC may not be the most realistic approach. A
% burst would represent the AC component only. As currently implemented,
% the burst is taking spikes from the background. might be better to
% constrain Nspikes to the AC component only.

dt=.01/1000;  % sec
tdur=2;       % sec
f=3;          % Hz, modulation frequency
Tburst=.2;    % sec, width of burst
% one burst cycle
yburst=sin(2*pi*(0:dt:Tburst)/Tburst);        % -1 to 1
t=0:dt:tdur;            % sec, time vector
y=zeros(size(t));
% burst times
Tibi=(1/f)-Tburst; % sec, time b/w end of one burst and beginning of the next
Tosc=Tburst+Tibi; % time b/w start of one burst and start of the next (i.e., period of effective modulation frequency)
tstart=0:Tosc:tdur; % start times for each burst
% insert bursts separated by Tibi zero's
for i=1:length(tstart)
  ton=tstart(i);
  ind=nearest(t,ton):nearest(t,ton+Tburst);
  y(ind)=yburst(1:length(ind));
end
Nspikes=2000;  % # spikes per burst
DC=0;
AC=(Nspikes-DC*(Tburst/2))*(pi/(Tburst*dt))/100000;
lambda=max(DC+AC*y,0);  % kHz, nonhomogeneous poisson rate, 0 to (DC+AC)
P=poissrnd(lambda*dt);  % poisson process
[sum(P(1:nearest(t,Tburst/2))) round(sum(P)/(tdur*f))]

kwidth=.01;              % s, width of gaussian for kernel regression
Ts=.005;                   % s, set to this effective time step for rate process before regression
s=find(P); % when spikes occur
raster=ones(length(s),2); % raster(:,1) -> spike times, raster(:,2) -> cell index
raster(:,1)=t(s);
[iFR,itime]=NWgaussKernelRegr(t,raster,1,kwidth,Ts);
figure;
subplot(2,1,1); plot(t,lambda/1000); xlabel('time'); ylabel('lambda [kHz]'); line(xlim,[DC DC]/1000);
subplot(2,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]'); line(xlim,[DC DC]/1000);

% POPULATION
Ninputs=100;
lambda=max(DC+AC*y,0); % Hz, Poisson rate modulation
lambda=repmat(lambda/Ninputs,[Ninputs 1]);
P=poissrnd(lambda*dt); % poisson process
[sum(sum(P(:,1:nearest(t,Tburst/2)))) round(sum(P(:))/(tdur*f))]
% collect spikes in raster
raster=[];
for i=1:Ninputs
  spks=find(P(i,:));
  raster=cat(1,raster,[t(spks)' i*ones(length(spks),1)]);
end
[iFR,itime]=NWgaussKernelRegr(t,raster,1:Ninputs,kwidth,Ts);
figure('position',[200 80 1350 900]);
subplot(3,1,1); plot(t,lambda); xlabel('time'); ylabel('lambda [kHz]');
subplot(3,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]');
subplot(3,1,3); imagesc(t,1:Ninputs,P); colormap(1-gray); axis xy

% note: oscillatory freq = 1/(Tburst+Tibi)

%% sinusoidal poisson spike bursts with exponentially-distributed interburst-intervals

% note: interburst interval = the time between the beginning of
% successive bursts (not the end of one and the beginning of the next).
% eliminate all IBIs less than the max Tburst used across simulations; this
% ensures that a comparable number of bursts will occur for different
% Tburst levels.

dt=.01/1000;  % sec
tdur=2;       % sec
Tburst=.01;    % sec, width of burst
maxTburst=.1;
% one burst cycle
yburst=sin(2*pi*(0:dt:Tburst)/Tburst);        % -1 to 1
t=0:dt:tdur;            % sec, time vector
y=zeros(size(t));
% burst times
meanIBI=.1; % sec, mean interburst interval
Tibi=exprnd(meanIBI,size(t));
Tibi=Tibi(Tibi>maxTburst);
tstart=[0 cumsum(Tibi)];
tstart=tstart(tstart<tdur);
% insert bursts separated by Tibi zero's
for i=1:length(tstart)
  ton=tstart(i);
  ind=nearest(t,ton):nearest(t,ton+Tburst);
  y(ind)=yburst(1:length(ind));
end
Nspikes=2000;  % # spikes per burst
DC=0;
AC=(Nspikes-DC*(Tburst/2))*(pi/(Tburst*dt))/100000;
lambda=max(DC+AC*y,0);  % kHz, nonhomogeneous poisson rate, 0 to (DC+AC)
P=poissrnd(lambda*dt);  % poisson process

kwidth=.01;              % s, width of gaussian for kernel regression
Ts=.005;                   % s, set to this effective time step for rate process before regression
s=find(P); % when spikes occur
raster=ones(length(s),2); % raster(:,1) -> spike times, raster(:,2) -> cell index
raster(:,1)=t(s);
[iFR,itime]=NWgaussKernelRegr(t,raster,1,kwidth,Ts);
figure;
subplot(2,1,1); plot(t,lambda/1000); xlabel('time'); ylabel('lambda [kHz]'); line(xlim,[DC DC]/1000);
subplot(2,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]'); line(xlim,[DC DC]/1000);

% POPULATION
Ninputs=100;
lambda=max(DC+AC*y,0); % Hz, Poisson rate modulation
lambda=repmat(lambda/Ninputs,[Ninputs 1]);
P=poissrnd(lambda*dt); % poisson process

% collect spikes in raster
raster=[];
for i=1:Ninputs
  spks=find(P(i,:));
  raster=cat(1,raster,[t(spks)' i*ones(length(spks),1)]);
end
[iFR,itime]=NWgaussKernelRegr(t,raster,1:Ninputs,kwidth,Ts);
figure('position',[200 80 1350 900]);
subplot(3,1,1); plot(t,lambda); xlabel('time'); ylabel('lambda [kHz]');
subplot(3,1,2); plot(itime,iFR/1000); xlabel('time'); ylabel('inst. rate of input [kHz]');
subplot(3,1,3); imagesc(t,1:Ninputs,P); colormap(1-gray); axis xy

%% two nonhomo poisson inputs w/ different degrees of synchrony and ~equal # of spikes:
% note: integral(a*sin(2*pi*f*t)dt,[0 1/(2f)]) = a/(pi*f)
% therefore, to match # of spikes b/w inputs w/ f1 and f2, and the first
% having a lambda AC component with amplitude AC=a1, then set input two
% AC=a2=a1*f2/f1.



S_ini=zeros(N,1); tau=2; kick=1;
S=nonhomPoissonGeneratorSpikeTimes(S_ini,lambda,tau,kick,Ninputs,tdur,dt);


