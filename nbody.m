close all;
clear all;
clc;

%% Create Your Own N-body Simulation (With Matlab/Octave)
% Philip Mocz (2020) Princeton Univeristy, @PMocz

% Simulate orbits of stars interacting due to gravity
% Code calculates pairwise forces according to Newton's Law of Gravity


%% Simulation parameters
N         = 100;   % Number of particles
t         = 0;     % current time of the simulation
tEnd      = 10;    % time at which simulation ends
dt        = 0.01;  % timestep
softening = 0.1;   % softening length
G         = 1;     % Newton's Gravitational Constant
plotRealTime = 1;  % switch on (1) for plotting as the simulation goes along


%% Generate Initial Conditions
rng(42);                % set the random number generator seed

mass = 20*ones(N,1)/N;  % total mass of particles is 20
pos = randn(N,3);       % randomly selected positions and velocities
vel = randn(N,3);

% Convert to Center-of-Mass Frame
vel = vel - mean((mass*[1 1 1]) .* vel) / mean(mass);

% calculate initial gravitational accelerations
acc = getAcc( pos, mass, G, softening );

% calculate initial energy of system
[ KE, PE ] = getEnergy( pos, vel, mass, G );

% number of timesteps
Nt = ceil(tEnd/dt);


%% save energies, particle orbits for plotting trails
pos_save = zeros(N,3,Nt+1);
pos_save(:,:,1) = pos;
KE_save = zeros(Nt+1,1);
KE_save(1) = KE;
PE_save = zeros(Nt+1,1);
PE_save(1) = PE;
t_all = (0:Nt)*dt;


%% Simulation Main Loop
fh = figure('position',[0 0 600 800]);

for i = 1:Nt
    
    % (1/2) kick
    vel = vel + acc * dt/2;
    
    % drift
    pos = pos + vel * dt;
    
    % update accelerations
    acc = getAcc( pos, mass, G, softening );
    
    % (1/2) kick
    vel = vel + acc * dt/2;
    
    % update time
    t = t + dt;
    
    % get energy of system
    [ KE, PE ] = getEnergy( pos, vel, mass, G );
    
    % save energies, positions for plotting trail
    pos_save(:,:,i+1) = pos;
    KE_save(i+1) = KE;
    PE_save(i+1) = PE;
    
    % plot in real time
    if (plotRealTime) || (i==Nt)
        subplot(3,1,1:2)
        xx = pos_save(:,1,max(i-50,1):i);
        yy = pos_save(:,2,max(i-50,1):i);
        plot(xx(:),yy(:),'.','color',[.7 .7 1]);
        hold on
        plot(pos(:,1),pos(:,2),'b.','markersize',14);
        hold off
        axis square
        axis([-2 2 -2 2])
        
        subplot(3,1,3)
        plot(t_all,KE_save,'r.')
        hold on
        plot(t_all,PE_save,'b.')
        plot(t_all,KE_save+PE_save,'k.')
        hold off
        axis([0 tEnd -300 300])
        
        drawnow
    end
    
end


%% add labels/legend
subplot(3,1,3)
xlabel('time')
ylabel('energy')
lh = legend('KE','PE','Etot');
set(lh,'location','northeast');

%% Save figure
saveas(fh,'nbody.png')

