%------- The FDTD formulation for demonstrating Raman-Nath STM guidestar
%        embedded between two disordered slabs
%------- Code performs the following :
% 1) FDTD wave interaction from a time dependent dielectric STM grating 
%    operating in the Raman-Nath regime placed in between two disordered 
%    slabs. 
% 2) Performs phase sensitive detection (PSD) of the modulated wavefronts
%    emerging from the Raman-Nath grating transmitted through the
%    disordered slab.
% 3) Performs IOPC as described in the main paper from both sides of the
%    disorder, in transmission
%---------------------- Need to ensure that T+R=1 -------------------------
clear 
close all

estimate_disorder_transmission_flag=1; % Set the flag 1 if the transport 
                       %               regime of the disorder needs to be
                       %               profiled. 

GPU_flag=0;            % If the GPU is used, set the flag to 1 for 
                       %      faster execution in case of a large disorder. 
                       % If no GPU is available, set the flag=0.
visualisation_flag=0;  % 1 for switching on the field visualisation 
                       %    (slower code execution)
                       % 0 for switching off the field visualisation
                       %    (faster code execution) 
no_of_roundtrips=10;                         
data_saving_flag=1;    % 1 for saving the data and 0 for not saving.
%------------------- initialisation : fundamental parameters --------------
NT_steady_state_forward_propagation=40000; % No of time steps for forward 
                                         %   propagation
NT_steady_state_PC=40000; % No of time steps for phase conjugation
                       % No of time steps to attain the steady state and
                       %    please adjust the value according to the 
                       %    size of the slab to ensure steady state.
% For quick computation of results (to reach steady-state quickly),
% it is suggested to use the slab geometry where krefW >> krefL, which 
% reduces the transverse oscillation of the source along the y axis, due to
% the mirror reflectance boundary condition at the transverse boundaries.
krefL=50;              % Dimensionless longitudinal thickness of the 
                        %           scattering slab along z, ie kref*L
aspect_ratio=10;         % Aspect ratio of the slab size = W/L;                           
z_domain_size_factor=2; % kref*Z=z_domain_size_factor*krefL
                        % Number of times krefL, is the domain z size.
                        % Take larger than or equal to 3.
krefW= ...
    aspect_ratio*krefL;  % Initial value for the dimensionless transverse 
                         %         length of the slab along y, ie kref*W    
nref=1.8;                % Background(reference medium) refractive index 
sigma=0.95;               % A parameter, usually < 1 used to control the  
                         %       degree of random scattering.
no_of_PSD_cycles_forward=10; % Number of Phase sensitive detection (PSD)
                         %  wave cycles used for averaging the amplitde
                         %  and phase values upon achieving steady-state. 
no_of_PSD_cycles_STM=1;  % For estimating the modulated wavefronts
                         % Setting 1 for fast estimation
% Call initialisation script for initialising rest of the parameters and 
%    return an initialisation data structure.
init_data=initialisation(GPU_flag,NT_steady_state_forward_propagation,krefL,...
           krefW,nref,sigma,z_domain_size_factor,no_of_PSD_cycles_forward);
                         
%------------- Define the source  -----------------------------------------
[source_data] = define_dipole_beam_source(GPU_flag,init_data); 
source_data.clip_lim=2.5; % Colorbar clip limit for better contrast

source_data.t_sign=+1; % +1 for forward propagating source and 
                       % -1 for time reversed version
src_z_start=floor(init_data.start_slab1/2);  % Half way between slab and the
                                            % left boundary
source_data.src_z=[src_z_start src_z_start+1];     % Dipole source location
source_data.TFSF_Nz=floor(source_data.src_z(1)/2); % TFSF boundary location 
                        % Half way between the source and the left boundary 
source_data.TFSF_region=1:source_data.TFSF_Nz;

if estimate_disorder_transmission_flag==1
estimate_disorder_transmission % Estimate the disorder transmission 
                               % conditionally
figure(2)                               
title('$Forward~propagation~without~STM$','Interpreter','Latex')                               
end
                      
%--- Section 2 : Add an STM region and estimate the modulated wavefronts ---
%--------------------------- STM initialisation ---------------------------
del_eps=0.5; % dielectric constant modulation by a space-time modulated curve  
spat_freq_fact=0.4; % spatial frequency of modulation
f_mod=0.2*(init_data.fsrc);  % temporal frequency of the modulation
grat_width=ceil(3.0*pi/(init_data.kref*init_data.dz)); % space-time modulated grating width
[STM_data]=define_STM_parameters(GPU_flag,init_data,source_data,del_eps,...
    spat_freq_fact,f_mod,grat_width);
init_data.no_of_PSD_cycles=no_of_PSD_cycles_STM;
forward_propagation_with_STM_Left_to_Right; 
title('$Forward~propagation~with~STM$','Interpreter','Latex')                               
% The script used to perform forward propagation from left to right 
% in presence of the STM modulation and detect the n=+1 sideband (upper).

%-------------------------- Perform IOPC ----------------------------------

for r_count=1:no_of_roundtrips
%---- Section 3 : Phase conjugation of the upper sideband from the right --
nchosen_R2L=+1; % Choose the upconverted wavefront 
                % for the phase conjugation from the right
phase_conjugation_STM_Right_to_left 
% The script used to perform phase conjugation from right to left, in presence
% of the STM modulation and detect the n=+1 sideband (upper).

%----- Section 4 : Phase conjugation of the omega_0 wavefront from the left
nchosen_L2R=0; % Choose the center frequency omega_0 corresponding to 
           % n=0, to be phase conjugated from the left. 
           % Such a wavefront at omega_0, actually emerged on the left,
           % due the STM modulation of the phase conjugated upperside band.
           % The STM region, downconverted the phase conjugated upperside band from the
           % right, into the center frequency, while passing through the STM region.
phase_conjugation_STM_Left_to_Right;
end

%--------------------------- Analysis -------------------------------------
figure
for r_count=1:no_of_roundtrips 
filename=sprintf('Efield_array_collection_PC_L2R_%d.mat',r_count);
load(filename)
field_mag_at_STM=abs(Efield_array_collection.EnTFSF(:,STM_data.STM_start_z));
plot(mag2db(field_mag_at_STM./max(field_mag_at_STM)),...
    (0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',1)
set(gca,'xdir','reverse')
xlabel('$dB$','Interpreter','Latex')
title('$Normalised~Intensity~at~the~STM~plane$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
set(gca,'FontSize',22)
drawnow;
hold on
end

%------------------------- Focus profile on focussing ---------------------

chosen_round_trip_nos=[10];
r_count=1;
STM_y_index=(STM_data.STM_start_y:STM_data.STM_end_y);

figure('position',[0 0 1200 600])
subplot(1,2,1)
filename=sprintf('Efield_array_collection_PC_R2L_%d.mat',...
    chosen_round_trip_nos(r_count));
load(filename)
field_mag_at_STM=abs(Efield_array_collection.EnTFSF(:,STM_data.STM_start_z));
norm_dB_intensity_at_STM=mag2db(field_mag_at_STM./max(field_mag_at_STM));
plot(norm_dB_intensity_at_STM,...
    (0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',1)
set(gca,'xdir','reverse')
xlabel('$dB$','Interpreter','Latex')
title('$Normalised~Intensity~at~the~STM~plane$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
set(gca,'FontSize',22)
xlim([-20 0])
hold on
line([0 -40],...
    [(STM_data.STM_start_y).*init_data.kref*init_data.dy ...
    (STM_data.STM_start_y).*init_data.kref*init_data.dy], ...
    'Color','magenta','LineStyle','--','LineWidth',1)
hold on
line([0 -40],...
    [(STM_data.STM_end_y).*init_data.kref*init_data.dy ...
    (STM_data.STM_end_y).*init_data.kref*init_data.dy], ...
    'Color','magenta','LineStyle','--','LineWidth',1)
 db_background=-6.021; %
 
 line([db_background db_background],...
[(0).*init_data.kref*init_data.dy (length(field_mag_at_STM)-1).*init_data.kref*init_data.dy], ...
 'Color','magenta','LineStyle','--','LineWidth',1)
%--------------------------------------------------------------------------
subplot(1,2,2)
eta_STM_to_plot=sqrt(STM_data.eps_ij_n(1));
plot(eta_STM_to_plot(:,STM_data.STM_start_z),...
    (0:init_data.Ny-1).*init_data.kref*init_data.dy,'-*','MarkerSize',1)
axis tight
set(gca,'xdir','reverse')
set(gca,'FontSize',22)
xlabel('$\eta_y~at~the~STM~plane$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
line([max(eta_STM_to_plot(:,STM_data.STM_start_z)) ...
    min(eta_STM_to_plot(:,STM_data.STM_start_z))],...
    [(STM_data.STM_start_y).*init_data.kref*init_data.dy ...
    (STM_data.STM_start_y).*init_data.kref*init_data.dy], ...
    'Color','magenta','LineStyle','--','LineWidth',1)
hold on
line([max(eta_STM_to_plot(:,STM_data.STM_start_z)) ...
    min(eta_STM_to_plot(:,STM_data.STM_start_z))],...
    [(STM_data.STM_end_y).*init_data.kref*init_data.dy ...
    (STM_data.STM_end_y).*init_data.kref*init_data.dy], ...
    'Color','magenta','LineStyle','--','LineWidth',1)
drawnow;
