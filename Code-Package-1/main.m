%------- The FDTD formulation of the second order scalar Helmholtz wave  
%        equation for the continuous dipole source excitation.
%------- Code performs the following :
% 1) FDTD wave scattering of a random slab medium for a continuous scalar 
%    dipole source.
% 2) Establishes a Total-field/Scattered field (TF/SF) boundary on the source
%    side for estimating the reflectance.
% 3) Phase Sensitive Lock-in Detection (PSD) at the transmitting and 
%    reflecting boundaries for phase and amplitude measurement of outgoing
%    steady state waves.
% 4) Evaluates wave current J(r,t) and estimates the total transmission(T) 
%    and reflection(R). 
% 5) Performs phase conjugation and estimate the field profile at the focus 
%    together with the estimation of the total transmission and reflection.

clear 
close all

GPU_flag=0;            % If the GPU is used, set the flag to 1 for 
                       %      faster execution in case of a large disorder. 
                       % If no GPU is available, set the flag=0.
visualisation_flag=0;  % 1 for switching on the field visualisation 
                       %    (slower code execution)
                       % 0 for switching off the field visualisation
                       %    (faster code execution) 
data_saving_flag=1;    % 1 for saving the data and 0 for not saving.
%------------------- initialisation : fundamental parameters --------------
NT_steady_state_forward_propagation=35000; % No of time steps for forward 
                                         %   propagation
NT_steady_state_PC=35000; 
                     % No of time steps for phase conjugation
                     % These are the no of time steps to attain the steady 
                     % state. Please increase the value according to the 
                     % size of the slab to ensure steady state.
                     
% For quick computation of results (to reach steady-state quickly),
% it is suggested to use the slab geometry where krefW >> krefL, which 
% reduces the transverse oscillation of the source along the y axis, due to
% the mirror reflectance boundary condition at the transverse boundaries.
krefL=100;              % Dimensionless longitudinal thickness of the 
                        %           scattering slab along z, ie kref*L
aspect_ratio=5;  % Aspect ratio of the slab size = W/L; 
                 % Recommended to take the aspect ratio value larger than 3
z_domain_size_factor=3; % kref*Z=z_domain_size_factor*krefL
                          % Number of times krefL, is the domain z size.
                          
krefW= ...
    aspect_ratio*krefL;  % Initial value for the dimensionless transverse 
                         %         length of the slab along y, ie kref*W    
nref=1.9;                % Background(reference medium) refractive index 
sigma=0.98;              % A parameter, usually < 1 used to control the  
                         %       degree of random scattering.
no_of_PSD_cycles=10;     % Number of Phase sensitive detection (PSD)
                         %  wave cycles used for averaging the amplitde
                         %  and phase values upon achieving steady-state. 
                         
% Call initialisation script for initialising rest of the parameters and 
%    return an initialisation data structure.
init_data=initialisation(GPU_flag,NT_steady_state_forward_propagation,krefL,...
           krefW,nref,sigma,z_domain_size_factor,no_of_PSD_cycles);
%------------- Define the source  -----------------------------------------
[source_data] = define_dipole_point_source(GPU_flag,init_data); 
                       % Returns the source data structure.       
source_data.clip_lim=0.2; % Colorbar clip limit for 
                          % better contrast to begin with                      
source_data.t_sign=+1; 
          % +1 for forward propagating source and 
          % -1 for time reversed version while performing phase conjugation
src_z_start=floor(init_data.start_slab1/2); 
   % Source position z value. It is set as half way between slab and the
   %                          left boundary.
source_data.src_z=[src_z_start src_z_start+1];     % Dipole source location
source_data.TFSF_Nz=floor(source_data.src_z(1)/2); 
 % Total-field Scattered field (TF/SF) boundary location 
 % Set as half way between the source and the left boundary 
source_data.TFSF_region=1:source_data.TFSF_Nz; % Defining the TF/SF region
%-- Forward FDTD wave propagation : Left to right -------------------------
tic
[Efield_array_collection] = FDTD_main_loop(GPU_flag,visualisation_flag,...
    init_data,source_data); % Main FDTD time-stepping loop
%---- Phase sensitive detection (PSD) for the forward propagation ---------
init_data.detec_trans_position=...
 floor(((init_data.start_slab1 + init_data.scat_width) + init_data.Nz)/2);
  % Transmission measurement z index location
  % Half way between the end of the slab and the right boundary                              
init_data.detec_refl_position=floor(source_data.TFSF_Nz/2); 
  % Reflection measurement z index location
  % Half way between the TF/SF boundary and the left boundary                                                            
init_data.detec_src_position=init_data.start_slab1;
          % Source transmission measurement z index location
          % Placed at the slab surface incidence

% Next is the PSD of the amplitude and phase of the scattered and 
%    the source wave, followed by the estimation of the total transmission
%    and reflection
PSD_array_collection= ...
  phase_sensitive_detection(GPU_flag,visualisation_flag,Efield_array_collection,...
                            init_data,source_data);
                   % Returns a data structure containing amplitude
                   % and phase of transmitted, reflected and source
                   % waves along with total transmission and reflection.
sprintf('Total time taken for forward propagation = %f min',toc/60)
%--------------------------------------------------------------------------
figure(2)
title('$Forward~propagation$','Interpreter','Latex')
drawnow
%save data for the paper
if data_saving_flag==1
save('PSD_array_collection_forward.mat','PSD_array_collection')
save('Efield_array_collection_forward.mat','Efield_array_collection')
save('init_data_forward.mat','init_data')
%------------- save graphic images --------------------------------------
%save image 
% figure(2)
% set(gca,'FontSize',12);
% exportgraphics(gcf,'forwardFDTD.eps','ContentType','vector')
end


tic
%------ phase conjugation : Right to Left ---------------------------------
init_data.NT_steady_state=NT_steady_state_PC;
source_data.t_sign=-1; % Setting the sign flag for phase conjugation.
                       % -1 implies time-reversal
source_data.src_z=[init_data.detec_trans_position init_data.detec_trans_position+1];  
                  % Defining the dipole source location
source_data.src_y=2:init_data.Ny+1; % Setting the active y region of the 
                                    % source
source_data.TFSF_Nz=floor((source_data.src_z(2)+init_data.Nz)/2);  
             % Setting the TFSF boundary location for the new source, 
             % this time for phase conjugation.
source_data.TFSF_region=source_data.TFSF_Nz:init_data.Nz; % Setting the 
                % TFSF region while performing phase conjugation
source_data.wave_amp= ...
    repmat(PSD_array_collection.E_trans_sig(2:end-1),1,2);
   % Setting the dipole source amplitude for phase conjugation, which is
   % based on the transmitted wave itself
source_data.wave_amp= ...
    0.1.*source_data.wave_amp./max(source_data.wave_amp); % Normalising for
                                 % better contrast.
% Next, set the dipole source phase for phase conjugation, which is
% based on the transmitted wave itself                                 
source_data.wave_phase= ...
    repmat(PSD_array_collection.theta_trans_sig(2:end-1),1,2);
source_data.wave_phase(:,2)=source_data.wave_phase(:,2)+pi; 
            % dipole phase factor pi is added to the second component of
            % the dipole source for actually making it a dipole source. 
source_data.clip_lim=gather(max(source_data.wave_amp(:)))+0.8;
          %  clip_lim is the limit value for the colour map clip to 
          % adjust the contrast of the field images. 

[Efield_array_collection] = ...
    FDTD_main_loop(GPU_flag,visualisation_flag,init_data,source_data);
                                       % Main FDTD time-stepping loop
%------------ Phase sensitive detection -----------------------------------
init_data.detec_trans_position=src_z_start;
  % Setting the transmission detection location
  % same as that of the original point source location
init_data.detec_refl_position=floor((source_data.TFSF_Nz+init_data.Nz)/2);
  % Reflection detection location
  % Set as half way between the TF/SF boundary and the right boundary
init_data.detec_src_position= ...
  init_data.start_slab1+init_data.scat_width+1; 
   % Source transmission detection location
   % Set as the RHS incident surface of the slab
   
% Next, obtain the PSD for the phase conjugation process and estimate 
% the associated total transmission and reflection
PSD_array_collection=phase_sensitive_detection(GPU_flag,visualisation_flag,...
    Efield_array_collection,init_data,source_data);
sprintf('Total time taken for phase conjugation = %f min',toc/60)

figure(6)
title('$Phase~conjugation$','Interpreter','Latex')
drawnow 
%save data for the thesis/paper
if data_saving_flag==1
save('PSD_array_collection_PC.mat','PSD_array_collection')
save('Efield_array_collection_PC.mat','Efield_array_collection')
% figure(6)
% xlabel('$k_{ref}z$','Interpreter','Latex')
% ylabel('$k_{ref}y$','Interpreter','Latex')
% set(gca,'FontSize',12);
% exportgraphics(gcf,'PCFDTD.eps','ContentType','vector')
end 

%---------------- plot the focus ------------------------------------------
foc_reg= ...
find(PSD_array_collection.E_trans_sig==max(PSD_array_collection.E_trans_sig));
       % finding the peak y location
E_trans_normalised= ...
(PSD_array_collection.E_trans_sig)./max(PSD_array_collection.E_trans_sig(:));  
                     % Normalise the transmission field data with
                     %                    its peak value for the dB plot.                                                     
figure('Position', [0 0 1500 800])
subplot(1,3,1) 
plot(E_trans_normalised.^2,(0:init_data.Ny+1).*init_data.kref*init_data.dy,...
    '-*','MarkerSize',4); % Plotting normalised intensity
ylabel('$k_{ref}y$','Interpreter','Latex')
set(gca,'xdir','reverse')
title('$a)Normalized~E^2~at~focus$','Interpreter','Latex')
axis tight;
set(gca,'FontSize',24)                                                      
                                                      
subplot(1,3,2)
plot(mag2db(E_trans_normalised),(0:init_data.Ny+1).*init_data.kref*init_data.dy,...
    '-*','MarkerSize',4);
ylabel('$k_{ref}y$','Interpreter','Latex')
set(gca,'xdir','reverse')
xlabel('$dB$','Interpreter','Latex')
title('$b)In~dB$','Interpreter','Latex')
axis tight;
set(gca,'FontSize',24)

subplot(1,3,3)
Ny_uplim  =(foc_reg + floor(foc_reg/2));   % For cropping the focus for
Ny_lowlim =(foc_reg - floor(foc_reg/2));   %    zooming in 
foc_reg_lamda_ratio=...
((((Ny_lowlim-foc_reg:Ny_uplim-foc_reg)-1).*init_data.dy)./init_data.lambda_ref);
plot(mag2db(E_trans_normalised(Ny_lowlim:Ny_uplim)),foc_reg_lamda_ratio,...
    '-*','MarkerSize',4);
set(gca,'xdir','reverse')
ylabel('$y/\lambda_{ref}$','Interpreter','Latex')
xlabel('$dB$','Interpreter','Latex')
title('$c)In~dB$','Interpreter','Latex')
ylim([-2 2])
set(gca,'FontSize',24)
%--------------------------------------------------------------------------
