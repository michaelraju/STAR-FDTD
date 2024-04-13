%------- The FDTD formulation for demonstrating Raman-Nath STM operation
%        without any disorder present  
%------- Code performs the following :
% 1) FDTD wave diffraction from a time dependent dielectric STM grating 
%    operating in the Raman-Nath regime  
% 2) Performs phase sensitive detection (PSD) of the modulated wavefronts
%    emerging from the Raman-Nath grating. 
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
NT_steady_state=12000;% No of time steps to attain the steady state and
                       %    please adjust the value according to the 
                       %    size of the slab to ensure steady state.

nref=1.0;                % Background(reference medium) refractive index 
sigma=0.0;               % A parameter, usually < 1 used to control the  
                         %       degree of random scattering.                         
% sigma is set to zero as there is no disorder slab in this case. 

krefL=90;              % Dimensionless longitudinal thickness of the 
                        %           scattering slab along z, ie kref*L
% krefL is not applicable in this case as there is no slab disorder
% Nevertheless, gave a value to reuse the previously developed codes 

aspect_ratio=12;         % Aspect ratio of the slab size = W/L;   
z_domain_size_factor=6; % kref*Z=z_domain_size_factor*krefL
                        % Number of times krefL, is the domain z size.
                        % Take larger than or equal to 3.
krefW= ...
    aspect_ratio*krefL;  % Initial value for the dimensionless transverse 
                         %         length of the slab along y, ie kref*W    
  

no_of_PSD_cycles=1;     % Number of Phase sensitive detection (PSD)
                         %  beat-cycles used for estimating the amplitde
                         %  and phase values upon achieving steady-state. 
                         
% Call initialisation script for initialising rest of the parameters and 
%    return an initialisation data structure.
init_data=initialisation(GPU_flag,NT_steady_state,krefL,krefW,nref,sigma,...
    z_domain_size_factor,no_of_PSD_cycles);
%--------------------------- STM initialisation ----------------------------
del_eps=0.5; % dielectric constant modulation by a space-time modulated curve  
spat_freq_fact=0.4;          % spatial frequency of modulation
f_mod=0.2*(init_data.fsrc);  % temporal frequency of the modulation
                             % keep the fraction less than or equal to 0.2
grat_width=ceil(3.0*pi/(init_data.kref*init_data.dz)); % space-time modulated 
                                                       %      grating width
%------------- Define the source  -----------------------------------------
[source_data] = define_dipole_beam_source(GPU_flag,init_data); 
                                       % Returns the source data structure.
source_data.clip_lim=3;  % Colorbar clip limit for better contrast

source_data.t_sign=+1; % +1 for forward propagating source and 
                       % -1 for time reversed version
src_z_start=floor(init_data.start_slab1/2);  % Half way between slab and the
                                            % left boundary
source_data.src_z=[src_z_start src_z_start+1];     % Dipole source location
source_data.TFSF_Nz=floor(source_data.src_z(1)/2); % TFSF boundary location 
                        % Half way between the source and the left boundary 
source_data.TFSF_region=1:source_data.TFSF_Nz;
%------------------- STM parameters ----------------------------------------
[STM_data]=define_STM_parameters(GPU_flag,init_data,source_data,del_eps,...
                                spat_freq_fact,f_mod,grat_width);
%------------ Forward FDTD wave propagation -------------------------------
tic
[Efield_array_collection] = FDTD_main_loop_STM(GPU_flag,visualisation_flag,...
                             init_data,source_data,STM_data);

%---- Phase sensitive detection for the harmonic (sideband) wavefronts ----
init_data.detec_trans_position=...
 floor(((init_data.start_slab1 + init_data.scat_width) + init_data.Nz)/2);
                        % Transmission measurement z index location
                        % Half way between the slab and the right boundary                              
init_data.detec_refl_position=floor(source_data.TFSF_Nz/2); 
                         % Reflection measurement z index location
                % Half way between the TF/SF boundary and the left boundary                                                            
init_data.detec_src_position=STM_data.STM_start_z;
                         % Source transmission measurement z index location

norder_array=[+3,+2,+1,0,-1];  % Various modulated temporal order number
for ncount=1:length(norder_array)
nth=norder_array(ncount);
PSD_array_collection=PSD_nth_SideBand(GPU_flag,visualisation_flag,...
    Efield_array_collection,init_data,source_data,STM_data,nth);
title(sprintf('$Phase~sensitive~detection~: \\bf n=%d$',nth),'Interpreter','Latex')
PSD_structure_collection(ncount)=PSD_array_collection;
end
%--------------------------------------------------------------------------
%save data
if data_saving_flag==1
save('PSD_structure_collection.mat','PSD_structure_collection')
save('Efield_array_collection_STM.mat','Efield_array_collection')
save('STM_data.mat','STM_data')
save('init_data_forward.mat','init_data')
end

% figure(3)
% xlabel('$k_{ref}z$','Interpreter','Latex')
% ylabel('$k_{ref}y$','Interpreter','Latex')
% set(gca,'FontSize',24);
% title('$Raman-Nath~regime$','Interpreter','Latex')
% exportgraphics(gcf,'Bragg_forward.eps','ContentType','vector')
saveas(gcf,'Raman-Nath-forward.fig')
%--------- plot the extracted sideband wavefronts -------------------------
krefy_detec=(0:init_data.Ny+1).*(init_data.kref*init_data.dy);

%--- find the peak values on the modulated wavefronts for normalisation
peak_val_store=zeros(1,length(norder_array));
for ncount=1:length(norder_array)
peak_val_store(ncount)=max(PSD_structure_collection(ncount).E_trans_sig);
end

figure('position',[70 70 500 800],'color','white')
xlabel('Transmitted~intensity~(dB)','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
set(gca,'FontSize',24)
for ncount=1:length(norder_array)
nth=norder_array(ncount);
plot(mag2db(PSD_structure_collection(ncount).E_trans_sig./max(peak_val_store)),...
    krefy_detec,'-*','MarkerSize',1)
title(sprintf('$Temporal~order~no,~n = %d$',nth),'Interpreter','Latex')
xlabel('$Transmitted~intensity~(dB)$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
ylim([0 init_data.kref*init_data.dy*(init_data.Ny-1)])
set(gca,'FontSize',22)
hold on
pause(0.5)
end

for lcount=1:numel(norder_array)
    if norder_array(lcount)>0
    legend_array{lcount}=sprintf('$n=+%d$',norder_array(lcount));
    elseif norder_array(lcount)<=0
    legend_array{lcount}=sprintf('$n=%d$',norder_array(lcount));
    end
end
legend(legend_array,'Interpreter','Latex');
title('');
%----------------------- STM parameters for the paper ---------------------
kmodL=(STM_data.grat_width-1)*STM_data.kmod*init_data.dz;
krefL=(STM_data.grat_width-1)*init_data.kref*init_data.dz;
STM_data.kmod/init_data.kref;
Q=(STM_data.kmod/init_data.kref)*(kmodL);
nu=krefL*STM_data.del_eps;
STM_data.omega_mod/(2*pi*init_data.fsrc);
%---------------------- Annotate --------------------
open('Raman-Nath-forward.fig')
side_description={sprintf('$k_{mod}L=%.2f$',kmodL),...
                  sprintf('$k_{ref}L=%.2f$',krefL),...
                  sprintf('$k_{mod}/k_{ref}=%.2f$',STM_data.kmod/init_data.kref),...
                  sprintf('$Q=%.3f$',Q),...
                  sprintf('$\\nu=%.3f$',nu),...
                  sprintf('$f_{mod}/f_{src}=%.3f$',STM_data.f_mod/init_data.fsrc)};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 
pause(1);