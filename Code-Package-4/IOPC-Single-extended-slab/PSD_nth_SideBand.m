function [PSD_array_collection] = PSD_nth_SideBand(GPU_flag,visualisation_flag,...
    Efield_array_collection,init_data,source_data,STM_data,nth)

Enm1=Efield_array_collection.Enm1;
En=Efield_array_collection.En;
Enp1=Efield_array_collection.Enp1;
Esrc_nm1=Efield_array_collection.Esrc_nm1;
Esrc_n=Efield_array_collection.Esrc_n;
Esrc_np1=Efield_array_collection.Esrc_np1;
eps_avg=init_data.eps_avg;
%--------------------------------------------------------------------------
start_slab1=init_data.start_slab1;
start_slab2=init_data.start_slab2;
%--------------------------------------------------------------------------
c0=init_data.c0;
Ny=init_data.Ny;
Nz=init_data.Nz;
kref=init_data.kref;
dz=init_data.dz;
dy=init_data.dy;
fsrc=init_data.fsrc;
dt=init_data.dt;
mirror_pos=init_data.mirror_pos;
clip_lim=source_data.clip_lim;
NT_steady_state=init_data.NT_steady_state;
wsrc0=init_data.wsrc0;
epsilon_ij=init_data.epsilon_ij;
scat_width=init_data.scat_width;
fsrc=init_data.fsrc;
%---------------------  source function -----------------------------------
src_y=source_data.src_y;         % Source y location  
src_z=source_data.src_z;         % Source z location
source_function=source_data.dipole_source;
wave_amp=source_data.wave_amp;
wave_phase=source_data.wave_phase;
t_sign=source_data.t_sign;
TFSF_Nz=source_data.TFSF_Nz;
TFSF_region=source_data.TFSF_region;

no_of_PSD_cycles=init_data.no_of_PSD_cycles;
%------------- For the Absorbing BC ---------------------------------------
cref=c0/sqrt(eps_avg);
p0=1.00023;
p2=0.51555;
if GPU_flag==1
term1=gpuArray((p2*(dz*(cref*dt)^2)/((dy^2)*(cref*dt+p0*dz))));
term2=gpuArray(((2*p0*dz)/(cref*dt+p0*dz)));
term3=gpuArray((cref*dt-p0*dz)/(cref*dt+p0*dz));
else
term1=((p2*(dz*(cref*dt)^2)/((dy^2)*(cref*dt+p0*dz))));
term2=(((2*p0*dz)/(cref*dt+p0*dz)));  
term3=(cref*dt-p0*dz)/(cref*dt+p0*dz);
end
%----------------------- STM parameters -----------------------------------
if GPU_flag==1
f_mod=gather(STM_data.f_mod);
else
f_mod=STM_data.f_mod;
end
grat_width=STM_data.grat_width;
STM_start_z=STM_data.STM_start_z;
STM_start_y=STM_data.STM_start_y;
STM_end_y=STM_data.STM_end_y;
eps_ij_n=STM_data.eps_ij_n;
Dneps_ij_n=STM_data.Dneps_ij_n;
D2neps_ij_n=STM_data.D2neps_ij_n;

detec_trans_position=init_data.detec_trans_position;
detec_refl_position=init_data.detec_refl_position;
%--------  step 2 : Phase sensitive detection (PSD) ----------------------- 
fsrc_SB=fsrc+nth*f_mod;  % side band frequency taken 
Tb=no_of_PSD_cycles*(1/(f_mod)); % Beat cycles
time_steps_per_beat=ceil(Tb/dt);
t_data=zeros(1,time_steps_per_beat);
Eref_0=zeros(1,time_steps_per_beat);  % In phase component
Eref_90=zeros(1,time_steps_per_beat); % Quadrature component

if GPU_flag==1
En_trans_PSD=gpuArray(zeros(Ny+2,time_steps_per_beat));% Transmitted wave stored
En_refl_PSD=gpuArray(zeros(Ny+2,time_steps_per_beat));% Reflected wave stored
else
En_trans_PSD=zeros(Ny+2,time_steps_per_beat);  % Transmitted wave stored
En_refl_PSD=zeros(Ny+2,time_steps_per_beat);   % Reflected wave stored
end

cy_count=0;
Eref=1;                                 % Reference amplitude
theta_ref=0;                            % Reference phase

if visualisation_flag==0
figure('Position', [50 50 1200 900],'color','white');
colormap wavecolormap   
hIm=imagesc([0 kref*dz*(Nz-1)],[0 kref*dy*(Ny-1)],En(2:end-1,2:end-1));
axis xy equal tight;
xlabel('$k_{ref} z$','Interpreter','Latex')
ylabel('$k_{ref} y$','Interpreter','Latex')
colorbar
caxis([-clip_lim clip_lim])
hold on
line([start_slab1*kref*dz (start_slab1+scat_width-1)*kref*dz (start_slab1+ ...
    scat_width-1)*kref*dz  start_slab1*kref*dz  start_slab1*kref*dz ],[0 0 ...
    Ny*kref*dy  Ny*kref*dy  0],'LineWidth',2);
    % Plot the boundary of the disorder
hold on
line([start_slab2*kref*dz (start_slab2+scat_width-1)*kref*dz (start_slab2+ ...
    scat_width-1)*kref*dz  start_slab2*kref*dz  start_slab2*kref*dz ],[0 0 ...
    Ny*kref*dy  Ny*kref*dy  0],'LineWidth',2);
    % Plot the boundary of the disorder
hold on
line([STM_start_z*kref*dz (STM_start_z+grat_width-1)*kref*dz ...
    (STM_start_z+grat_width-1)*kref*dz  STM_start_z*kref*dz ...
    STM_start_z*kref*dz ],[STM_start_y*kref*dy STM_start_y*kref*dy ...
    STM_end_y*kref*dy  STM_end_y*kref*dy  STM_start_y*kref*dy],'LineWidth',2)
hold on
line([TFSF_Nz*kref*dz  TFSF_Nz*kref*dz ],[0 Ny*kref*dy], ...
    'color','w','LineWidth',2);
    % This is the TFSF boundary plotted on the left side of the source.
hold on
line([(src_z(1)-1)*kref*dz  (src_z(1)-1)*kref*dz ],...
    [0 Ny*kref*dy],'color','w','LineWidth',2)    
set(gca,'FontSize',22)
end

for t_count=(NT_steady_state+1) : (NT_steady_state+1) + time_steps_per_beat-1
t_count
%---------------------- inject wave source ---------------------------------
En(src_y,src_z)=En(src_y,src_z)+ ...
    source_function(t_count,wave_amp,wave_phase,t_sign,fsrc); % Soft dipole source
Esrc_n(src_y,src_z)=Esrc_n(src_y,src_z)+ ...
    source_function(t_count,wave_amp,wave_phase,t_sign,fsrc); % soft dipole source
% Enforce the transverse mirror boundary condition
En(mirror_pos,:)=0;    
En(end-mirror_pos,:)=0;   
Esrc_n(mirror_pos,:)=0;
Esrc_n(end-mirror_pos,:)=0;
%               
%--------------------field update equations : new version -----------------
% Field at next time step is a function of the field at the current time
% step and the field at the previous time step
 W0=eps_ij_n(t_count).*(dz/(c0*dt))^2;
 W1=Dneps_ij_n(t_count).*((dz/c0)^2).*(1/dt);
 W2=D2neps_ij_n(t_count).*((dz/c0)^2);


Enp1(2:end-1,2:end-1)= (1./(W0+W1)).*(En(2:end-1,3:end)-2.*En(2:end-1,2:end-1)+ En(2:end-1,1:end-2) + ...
                                      En(3:end,2:end-1)-2.*En(2:end-1,2:end-1)+ En(1:end-2,2:end-1)) + ...
              ((2.*W0-W2)./(W0+W1)).*En(2:end-1,2:end-1)+ ...
                -((W0-W1)./(W0+W1)).*Enm1(2:end-1,2:end-1);


Esrc_np1(2:end-1,2:end-1)=2.*Esrc_n(2:end-1,2:end-1)-Esrc_nm1(2:end-1,2:end-1) + ...
       wsrc0.*(Esrc_n(2:end-1,3:end)-2.*Esrc_n(2:end-1,2:end-1)+ Esrc_n(2:end-1,1:end-2)) + ...
       wsrc0.*(Esrc_n(3:end,2:end-1)-2.*Esrc_n(2:end-1,2:end-1)+ Esrc_n(1:end-2,2:end-1));
%------------------ Mur ABC -----------------------------------------------
[Enp1,Esrc_np1] = ...
MUR_abc_2nd_order(Enm1,En,Enp1,Esrc_nm1,Esrc_n,Esrc_np1,term1,term2,term3);
%---------- Estimate the Total-field/Scatter field representation by 
%           subtracting the source on the left side of the boundary -------
EnTFSF=En;
EnTFSF(:,TFSF_region)=En(:,TFSF_region)-Esrc_n(:,TFSF_region);
%------------------------ plot -------------------------------------------
hIm.CData = gather(EnTFSF(2:end-1,2:end-1));
drawnow   
%---------------------------- Update the fields ---------------------------
Enm1=En;
En=Enp1;
Esrc_nm1=Esrc_n;
Esrc_n=Esrc_np1;
%---------------------- PSD detection ---------------------------------
cy_count=cy_count+1;

t_data(cy_count)=(t_count-1)*dt;
%----- Note that the in-phase and the quadrature phase components
% have fsrc_SB as its frequency ------------------------------------------

Eref_0(cy_count)=Eref*cos(2*pi*fsrc_SB*t_data(cy_count)+theta_ref);
Eref_90(cy_count)=Eref*cos(2*pi*fsrc_SB*t_data(cy_count)+theta_ref+pi/2);

En_trans_PSD(:,cy_count)=(En(:,detec_trans_position));  
%%% detection at detec_trans_position
%---------------------- Record the reflection for 1 beat period ---------
En_refl_PSD(:,cy_count)=(En(:,detec_refl_position));  
%%% detection at detec_refl_position
%--------------------------------------------------------------------------
end
%--------------- PSD amplitude and phase estimation -----------------------
% Estimate transmission
Xtrans=(1/Tb)*(sum(En_trans_PSD.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Ytrans=(1/Tb)*(sum(En_trans_PSD.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_trans_sig=2*sqrt(Xtrans.^2+Ytrans.^2)/Eref;
E_theta_trans_sig=atan2(Ytrans,Xtrans) + theta_ref;
%------------ PSD reflection amplitude and phase estimation ---------------

%-------------------- PSD row-wise : transmission  ------------------------
Xrefl=(1/Tb)*(sum(En_refl_PSD.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Yrefl=(1/Tb)*(sum(En_refl_PSD.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_refl_sig=2*sqrt(Xrefl.^2+Yrefl.^2)/Eref;
E_theta_refl_sig=atan2(Yrefl,Xrefl) + theta_ref;


if GPU_flag==1
En=gather(En);
end
save('En_L2R.mat','En');
save('epsilon_ij.mat','epsilon_ij');

%------------------- define a structure to return the data ----------------
if GPU_flag==1
PSD_array_collection.E_trans_sig=gather(E_trans_sig);
PSD_array_collection.E_theta_trans_sig=gather(E_theta_trans_sig);
PSD_array_collection.E_refl_sig=gather(E_refl_sig);
PSD_array_collection.E_theta_refl_sig=gather(E_theta_refl_sig);
else 
PSD_array_collection.E_trans_sig=E_trans_sig;
PSD_array_collection.E_theta_trans_sig=E_theta_trans_sig;
PSD_array_collection.E_refl_sig=E_refl_sig;
PSD_array_collection.E_theta_refl_sig=E_theta_refl_sig;
end
%-------------------------------------------------------------------------
end

