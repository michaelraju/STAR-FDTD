function [PSD_array_collection] = phase_sensitive_detection(GPU_flag,...
    visualisation_flag,Efield_array_collection,init_data,source_data)

%-------------- Extract the data from the structures ----------------------
Enm1=Efield_array_collection.Enm1;
En=Efield_array_collection.En;
Enp1=Efield_array_collection.Enp1;
Esrc_nm1=Efield_array_collection.Esrc_nm1;
Esrc_n=Efield_array_collection.Esrc_n;
Esrc_np1=Efield_array_collection.Esrc_np1;
EnTFSF=Efield_array_collection.EnTFSF;
c0=init_data.c0;
Ny=init_data.Ny;
Nz=init_data.Nz;
kref=init_data.kref;
dz=init_data.dz;
dy=init_data.dy;
fsrc=init_data.fsrc;
dt=init_data.dt;
start_slab1=init_data.start_slab1;
scat_width=init_data.scat_width;
mirror_pos=init_data.mirror_pos;
clip_lim=source_data.clip_lim;
NT_steady_state=init_data.NT_steady_state;
w0=init_data.w0;
wsrc0=init_data.wsrc0;
eps_avg=init_data.eps_avg;
detec_trans_position=init_data.detec_trans_position;
detec_refl_position=init_data.detec_refl_position;
detec_src_position=init_data.detec_src_position;
no_of_PSD_cycles=init_data.no_of_PSD_cycles;
%------------- Data for the generalised Mur absorbing BC ------------------
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
%-------------------------- Source ----------------------------------------
source_function=source_data.dipole_source;
wave_amp=source_data.wave_amp;
wave_phase=source_data.wave_phase;
t_sign=source_data.t_sign;
TFSF_Nz=source_data.TFSF_Nz;
src_y=source_data.src_y;
src_z=source_data.src_z;
TFSF_region=source_data.TFSF_region;
%--------  step 2 : Phase sensitive detection (PSD) ----------------------- 
Tb=no_of_PSD_cycles*(1/(fsrc));       % many cycles of one Time period
time_steps_per_cycle=ceil(Tb/dt);
t_data=zeros(1,time_steps_per_cycle);
Eref_0=zeros(1,time_steps_per_cycle);  % In-phase component
Eref_90=zeros(1,time_steps_per_cycle); % Quadrature component

if GPU_flag==1
En_trans_PSD=gpuArray(zeros(Ny+2,time_steps_per_cycle));  
                                % Transmitted wave stored
En_refl_PSD=gpuArray(zeros(Ny+2,time_steps_per_cycle));   
                                 % Reflected wave stored
En_source_PSD=gpuArray(zeros(Ny+2,time_steps_per_cycle)); 
                                 % Source wave stored

En_trans_PSD_plus_1=gpuArray(zeros(Ny+2,time_steps_per_cycle));  
                                  % Transmitted wave stored
En_refl_PSD_plus_1=gpuArray(zeros(Ny+2,time_steps_per_cycle));   
                                  % Reflected wave stored
En_source_PSD_plus_1=gpuArray(zeros(Ny+2,time_steps_per_cycle)); 
                                  % Source wave stored

En_trans_PSD_minus_1=gpuArray(zeros(Ny+2,time_steps_per_cycle));  
                                   % Transmitted wave stored
En_refl_PSD_minus_1=gpuArray(zeros(Ny+2,time_steps_per_cycle));   
                                   % Reflected wave stored
En_source_PSD_minus_1=gpuArray(zeros(Ny+2,time_steps_per_cycle)); 
                                   % Source wave stored

else
En_trans_PSD=zeros(Ny+2,time_steps_per_cycle);  % Transmitted wave stored
En_refl_PSD=zeros(Ny+2,time_steps_per_cycle);   % Reflected wave stored
En_source_PSD=zeros(Ny+2,time_steps_per_cycle); % Source wave stored

En_trans_PSD_plus_1=(zeros(Ny+2,time_steps_per_cycle));  
                                                % Transmitted wave stored
En_refl_PSD_plus_1=(zeros(Ny+2,time_steps_per_cycle));   
                                                % Reflected wave stored
En_source_PSD_plus_1=(zeros(Ny+2,time_steps_per_cycle)); 
                                                % Source wave stored

En_trans_PSD_minus_1=(zeros(Ny+2,time_steps_per_cycle));  
                                               % Transmitted wave stored
En_refl_PSD_minus_1=(zeros(Ny+2,time_steps_per_cycle));   
                                               % Reflected wave stored
En_source_PSD_minus_1=(zeros(Ny+2,time_steps_per_cycle)); 
                                               % Source wave stored
end

cy_count=0;
Eref=1;                                 % Reference amplitude for PSD
theta_ref=0;                            % Reference phase for PSD

if visualisation_flag==0
figure('Position', [0 0 1500 1200]);
end
colormap wavecolormap   
hIm=imagesc([0 kref*dz*(Nz-1)],[0 kref*dy*(Ny-1)],EnTFSF(2:end-1,2:end-1));
axis xy; axis equal; axis tight
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
colorbar
caxis([-clip_lim clip_lim])
hold on
line([start_slab1*kref*dz (start_slab1+scat_width-1)*kref*dz  (start_slab1+ ...
    scat_width-1)*kref*dz  start_slab1*kref*dz  start_slab1*kref*dz ],[0 0 ...
    Ny*kref*dy  Ny*kref*dy  0],'LineWidth',2)
hold on
line([(TFSF_Nz-1)*kref*dz  (TFSF_Nz-1)*kref*dz ],[0 Ny*kref*dy],...
    'color','w','LineWidth',2)
hold on
line([(src_z(1)-1)*kref*dz  (src_z(1)-1)*kref*dz ],[0 Ny*kref*dy],...
    'color','w','LineWidth',2)
set(gca,'FontSize',22)

%------------- FDTD main Loop after steady state --------------------------
for t_count=(NT_steady_state+1) : (NT_steady_state+1) + time_steps_per_cycle-1
t_count
%---------------------- inject PC wave source -----------------------------
En(src_y,src_z)=En(src_y,src_z)+ ...
    source_function(t_count,wave_amp,wave_phase,t_sign); % Soft dipole source
Esrc_n(src_y,src_z)=Esrc_n(src_y,src_z)+...
    source_function(t_count,wave_amp,wave_phase,t_sign); % soft dipole source
%--------------------------------------------------------------------------
% Enforce the transverse mirror boundary condition
En(mirror_pos,:)=0;    
En(end-mirror_pos,:)=0;   
Esrc_n(mirror_pos,:)=0;
Esrc_n(end-mirror_pos,:)=0;             
%--------------------field update equations -------------------------------
% Field at next time step is a function of the field at the current time
% step and the field at the previous time step
Enp1(2:end-1,2:end-1)=2.*En(2:end-1,2:end-1)-Enm1(2:end-1,2:end-1) + ...
       w0.*(En(2:end-1,3:end)-2.*En(2:end-1,2:end-1)+ En(2:end-1,1:end-2)) + ...
       w0.*(En(3:end,2:end-1)-2.*En(2:end-1,2:end-1)+ En(1:end-2,2:end-1));
   
Esrc_np1(2:end-1,2:end-1)=2.*Esrc_n(2:end-1,2:end-1)-Esrc_nm1(2:end-1,2:end-1) + ...
       wsrc0.*(Esrc_n(2:end-1,3:end)-2.*Esrc_n(2:end-1,2:end-1)+ Esrc_n(2:end-1,1:end-2)) + ...
       wsrc0.*(Esrc_n(3:end,2:end-1)-2.*Esrc_n(2:end-1,2:end-1)+ Esrc_n(1:end-2,2:end-1));
%------------------ Mur ABC -----------------------------------------------          
[Enp1,Esrc_np1] = MUR_abc_2nd_order(Enm1,En,Enp1,Esrc_nm1,Esrc_n,Esrc_np1,term1,term2,term3);
%---------- Estimate the Total-field/Scatter field representation by 
%  subtracting the source on the left side of the boundary ----------------
EnTFSF=En;
EnTFSF(:,TFSF_region)=En(:,TFSF_region)-Esrc_n(:,TFSF_region);
%------------------------ plot --------------------------------------------
hIm.CData = gather(EnTFSF(2:end-1,2:end-1));  % Plot the TF/SF
%hIm.CData = gather(Esrc_n(2:end-1,2:end-1)); % Plot the source field
drawnow
%---------------------- PSD detection -------------------------------------
cy_count=cy_count+1;
t_data(cy_count)=(t_count-1)*dt;
Eref_0(cy_count)=Eref*cos(2*pi*fsrc*t_data(cy_count)+theta_ref);
Eref_90(cy_count)=Eref*cos(2*pi*fsrc*t_data(cy_count)+theta_ref+pi/2);
%---------------------- Record the transmission ---------------------------
En_trans_PSD(:,cy_count)=(EnTFSF(:,detec_trans_position));% detection at 
                              %                        detec_trans_position 
En_trans_PSD_plus_1(:,cy_count)=(EnTFSF(:,detec_trans_position+1)); 
                              % detection at 
                              %   detec_trans_position + 1 (for evaluating
                              %   gradient for estimating transmission and 
                              %   reflection later)
En_trans_PSD_minus_1(:,cy_count)=(EnTFSF(:,detec_trans_position-1));  
                              % detection at 
                              %       detec_trans_position - 1
%---------------------- Record the reflection -----------------------------
En_refl_PSD(:,cy_count)=(EnTFSF(:,detec_refl_position));  
En_refl_PSD_plus_1(:,cy_count)=(EnTFSF(:,detec_refl_position+1));   
En_refl_PSD_minus_1(:,cy_count)=(EnTFSF(:,detec_refl_position-1));  
%---------------------- Record the source transmission --------------------
En_source_PSD(:,cy_count)=(Esrc_n(:,detec_src_position));  %%% detection at
                                     % 1 cell after the TF/SF source boundary 
En_source_PSD_plus_1(:,cy_count)=(Esrc_n(:,detec_src_position+1));   
En_source_PSD_minus_1(:,cy_count)=(Esrc_n(:,detec_src_position-1));  
%------------------ Update the fields -------------------------------------
Enm1=En;
En=Enp1;
Esrc_nm1=Esrc_n;
Esrc_n=Esrc_np1;
end
hfig=gcf;
%----- Phase sensitive detection of amplitude and phase -------------------
% Estimate transmitted wave amplitude and phase
Xtrans=(1/Tb)*(sum(En_trans_PSD.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Ytrans=(1/Tb)*(sum(En_trans_PSD.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_trans=2*sqrt(Xtrans.^2+Ytrans.^2)/Eref;
theta_trans=atan2(Ytrans,Xtrans) + theta_ref;

Xtrans=(1/Tb)*(sum(En_trans_PSD_plus_1.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Ytrans=(1/Tb)*(sum(En_trans_PSD_plus_1.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_trans_plus_1=2*sqrt(Xtrans.^2+Ytrans.^2)/Eref;
theta_trans_plus_1=atan2(Ytrans,Xtrans) + theta_ref;

Xtrans=(1/Tb)*(sum(En_trans_PSD_minus_1.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Ytrans=(1/Tb)*(sum(En_trans_PSD_minus_1.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_trans_minus_1=2*sqrt(Xtrans.^2+Ytrans.^2)/Eref;
theta_trans_minus_1=atan2(Ytrans,Xtrans) + theta_ref;
%--------------------------------------------------------------------------
figure('position',[100 100 1200 900])
subplot(1,2,1)    
plot(E_trans,(0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',0.5)
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
title('$Transmission~Magnitude$','Interpreter','Latex')
set(gca,'FontSize',22)

subplot(1,2,2)    
plot(theta_trans,(0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',0.5)
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
title('$Transmission~Phase$','Interpreter','Latex')
set(gca,'FontSize',22)

% Estimate reflected wave amplitude and phase
Xrefl=(1/Tb)*(sum(En_refl_PSD.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Yrefl=(1/Tb)*(sum(En_refl_PSD.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_refl=2*sqrt(Xrefl.^2+Yrefl.^2)/Eref;
theta_refl=atan2(Yrefl,Xrefl) + theta_ref;

Xrefl=(1/Tb)*(sum(En_refl_PSD_plus_1.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Yrefl=(1/Tb)*(sum(En_refl_PSD_plus_1.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_refl_plus_1=2*sqrt(Xrefl.^2+Yrefl.^2)/Eref;
theta_refl_plus_1=atan2(Yrefl,Xrefl) + theta_ref;

Xrefl=(1/Tb)*(sum(En_refl_PSD_minus_1.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Yrefl=(1/Tb)*(sum(En_refl_PSD_minus_1.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_refl_minus_1=2*sqrt(Xrefl.^2+Yrefl.^2)/Eref;
theta_refl_minus_1=atan2(Yrefl,Xrefl) + theta_ref;

figure('position',[100 100 1200 900])
subplot(1,2,1)    
plot(E_refl,(0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',0.5)
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
title('$Reflection~Magnitude$','Interpreter','Latex')
set(gca,'FontSize',22)

subplot(1,2,2)    
plot(theta_refl,(0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',0.5)
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
title('$Reflection~Phase$','Interpreter','Latex')
set(gca,'FontSize',22)

% Estimate source wave amplitude and phase
Xsource=(1/Tb)*(sum(En_source_PSD.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Ysource=(1/Tb)*(sum(En_source_PSD.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_source=2*sqrt(Xsource.^2+Ysource.^2)/Eref;
theta_source=atan2(Ysource,Xsource) + theta_ref;

Xsource=(1/Tb)*(sum(En_source_PSD_plus_1.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Ysource=(1/Tb)*(sum(En_source_PSD_plus_1.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_source_plus_1=2*sqrt(Xsource.^2+Ysource.^2)/Eref;
theta_source_plus_1=atan2(Ysource,Xsource) + theta_ref;

Xsource=(1/Tb)*(sum(En_source_PSD_minus_1.*repmat(Eref_0,[Ny+2,1]).*dt,2));
Ysource=(1/Tb)*(sum(En_source_PSD_minus_1.*repmat(Eref_90,[Ny+2,1]).*dt,2));
E_source_minus_1=2*sqrt(Xsource.^2+Ysource.^2)/Eref;
theta_source_minus_1=atan2(Ysource,Xsource) + theta_ref;

%--------------------------------------------------------------------------
figure('position',[100 100 1200 900])
subplot(1,2,1)    
plot(E_source,(0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',0.5)
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
title('$Source~Magnitude$','Interpreter','Latex')
set(gca,'FontSize',22)

subplot(1,2,2)    
plot(theta_source,(0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',0.5)
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
title('$Source~Phase$','Interpreter','Latex')
set(gca,'FontSize',22)
%------ Integrating the wave current for one PSD cycle for estimating 
%       the total transmission and reflection

%------------ Define wave current -----------------------------------------
omega0=2*pi*fsrc;
Jz_boundary= @(t,Ez0,Ez_plus_1,Ez_minus_1,theta_z0,theta_zplus_1,theta_zminus_1) ...
    (omega0/2*dz).*(Ez0.*sin(omega0.*t+theta_z0).*(Ez_plus_1.*cos(omega0.*t+theta_zplus_1)-...
    Ez_minus_1.*cos(omega0.*t+theta_zminus_1)));

Jz_trans_time_series=Jz_boundary(t_data,E_trans,E_trans_plus_1,E_trans_minus_1,...
    theta_trans,theta_trans_plus_1,theta_trans_minus_1);
Jz_trans_integrated=zeros(Ny,1);

Jz_refl_time_series=Jz_boundary(t_data,E_refl,E_refl_plus_1,E_refl_minus_1,...
    theta_refl,theta_refl_plus_1,theta_refl_minus_1);
Jz_refl_integrated=zeros(Ny,1);

Jz_source_time_series=Jz_boundary(t_data,E_source,E_source_plus_1,E_source_minus_1,...
    theta_source,theta_source_plus_1,theta_source_minus_1);
Jz_source_integrated=zeros(Ny,1);

for ycount=1:Ny % Integrating over time to obtain spatially resolved average current
Jz_trans_integrated(ycount)=(1/Tb)*trapz(t_data,Jz_trans_time_series(ycount,:));
Jz_refl_integrated(ycount)=-(1/Tb)*trapz(t_data,Jz_refl_time_series(ycount,:));
Jz_source_integrated(ycount)=(1/Tb)*trapz(t_data,Jz_source_time_series(ycount,:));
end
%----------- Estimate total transmission and reflection -------------------
y_detec=(0:Ny-1).*dy;
T=trapz(y_detec,Jz_trans_integrated)/trapz(y_detec,Jz_source_integrated);% normalising wrt source 
R=trapz(y_detec,Jz_refl_integrated)/trapz(y_detec,Jz_source_integrated);% normalising wrt source 
sprintf('Transmission = %f, reflection = %f and T+R= %f',T,R,T+R)
%------------------- define a structure to return the data ----------------
PSD_array_collection.E_trans_sig=E_trans;
PSD_array_collection.theta_trans_sig=theta_trans;
PSD_array_collection.E_refl_sig=E_refl;
PSD_array_collection.theta_refl_sig=theta_refl;
PSD_array_collection.E_source_sig=E_source;
PSD_array_collection.theta_source_sig=theta_source;
PSD_array_collection.R=R;
PSD_array_collection.T=T;
%------------------------------ Annotate ----------------------------------
figure(hfig)
side_description={sprintf('$T=%.5f$',PSD_array_collection.T),...
                  sprintf('$R=%.5f$',PSD_array_collection.R),...
                  sprintf('$T+R=%.5f$',PSD_array_collection.T+PSD_array_collection.R),...
                  sprintf('$Boundary~z~positions$'),...
                  sprintf('$k_{ref}z|_{Refl}=%.3f$',(init_data.kref*init_data.dz*(detec_refl_position-1))),...
                  sprintf('$k_{ref}z|_{TFSF}=%.3f$',(init_data.kref*init_data.dz*(TFSF_Nz-1))),...
                  sprintf('$k_{ref}z|_{src}=%.3f$',(init_data.kref*init_data.dz*(src_z(1)-1))),...
                  sprintf('$k_{ref}z|_{Trans}=%.3f$',(init_data.kref*init_data.dz*(detec_trans_position-1))),...
                  };
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 
drawnow;
end

