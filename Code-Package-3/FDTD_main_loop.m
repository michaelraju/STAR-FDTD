function [Efield_array_collection] = FDTD_main_loop(GPU_flag,visualisation_flag,init_data,source_data)
%----- Extract the parameters from the structures -------------------------
Ny=init_data.Ny;
Nz=init_data.Nz;
kref=init_data.kref;
dz=init_data.dz;
dy=init_data.dy;
dt=init_data.dt;
start_slab1=init_data.start_slab1;
scat_width=init_data.scat_width;
mirror_pos=init_data.mirror_pos;
clip_lim=source_data.clip_lim;
NT_steady_state=init_data.NT_steady_state;
w0=init_data.w0;
wsrc0=init_data.wsrc0;
c0=init_data.c0;
eps_avg=init_data.eps_avg;
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

%---------------------  Initialise field matrices -------------------------
if GPU_flag==1
Enp1=gpuArray(zeros(Ny+2,Nz+2));  % Field at n+1 th time step
En=gpuArray(zeros(Ny+2,Nz+2));    % Field at nth time step
Enm1=gpuArray(zeros(Ny+2,Nz+2));  % Field at n-1 th time step
Esrc_np1=gpuArray(zeros(Ny+2,Nz+2));% Source only field at n+1 th time step
Esrc_n=gpuArray(zeros(Ny+2,Nz+2));  % Source only field at n th time step
Esrc_nm1=gpuArray(zeros(Ny+2,Nz+2));% Source only field at n-1 th time step
                    % to be used in the source only field update equations
else 
% defining the same matrices but without using gpuArray
Enp1=(zeros(Ny+2,Nz+2));
En=(zeros(Ny+2,Nz+2));
Enm1=(zeros(Ny+2,Nz+2));
Esrc_np1=(zeros(Ny+2,Nz+2));
Esrc_n=(zeros(Ny+2,Nz+2));
Esrc_nm1=(zeros(Ny+2,Nz+2));
end
%------------- FDTD main Loop to attain steady state ----------------------
% If visualisation is switched on, the code execution will be slower
if visualisation_flag==1
figure('Position', [0 0 800 1500]);
colormap wavecolormap   % Colormap for visualising wave fields
hIm=imagesc([0 kref*dz*(Nz-1)],[0 kref*dy*(Ny-1)],En(2:end-1,2:end-1));
axis xy equal tight; 
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
colorbar
caxis([-clip_lim clip_lim])
hold on
line([start_slab1*kref*dz (start_slab1+scat_width-1)*kref*dz (start_slab1+ ...
    scat_width-1)*kref*dz  start_slab1*kref*dz  start_slab1*kref*dz ],[0 0 ...
    Ny*kref*dy  Ny*kref*dy  0],'LineWidth',2);
    % Plot the boundary of the disorder
hold on
line([(TFSF_Nz-1)*kref*dz  (TFSF_Nz-1)*kref*dz ],[0 Ny*kref*dy], ...
    'color','w','LineWidth',2);
    % This is the TFSF boundary plotted on the left side of the source.
hold on
line([(src_z(1)-1)*kref*dz  (src_z(1)-1)*kref*dz ],[0 Ny*kref*dy],...
'color','w','LineWidth',2)
set(gca,'FontSize',22)
end

for t_count=1:NT_steady_state
sprintf('Time count = %d/%d',t_count,NT_steady_state)
%---------------------- inject wave source ---------------------------------
En(src_y,src_z)=En(src_y,src_z)+ ...
    source_function(t_count,wave_amp,wave_phase,t_sign,fsrc); % Soft dipole source
Esrc_n(src_y,src_z)=Esrc_n(src_y,src_z)+ ...
    source_function(t_count,wave_amp,wave_phase,t_sign,fsrc); % soft dipole source
%------------- Enforce the transverse mirror boundary condition -----------
En(mirror_pos,:)=0;    
En(end-mirror_pos,:)=0;   
Esrc_n(mirror_pos,:)=0;
Esrc_n(end-mirror_pos,:)=0;
%--------------------field update equations -------------------------------
% Field at next time step is a function of the field at the current time
%           step and the field at the previous time step.
Enp1(2:end-1,2:end-1)=2.*En(2:end-1,2:end-1)-Enm1(2:end-1,2:end-1) + ...
       w0.*(En(2:end-1,3:end)-2.*En(2:end-1,2:end-1)+ En(2:end-1,1:end-2)) + ...
       w0.*(En(3:end,2:end-1)-2.*En(2:end-1,2:end-1)+ En(1:end-2,2:end-1));
   
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
% start commenting here for faster code run if visualisation is not
% required
if visualisation_flag==1   % If statement slows the code execution 
                           % nevertheless added it for pedagogical purposes  
hIm.CData = gather(EnTFSF(2:end-1,2:end-1));
drawnow                       
end
%--------------------------------------------------------------------------
% end commenting here for faster code run
%---------------------------- Update the fields ---------------------------
Enm1=En;            % Update the fields
En=Enp1;            % Update the fields
Esrc_nm1=Esrc_n;    % Update the fields
Esrc_n=Esrc_np1;    % Update the fields
end
%-------- Create a structure for storing and returning the field values ---
Efield_array_collection.Enm1=Enm1;
Efield_array_collection.En=En;
Efield_array_collection.Enp1=Enp1;
Efield_array_collection.Esrc_nm1=Esrc_nm1;
Efield_array_collection.Esrc_n=Esrc_n;
Efield_array_collection.Esrc_np1=Esrc_np1;
Efield_array_collection.EnTFSF=EnTFSF;
%--------------------------------------------------------------------------
end

