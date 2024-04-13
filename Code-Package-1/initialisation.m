function [init_data] = initialisation(GPU_flag,NT_steady_state,krefL,krefW,...
    nref,sigma,z_domain_size_factor,no_of_PSD_cycles)

c0=299792458;            % Velocity of the wave in free space
fsrc=3.2*(10^14);        % Temporal frequency of the source in Hz
lambda0=c0/fsrc;         % Wavelength in nm
k0=2*pi/lambda0;         % Free space wave vector
kref=2*pi*nref/lambda0;  % Background(reference) medium wave vector 
lambda_ref=(2*pi)/kref;  % Background(reference) medium wavelength
eps_avg=nref^2;          % Background(reference) medium dielectric constant
%------------------- Initilisation : grid parameters ----------------------
eta_min=sqrt(eps_avg-(sigma*(eps_avg-1))); % Minimum refractive index
eta_max=sqrt(eps_avg+(sigma*(eps_avg-1))); % Maximum refractive index
lam_fact=20;                 % A constant used for setting up grid size. 
                             % Take lam_fact >= 15 
                             % Larger the lam_fact value, finer the grid
                             % and larger memory needed
dz=(c0/(fsrc*eta_max))/lam_fact;
                             % dz is (1/lam_fact) of the minimum wavelength
dy=dz;                       %   implying a square grid cell. 
                             % Always keep dz=dy as it is an 
                             %       assumption taken for simplifying  
                             %       the derivations of analytical results.
   
time_fact=2;     % A constant >1 used for setting up the time step.
                 % Take 2 as a conservative value for safe time stepping 
                 %   as given below
dt=(dz/(sqrt(2)*c0/eta_min))/time_fact; % dt is (1/time_fact) of the time 
                 %  taken to physically cross a grid cell size
                 %  This is to satisy the Courant–Friedrichs–Lewy condition                               
%----------------- Define geometry of the problem -------------------------
mirror_pos=1;           % Define transverse mirror boundary position offset 
Ny=ceil(krefW/(kref*dy))+ 2;  % Transverse dimensionaless width kref*W
                              % Added 2 for the boundary mirrors
macropixel_size=8; % Particle macro pixel size in terms of no of grid 
                    %  points. This is to increase scattering efficiency.                            
Ny=ceil(Ny/macropixel_size)*macropixel_size; % Round off to the nearest 
                                      %  macropixel size for grid snapping.    
krefW=(Ny-1)*kref*dy; %Reevaluate the dimensionless transverse length 
scat_width=ceil(krefL/(kref*dz));   % Width of the scattering slab in terms 
                                    %   of no of grid points.
scat_width=ceil(scat_width/macropixel_size)*macropixel_size;  
                              % Round off to the nearest macropixel size  
krefL=(scat_width-1)*kref*dz; % Reevaluate the dimensionless thickness                                                  
Nz=floor(z_domain_size_factor*scat_width);   
                               % Longitudal dimensionless thickness kref*L
start_slab1=floor((Nz/2)-(scat_width/2)); 
                         % Starting z position of the random slab
epsilon_ij= ...          % Initialising epsilon_ij, the dielectric constant
 eps_avg.*ones(Ny+2,Nz+2); 
epsilon_ij(2:end-1,start_slab1:start_slab1+scat_width-1)= ...
   macropixels(scat_width,Ny,macropixel_size,macropixel_size,eps_avg,sigma);
              % The function macropixel returns a random grid of dielectric
              %                       constant as a function of z and y.
%------------------ Plot the disorder refractive index --------------------
figure('position',[100 100 600 600])
colormap jet
imagesc([0 kref*dz*(Nz-1)],[0 kref*dy*(Ny-1)],sqrt(epsilon_ij));  
% Refractive index = sqrt(dielectric constant)
title('$Refractive~index~(\eta_{z,y})~distribution$','Interpreter','Latex')
xlabel('$k_{ref}z$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis xy equal tight
colorbar
set(gca,'FontSize',22)
%------------ Some Checks for better performance --------------------------

% Courant–Friedrichs–Lewy (CFL) condition for stable operation of FDTD
if ((c0/eta_min)*dt/dz) > 1/sqrt(2)
    disp('Warning!! CFL condition not satisfied')
    pause;
else 
    sprintf('Good that the CFL constant = %f is less than 1/sqrt(2)=0.7071',((c0/eta_min)*dt/dz))
end

if (start_slab1 + scat_width) > Nz
    disp('Error : Sample size exceeds the grid boundaries')
    pause;
end

if (dy~=dz)
 disp('dz has to be equal to dy for the validity of the analytical derivations')
 pause;
end  

if eta_min < 1
   sprintf('Warning : eta_min < 1 !!!')  % Realistically refractive index
                                         %    preferred to be larger than 1
   pause;
end

if kref*dy > 0.3
sprintf('k_{ref}dy = %f',kref*dy)   % Recommended kref*dz=0.3 or less
sprintf('k_{ref}dz = %f',kref*dz)   % Recommended kref*dy=0.3 or less
disp('kref dy = kref dz is expected to be < 0.3')
pause;
end

% ------------- initialise GPU arrays -------------------------------------
if GPU_flag==1
w0=gpuArray(((c0*dt/dy)^2)./(epsilon_ij(2:end-1,2:end-1))); 
                   % FDTD factor to be used in the field update equations.
wsrc0=gpuArray(((c0*dt/dz)^2)./((eps_avg))); 
                   % FDTD factor to be used in the source only field 
                   %                update equations.
else 
% Defining the same matrices, but without using gpuArray
w0=(((c0*dt/dy)^2)./(epsilon_ij(2:end-1,2:end-1)));
wsrc0=(((c0*dt/dz)^2)./((eps_avg)));
end

%-- Create an initialisation structure to pass on to other functions ------
init_data.Ny=Ny;
init_data.Nz=Nz;
init_data.kref=kref;
init_data.dz=dz;
init_data.dy=dy;
init_data.dt=dt;
init_data.start_slab1=start_slab1;
init_data.scat_width=scat_width;
init_data.mirror_pos=mirror_pos;
init_data.NT_steady_state=NT_steady_state;
init_data.w0=w0;
init_data.wsrc0=wsrc0;
init_data.fsrc=fsrc;
init_data.c0=c0;
init_data.eps_avg=eps_avg;
init_data.epsilon_ij=epsilon_ij;
init_data.nref=nref;
init_data.lambda0=lambda0;
init_data.sigma=sigma;
init_data.lam_fact=lam_fact;
init_data.time_fact=time_fact;
init_data.start_slab1=start_slab1;
init_data.lambda_ref=lambda_ref;
init_data.k0=k0;
init_data.krefW=krefW;
init_data.krefL=krefL;
init_data.z_domain_size_factor=z_domain_size_factor;
init_data.no_of_PSD_cycles=no_of_PSD_cycles;
end

