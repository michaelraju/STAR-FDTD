function source_data = define_dipole_beam_source(GPU_flag,init_data)
%-------------- Extract the data from the structures ----------------------
Ny=init_data.Ny;
kref=init_data.kref;
dt=init_data.dt;
dz=init_data.dz;
%----------------------------------
src_y=floor(Ny/2); % Placing the source along Y axis
beam_width=floor(Ny/2);  % size of the incident beam 
incident_angle=0*pi/180; % Incident angle in randians
source_data.beam_width=beam_width;
source_data.incident_angle=incident_angle;
src_yindex=zeros(Ny,1);
src_yindex(1+src_y-floor(beam_width/2) : src_y+floor(beam_width/2),1)=1;    
E0src=0.2;                                 % amplitude of the source function
%------------------ define source as a phased array at an angle ----------
wave_phase=zeros(Ny,2);
wave_phase(1+src_y-floor(beam_width/2) : src_y+ceil(beam_width/2),1)= ...
                          (kref*(0:beam_width-1).*dz).*tan(incident_angle);
wave_phase(1+src_y-floor(beam_width/2) : src_y+ceil(beam_width/2),2)=...
    wave_phase(1+src_y-floor(beam_width/2) : src_y+ceil(beam_width/2),1)+pi;
wave_amp(:,1)=src_yindex.*E0src;
wave_amp(:,2)=src_yindex.*E0src;
%----------------- Definition of the dipole source function -----------------
clip_lim=3;     % Colorbar clip limit for better contrast
slowly_rise =@(t_count) (1-exp(-((abs(t_count)-1)*dt)/(200*dt)));
                   % A function used to slowly grow the source, to
                   %       ease-in the source without an abrupt jump to the 
                   %       peak value.

if GPU_flag==1
   wave_amp=gpuArray(wave_amp);
   wave_phase=gpuArray(wave_phase);
end
%---------- Defining the dipole source as an inline function --------------
dipole_source = ...
 @(t_count,wave_amp,wave_phase,t_sign,fsrc) ...
 (slowly_rise(t_count).*wave_amp.*cos(t_sign*2*pi*fsrc*(t_count-1)*dt + ...
 wave_phase)); %t=+t : forward propagation, as t_sign=+1.
%------------- Return the parameters as a structure -----------------------
source_data.src_y=find(wave_amp(:,1));
source_data.wave_amp=wave_amp(source_data.src_y,:);
source_data.wave_phase=wave_phase(source_data.src_y,:);
source_data.dipole_source=dipole_source;
%--------------------------------------------------------------------------
end

