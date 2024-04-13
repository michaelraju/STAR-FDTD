function source_data = define_dipole_point_source(GPU_flag,init_data)
%-------------- Extract the data from the structures ----------------------
Ny=init_data.Ny;
src_y=floor((Ny+2)/2);
fsrc=init_data.fsrc;
dt=init_data.dt;
%----------------- Definition of the dipole source function -----------------
wave_amp=[1 1];    % Dipole amplitude
wave_phase=[0 pi]; % Dipole phase
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
 @(t_count,wave_amp,wave_phase,t_sign) (slowly_rise(t_count).*wave_amp.*cos(t_sign*2*pi*fsrc*(t_count-1)*dt + ...
 wave_phase)); %t=+t : forward propagation, as t_sign=+1.
%------------- Return the parameters as a structure -----------------------
source_data.wave_amp=wave_amp;
source_data.wave_phase=wave_phase;
source_data.src_y=src_y;
source_data.dipole_source=dipole_source;
%--------------------------------------------------------------------------
end

