%------------ Forward FDTD wave propagation -------------------------------
tic
[Efield_array_collection] = FDTD_main_loop(GPU_flag,visualisation_flag,...
    init_data,source_data);
%------------ Phase sensitive detection for the forward propagation -------
init_data.detec_trans_position=...
 floor(((init_data.start_slab1 + init_data.scat_width) + init_data.Nz)/2);
                                % Transmission measurement z index location
  % Half way between the slab and the right boundary                              
init_data.detec_refl_position=floor(source_data.TFSF_Nz/2); 
                                % Reflection measurement z index location
% Half way between the TF/SF boundary and the left boundary                                                            
init_data.detec_src_position=init_data.start_slab1;
          % Source transmission measurement z index location
          % Placed at the slab surface incidence                        

PSD_array_collection= ...
  phase_sensitive_detection(GPU_flag,visualisation_flag,Efield_array_collection,...
                 init_data,source_data);
                   % Returns a data structure containing amplitude
                   % and phase of transmitted, reflected and source
                   % waves along with total transmission and reflection.
%--------------------------------------------------------------------------
sprintf('Total time taken for profiling the disorder without STM = %f min',toc/60)
%save data for the thesis/paper
if data_saving_flag==1
save('PSD_array_collection_diso_profiling.mat','PSD_array_collection')
save('Efield_array_collection_diso_profiling.mat','Efield_array_collection')
save('init_data_diso_profiling.mat','init_data')
%------------- save graphic images --------------------------------------
%save image for the thesis/paper
% figure(2)
% set(gca,'FontSize',12);
% exportgraphics(gcf,'forwardFDTD.eps','ContentType','vector')
end