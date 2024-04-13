tic
[Efield_array_collection] = FDTD_main_loop_STM(GPU_flag,visualisation_flag,...
                            init_data,source_data,STM_data);

%---- Phase sensitive detection for the harmonic (sideband) wavefronts ----
init_data.detec_trans_position=...
 floor((init_data.start_slab2 + init_data.scat_width + init_data.Nz)/2);
                                % Transmission measurement z index location
  % Half way between the slab and the right boundary                              
init_data.detec_refl_position=floor(source_data.TFSF_Nz/2); 
                                % Reflection measurement z index location
% Half way between the TF/SF boundary and the left boundary                                                            
init_data.detec_src_position=init_data.start_slab1;
          % Source transmission measurement z index location
          % Placed at the slab surface incidence                         

norder_array=[+1 0]; % Detect the upper side band and the center frequency
                     % out of which only the upper side band is needed for
                     % phase conjugation.
for ncount=1:length(norder_array)
nth=norder_array(ncount);
PSD_array_collection=PSD_nth_SideBand(GPU_flag,visualisation_flag,...
    Efield_array_collection,init_data,source_data,STM_data,nth);
title(sprintf('$Phase~sensitive~detection~: \\bf n=%d$',nth),'Interpreter','Latex')
if ncount<length(norder_array)
close; % Close some of the figures to prevent accumulation of lot of 
       % images for the entire IOPC
end
PSD_structure_collection(ncount)=PSD_array_collection;
end
%--------------------------------------------------------------------------
%save data for the thesis
save('PSD_array_collection_STM.mat','PSD_array_collection')
save('Efield_array_collection_STM.mat','Efield_array_collection')
save('STM_data.mat','STM_data')
save('init_data_forward.mat','init_data')

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
    krefy_detec,'-*',...
    'MarkerSize',1)
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
sprintf('Total time taken for forward propagation with STM = %f min',toc/60)
pause(1);
saveas(gcf,'STM_forward_L2R.fig')
