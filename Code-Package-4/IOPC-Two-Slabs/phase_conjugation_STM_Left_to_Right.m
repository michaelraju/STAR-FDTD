nchosen=nchosen_L2R;
PSD_array_collection=PSD_structure_collection(find(norder_array_redefined==nchosen));
init_pc_data=init_data;  % Initialisation data for phase conjugation
init_pc_data.NT_steady_state=NT_steady_state_PC;
init_pc_data.no_of_PSD_cycles=no_of_PSD_cycles_STM;
init_pc_data.fsrc=init_data.fsrc+nchosen*STM_data.f_mod;  
                                     % modulated frequency taken 
%---------------- phase conjugation : Right to Left -----------------------
source_data.t_sign=-1;   
%src_z_start=floor(init_data.start_slab1/2);  % Half way between slab and the
                                            % left boundary
source_data.src_z=[src_z_start src_z_start+1];     % Dipole source location
source_data.TFSF_Nz=floor(source_data.src_z(1)/2); % TFSF boundary location 
                        % Half way between the source and the left boundary 
source_data.TFSF_region=1:source_data.TFSF_Nz;

source_data.wave_amp= ...
    repmat(PSD_array_collection.E_trans_sig(2:end-1),1,2);% dipole source amplitude
source_data.wave_amp= ...
     0.1.*source_data.wave_amp./max(source_data.wave_amp); 
                                         % normalise the amplitude                                                                                                                                         
source_data.wave_phase= ...
    repmat(PSD_array_collection.E_theta_trans_sig(2:end-1),1,2); 
source_data.wave_phase(:,2)=source_data.wave_phase(:,2)+pi; % dipole phase
source_data.clip_lim=gather(max(source_data.wave_amp(:)))+0.6;



% [Efield_array_collection] = ...
%     FDTD_main_loop(GPU_flag,visualisation_flag,init_pc_data,source_data);
[Efield_array_collection] = FDTD_main_loop_STM(GPU_flag,visualisation_flag,...
    init_pc_data,source_data,STM_data);

%------------ Phase sensitive detection -----------------------------------
init_pc_data.detec_trans_position=...
 floor((init_data.start_slab2 + init_data.scat_width + init_data.Nz)/2);
                                % Transmission measurement z index location
init_pc_data.detec_refl_position=floor(source_data.TFSF_Nz/2); 
                                % Reflection measurement z index location
% init_pc_data.detec_src_position= ...
%   init_pc_data.detec_trans_position;   %Source transmission detection location
init_pc_data.detec_src_position=init_pc_data.start_slab1; 
   % Source transmission detection location
   % Set as the LHS incident surface of the slab


norder_array=[1 0];                      
norder_array_redefined=norder_array+nchosen; % The effective array notation values.
% One has to take care of the term norder_array, because of the way,
%                FDTD functions are defined, where norder_array is defined
%                 relative to the frequency of the source used.

for ncount=1:length(norder_array)
nth=norder_array(ncount);
PSD_array_collection=PSD_nth_SideBand(GPU_flag,visualisation_flag,...
    Efield_array_collection,init_pc_data,source_data,STM_data,nth);
if ncount < length(norder_array)
close; % Close some of the figures to prevent accumulation of lot of 
       % images for the entire IOPC
end
PSD_structure_collection(ncount)=PSD_array_collection;
end
%-------------------------- generated the focus ---------------------------
sprintf('Total time taken for phase conjugation = %f min',toc/60)

filename1=sprintf('PSD_array_collection_PC_L2R_%d.mat',r_count);
filename2=sprintf('Efield_array_collection_PC_L2R_%d.mat',r_count);

%save data for the thesis/paper
if data_saving_flag==1
save(filename1,'PSD_array_collection')
save(filename2,'Efield_array_collection')
% figure(6)
% xlabel('$k_{ref}z$','Interpreter','Latex')
% ylabel('$k_{ref}y$','Interpreter','Latex')
% set(gca,'FontSize',12);
% exportgraphics(gcf,'PCFDTD.eps','ContentType','vector')
end 

%----------------- To mark the position of the STM ------------------------
hold on
STM_start_z=STM_data.STM_start_z;
STM_start_y=STM_data.STM_start_y;
grat_width=STM_data.grat_width;
STM_end_y=STM_data.STM_end_y;
kref=init_data.kref;
dz=init_data.dz;
dy=init_data.dy;
line([STM_start_z*kref*dz (STM_start_z+grat_width-1)*kref*dz ...
    (STM_start_z+grat_width-1)*kref*dz ...
    STM_start_z*kref*dz  STM_start_z*kref*dz ], ...
    [STM_start_y*kref*dy STM_start_y*kref*dy STM_end_y*kref*dy ...
    STM_end_y*kref*dy  STM_start_y*kref*dy], ...
    'color','blue','LineWidth',2)

side_description={sprintf('$Phase~conjugation$'),...
                  sprintf('$from~the~left~to~the~right$'),...
                  sprintf('$Iteration~no~: %d$',r_count)};
annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
   side_description,'FontSize',16, ...
   'Interpreter','latex','FitBoxToText','on'); 
saveas(gcf,sprintf('STM_PC_L2R_%d.fig',r_count))
%--------------------------------------------------------------------------
%--------- plot the extracted sideband wavefronts -------------------------
krefy_detec=(0:init_data.Ny+1).*(init_data.kref*init_data.dy);

% peak_val_store=zeros(1,length(norder_array));
% for ncount=1:length(norder_array)
% peak_val_store(ncount)=max(PSD_structure_collection(ncount).E_trans_sig);
% end

% figure('position',[70 70 500 800],'color','white')
% xlabel('Transmitted~intensity~(dB)','Interpreter','Latex')
% ylabel('$k_{ref}y$','Interpreter','Latex')
% axis tight
% set(gca,'FontSize',24)
% for ncount=1:length(norder_array)
% nth=norder_array(ncount);
% plot(mag2db(PSD_structure_collection(ncount).E_trans_sig./max(peak_val_store)),krefy_detec,'-*',...
%     'MarkerSize',1)
% title(sprintf('$Temporal~order~no,~n = %d$',norder_array_redefined(ncount)),...
%     'Interpreter','Latex')
% xlabel('$Transmitted~intensity~(dB)$','Interpreter','Latex')
% ylabel('$k_{ref}y$','Interpreter','Latex')
% ylim([0 init_data.kref*init_data.dy*(init_data.Ny-1)])
% set(gca,'FontSize',22)
% hold on
% pause(0.5)
% end
% 
% for lcount=1:numel(norder_array)
%     if norder_array(lcount)>0
%     legend_array{lcount}=sprintf('$n=+%d$',norder_array_redefined(lcount));
%     elseif norder_array(lcount)<=0
%     legend_array{lcount}=sprintf('$n=%d$',norder_array_redefined(lcount));
%     end
% end
% legend(legend_array,'Interpreter','Latex');
% title('');
% 
% side_description={sprintf('$Phase~conjugation$'),...
%                   sprintf('$from~the~left~to~the~right$'),...
%                   sprintf('$Iteration~no~: %d$',r_count)};
% annotation('textbox', [0.005, 0.95, 0.001, 0.001], 'string', ...
%    side_description,'FontSize',16, ...
%    'Interpreter','latex','FitBoxToText','on'); 
% saveas(gcf,sprintf('STM_PSD_L2R_%d.fig',r_count))

tic
%---------- Intensity profile at the STM Y plane -------------------------
figure('position',[70 70 500 800],'color','white') 
field_magnitude_at_STM=abs(Efield_array_collection.EnTFSF(:,STM_data.STM_start_z));
plot(mag2db(field_magnitude_at_STM./max(field_magnitude_at_STM)),...
    (0:init_data.Ny+1).*init_data.kref*init_data.dy,'-*','MarkerSize',1)
set(gca,'xdir','reverse')
xlabel('$dB$','Interpreter','Latex')
title('$Normalised~Intensity~at~the~STM~plane$','Interpreter','Latex')
ylabel('$k_{ref}y$','Interpreter','Latex')
axis tight
set(gca,'FontSize',22)

saveas(gcf,sprintf('STM_focus_profile_L2R_%d.fig',r_count))
