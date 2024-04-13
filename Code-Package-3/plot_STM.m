%----------------- To mark the position of the STM ------------------------
if  source_data.t_sign==-1
    STM_plot_color='white';
else
    STM_plot_color='blue';
end
    
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
    'color',STM_plot_color,'LineWidth',2)
%--------------------------------------------------------------------------