function cmap = wavecolormap
                                                        % Code developed by
                                                         % Michael Raju
no_of_interp_colours=60;
colours  = [0.65   1  1;     % brighter cyan
             0    1  1;     % cyan
             0    0  1;     % blue
             0    0  0;     % black
             1    0  0;     % red
             1    1  0;     % yellow
             1    1  0.65];  % brighter yellow
no_of_base_colours=size(colours,1);       
step_per_color=ceil(no_of_interp_colours/no_of_base_colours);
cmap=zeros(step_per_color*(no_of_base_colours-1),3);
% Next interpolate between colours to get a smoother profile       
% The following code is very specific to the list of 7 colors.

% interpolate colour 1
icount=1;
index=1+(icount-1)*step_per_color:step_per_color*icount;    
cmap(index,1) = linspace(colours(icount,1),colours(icount+1,1),step_per_color) ;
cmap(index,2:3)=repmat(colours(icount,2:3),step_per_color,1); 
% interpolate colour 2
icount=2;
index=(icount-1)*step_per_color:step_per_color*icount;    
cmap(index,2) = linspace(colours(icount,2),colours(icount+1,2),step_per_color+1) ;
cmap(index,[1 3])=repmat(colours(icount,[1 3]),step_per_color+1,1); 
% interpolate colour 3
icount=3;
index=(icount-1)*step_per_color:step_per_color*icount;    
cmap(index,3) = linspace(colours(icount,3),colours(icount+1,3),step_per_color+1) ;
cmap(index,[1 2])=repmat(colours(icount,[1 2]),step_per_color+1,1); 
% interpolate colour 4
icount=4;
index=(icount-1)*step_per_color:step_per_color*icount;    
cmap(index,1) = linspace(colours(icount,1),colours(icount+1,1),step_per_color+1) ;
cmap(index,2:3)=repmat(colours(icount,2:3),step_per_color+1,1); 
% interpolate colour 5
icount=5;
index=(icount-1)*step_per_color:step_per_color*icount;    
cmap(index,2) = linspace(colours(icount,2),colours(icount+1,2),step_per_color+1) ;
cmap(index,[1 3])=repmat(colours(icount,[1 3]),step_per_color+1,1); 
% interpolate colour 6
icount=6;
index=(icount-1)*step_per_color:step_per_color*icount;    
cmap(index,3) = linspace(colours(icount,3),colours(icount+1,3),step_per_color+1) ;
cmap(index,[1 2])=repmat(colours(icount,[1 2]),step_per_color+1,1); 
end



