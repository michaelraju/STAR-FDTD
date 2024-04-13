% Generates random slab disorder (dielectric constant) with a given macropixel size
function [eps_profile] = macropixels(Nz,Ny,pNz,pNy,eps_avg,sigma)
rng shuffle;
if mod(Nz,pNz)~=0
sprintf('Error in macropixel dividing for snapping. Choose differently !!')    
                                        %--- Need to be exactly dividable
pause;
end
Mz=Nz/pNz;
My=Ny/pNy;
eps_M=eps_avg+(sigma*(eps_avg-1)).*(2.*rand(My,Mz)-1); %---- Macro-pixels
for jcount=1:My
   for icount=1:Mz
   eps_profile(pNy*(jcount-1)+1:pNy*jcount,pNz*(icount-1)+1:pNz*icount)= ...
                                                      eps_M(jcount,icount);
   end
end
end

