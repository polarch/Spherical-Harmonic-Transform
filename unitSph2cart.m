function [xyz] = unitSph2cart(aziElev)
%UNITSPH2CART Get coordinates of unit vectors from azimuth-elevation
%   Similar to sph2cart, assuming unit vectors and using matrices for the 
%   inputs and outputs, instead of separate vectors of coordinates. More
%   convenient for many of the SH functions found here.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 6/5/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(aziElev,2)~=2, aziElev = aziElev.'; end
    
[xyz(:,1), xyz(:,2), xyz(:,3)] = sph2cart(aziElev(:,1), aziElev(:,2), 1);

end

