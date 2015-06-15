function [aziElev] = unitCart2sph(xyz)
%UNITCART2SPH Get directly azimuth and elevation from unit vectors
%   Similar to cart2sph, assuming unit vectors and using matrices for the 
%   inputs and outputs, instead of separate vectors of coordinates. More
%   convenient for many of the SH functions found here.

if size(xyz,2)~=3, xyz = xyz.'; end
    
[aziElev(:,1), aziElev(:,2)] = cart2sph(xyz(:,1), xyz(:,2), xyz(:,3));

end

