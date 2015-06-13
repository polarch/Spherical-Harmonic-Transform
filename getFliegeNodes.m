function [vecs, dirs, weights] = getFliegeNodes(index)
%GETFLIEGENODES Returns the spherical coordinates of Fliege-Maier nodes
%
%   GETFLIEGENODES returns the unit vectors, the spherical coordinates
%   and the weights of Fliege-MAier sets of points on the sphere for low-error
%   integration of spherical functions. The points can be used for integration 
%   through direct summation of the function evaluated at these points, and 
%   weighted with the respective weights. Each set is appropriate for spherical  
%   harmonic transform of order N = order+1, where order is the number of
%   the set. Up to N = 29 SHT is supported. The spherical coordinates are 
%   given in the [azi1 elev1; azi2 elev2; ...; aziQ elevQ] convention.
%
%   The designs have been copied from:
%       http://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html
%   and should be referenced as:
%       "A two-stage approach for computing cubature formulae for the sphere.",
%       Jorg Fliege and Ulrike Maier, Mathematik 139T, Universitat Dortmund, 
%       Fachbereich Mathematik, Universitat Dortmund, 44221. 1996.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, archontis.politis@aalto.fi, 7/7/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('fliegeMaierNodes_1_30.mat');

if index>30
    error('Designs of order greater than 30 are not implemented.')
elseif index<2
    error('Order should be at least 2.')
end

vecs = fliegeNodes{index}(:,1:3);
[dirs(:,1), dirs(:,2)] = cart2sph(vecs(:,1), vecs(:,2), vecs(:,3));

weights = fliegeNodes{index}(:,4);

end
