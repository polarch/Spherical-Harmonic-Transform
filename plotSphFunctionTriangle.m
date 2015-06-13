function h_ax = plotSphFunctionTriangle(F, dirs, realComplex, h_ax)
%PLOTSPHFUNCTIONTRIANGLE Plots a spherical function on unstructured grid
%
%   F:  vector of K function values on the sampling points
%   dirs:   [azimuth1 inclination1; ...; azimuthK inclinationK] angles in 
%           rads for each evaluation point, where inclination is the polar 
%           angle from zenith: inclination = pi/2-elevation
%   realComplex: {'real','complex'} if the function is real then it is
%                plotted with one surface for its positive part and one for
%                its negative part. If it is complex, the magnitude
%                function is plotted, with its phase mapped on the colormap
%   h_ax: optional argument to define an axis handle for the plot,
%         otherwise the new axis handle is returned
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 20/02/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==3
    figure
    h_ax = axes;
elseif nargin==2
    figure
    h_ax = axes;   
    realComplex = 'complex';
else
    axes(h_ax);
end

% triangulate sampling points
aziIncl2aziElev = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];
dirs = aziIncl2aziElev(dirs);
[triangulated.vertices(:,1), triangulated.vertices(:,2), triangulated.vertices(:,3)] = ...
    sph2cart(dirs(:,1), dirs(:,2), 1);
triangulated.vertices = abs(F)*ones(1,3) .* triangulated.vertices;
triangulated.faces = sphDelaunay(dirs);

% construct real positive and negative colormap if real function
CData = zeros(length(F),3);
if isequal(realComplex, 'real')
    pos_idx = find(F>=0);
    neg_idx = find(F<0);
    for i=1:length(pos_idx)
        CData(pos_idx(i),:) = [0 0 255];
    end
    for i=1:length(neg_idx)
        CData(neg_idx(i),:) = [255 0 0];
    end
elseif isequal(realComplex, 'complex')
    CData = angle(F);
end

% plot 3d axes
maxF = max(max(abs(F)));
line([0 1.5*maxF],[0 0],[0 0],'color',[1 0 0])
line([0 0],[0 1.5*maxF],[0 0],'color',[0 1 0])
line([0 0],[0 0],[0 1.5*maxF],'color',[0 0 1])

% plot function
hold on
triangulated.Facecolor = 'interp';
triangulated.Edgecolor = 'k';
triangulated.FaceVertexCData = CData;
patch(triangulated)
xlabel('x')
ylabel('y')
zlabel('z')
light('Position',[0 0 1],'Style','infinite');
light('Position',[-1 -1 -1],'Style','infinite');
light('Position',[0 0 -1],'Style','infinite');
material shiny
axis equal
grid

end
