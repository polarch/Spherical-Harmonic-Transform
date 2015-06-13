function h_ax = plotSphFunctionGrid(F, aziRes, polarRes, realComplex, h_ax)
%PLOTSPHFUNCTIONGRID Plots a spherical function defined on a grid
%
%   F:  matrix of function values on the grid points
%   aziRes: grid resolution at azimuth (degrees)
%   polarRes: grid resolution in elevation (degrees)
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

if nargin==4
    figure
    h_ax = axes;
elseif nargin==3
    figure
    h_ax = axes;   
    realComplex = 'complex';
else
    axes(h_ax);
end

% construct grid
azi = 0:aziRes:360;
elev = 0:polarRes:180;
azi_rad = azi*pi/180;
elev_rad = elev*pi/180;
[Az, El] = meshgrid(azi_rad, elev_rad);

% construct real positive and negative parts if real function
D_x = cos(Az).*sin(El).*squeeze(abs(F));
D_y = sin(Az).*sin(El).*squeeze(abs(F));
D_z = cos(El).*squeeze(abs(F));

if isequal(realComplex, 'real')
    Dp_x = D_x.*(F>=0);
    Dp_y = D_y.*(F>=0);
    Dp_z = D_z.*(F>=0);
    Dn_x = D_x.*(F<0);
    Dn_y = D_y.*(F<0);
    Dn_z = D_z.*(F<0);
elseif isequal(realComplex, 'complex')
    Dm_x = D_x.* abs(F);
    Dm_y = D_y.* abs(F);
    Dm_z = D_z.* abs(F);
end

% plot 3d axes
maxF = max(max(abs(F)));
line([0 1.5*maxF],[0 0],[0 0],'color',[1 0 0])
line([0 0],[0 1.5*maxF],[0 0],'color',[0 1 0])
line([0 0],[0 0],[0 1.5*maxF],'color',[0 0 1])

% plot function
hold on
if isequal(realComplex, 'real')
    Hp = surf(Dp_x, Dp_y, Dp_z);
    Hn = surf(Dn_x, Dn_y, Dn_z);
    set(Hp, 'FaceColor', 'b')
    set(Hp, 'EdgeAlpha', 1)
    set(Hn, 'FaceColor', 'r')
    set(Hn, 'EdgeAlpha', 1)
elseif isequal(realComplex, 'complex')
    Hm = surf(Dm_x, Dm_y, Dm_z, angle(F));
    set(Hm, 'EdgeAlpha', 1)
end
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
