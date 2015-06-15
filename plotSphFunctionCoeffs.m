function h_ax = plotSphFunctionCoeffs(F_N, basisType, aziRes, polarRes, realComplex, h_ax)
%PLOTSPHFUNCTIONGRID Plots a spherical function defined on a grid
%
%   F_N:  (N+1)^2 vector of SH coefficients
%   basisType:  {'real','complex'} SH basis type for the coefficients
%   aziRes: grid resolution at azimuth (degrees)
%   polarRes: grid resolution in elevation (degrees)
%   realComplex: {'real','complex'} if the function is real then it is
%                plotted with blue for its positive part and red for
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

if nargin==5
    figure
    h_ax = axes;
elseif nargin==4
    figure
    h_ax = axes;
    realComplex = 'complex';   
elseif nargin==2
    figure
    h_ax = axes;
    realComplex = 'complex';
    aziRes = 5;
    polarRes = 5;
else
    axes(h_ax);
end

% get function values at grid by the inverse SHT
dirs = grid2dirs(aziRes, polarRes);
F = inverseSHT(F_N, dirs, basisType);
Fgrid = Fdirs2grid(F, aziRes, polarRes, 1);

% construct grid
azi = 0:aziRes:360;
elev = 0:polarRes:180;
azi_rad = azi*pi/180;
elev_rad = elev*pi/180;
[Az, El] = meshgrid(azi_rad, elev_rad);

% construct real positive and negative parts if real function
D_x = cos(Az).*sin(El).*squeeze(abs(Fgrid));
D_y = sin(Az).*sin(El).*squeeze(abs(Fgrid));
D_z = cos(El).*squeeze(abs(Fgrid));

if isequal(realComplex, 'real')
    Dp_x = D_x.*(Fgrid>=0);
    Dp_y = D_y.*(Fgrid>=0);
    Dp_z = D_z.*(Fgrid>=0);
    Dn_x = D_x.*(Fgrid<0);
    Dn_y = D_y.*(Fgrid<0);
    Dn_z = D_z.*(Fgrid<0);
elseif isequal(realComplex, 'complex')
    Dm_x = D_x;
    Dm_y = D_y;
    Dm_z = D_z;
end

% plot 3d axes
maxF = max(max(abs(Fgrid)));
line([0 1.1*maxF],[0 0],[0 0],'color',[1 0 0])
line([0 0],[0 1.1*maxF],[0 0],'color',[0 1 0])
line([0 0],[0 0],[0 1.1*maxF],'color',[0 0 1])

% plot function
hold on
if isequal(realComplex, 'real')
    Hp = surf(Dp_x, Dp_y, Dp_z);
    Hn = surf(Dn_x, Dn_y, Dn_z);
    set(Hp, 'FaceColor', 'b')
    set(Hp, 'EdgeAlpha', 0.2)
    set(Hn, 'FaceColor', 'r')
    set(Hn, 'EdgeAlpha', 0.2)
elseif isequal(realComplex, 'complex')
    Hm = surf(Dm_x, Dm_y, Dm_z, angle(Fgrid));
    set(Hm, 'EdgeAlpha', 0.2)
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
