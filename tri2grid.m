function [Xi,Yi,Zi,in] = tri2grid(MUA,Fin,res)

% This function takes the MUA structure and a field (Fin) and linearly maps
% that field to a structured grid defined by the extent of MUA and a user
% defined uniform grid resolution (res)
%
% Usage: [Xi,Yi,Zi,in,out] = tri2grid(MUA,Fin,res)
%
% Inputs: 
% - MUA: the Ua Mesh
% - Fin: the field to map (with dimensions MUA.Nnodes x 1
% - res: a scalar defining the uniform grid resolution
%
% Outputs:
% - Xi, Yi: the uniform grid onto which Fin is mapped
% - Zi: the variable Fin mapped onto the grid Xi, Yi
% - in: logical vector that is true if the grid square is within the mesh
%       NB: this is not the same as within the convex hull which is what
%       scatteredInterpolant would use

x = MUA.coordinates(:,1);
y = MUA.coordinates(:,2);

% xg and yg are the coordinates of the corners of each pixel in the image
xg = min(x):res:max(x)+res; % need to go beyond the domain for reasons
yg = min(y):res:max(y)+res;
[Xi,Yi] = meshgrid(xg,yg);

xx = reshape(Xi,[],1);
yy = reshape(Yi,[],1);

[ti,bc] = pointLocation(MUA.TR,[xx yy]);

in = ~isnan(ti);

ti2 = ti(in);
bc2 = bc(in,:);

Vals = Fin(MUA.connectivity);

triVals = Vals(ti2,:);

Vq = dot(bc2',triVals')';

Zi = xx*0;
Zi(in) = Vq;
Zi = reshape(Zi,size(Xi));
in = reshape(in,size(Xi));
