function [Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO(UserVar,CtrlVar,MUA,GF,h,rhoi,rhow,varargin)

%% PICO melt rate parameterisation v0.91 ==================================
% An Ua implementation of the Potsdam Ice-shelf Cavity mOdel (PICO),
% details of the model can be found in Reese et al. 2018
% https://www.the-cryosphere.net/12/1969/2018/
%
%
% This is an Ua module intended to be called from within 'DefineMassBalance.m'
% The main output is a melt rate, given in metres per year
%
% NOTE: This has been changed to use the same sign convention as Ua, i.e.
% a negative number indicates basal melting and positive indicates freezing
%
% -------------------------------------------------------------------------
% To run PICO, use the following syntax:
%
% ab = PICO(UserVar,CtrlVar,MUA,GF,h,rhoi,rhow,PICO_opts)
%
% where:
%   UserVar: User Variables (Ua structure)
%   CtrlVar: Control Variables (Ua structure)
%   MUA: Mesh (Ua structure)
%   GF: Grounded/floating mask (Ua structure)
%   h: ice thickness (MUA.Nnodes x 1)
%   rhoi: ice density (scalar)
%   rhow: water density (scalar)
%   PICO_opts: optional structure containing various PICO options (details below)
%
%% PICO_opts ==============================================================
%
% This is a structure similar to CtrlVar that contains all the options
% available when running the PICO model. It is not strictly necessary to
% define this structure when calling PICO but it is highly recommended,
% particularly for problem specific fields e.g. related to the ambient
% ocean fields. In particular, think carefully about your choices for the
% basins interpolant, temperature/salinity fields, C1 and gamTstar.
% -------------------------------------------------------------------------
% General options:
%
% - PICO_opts.algorithm: 'graph', 'watershed', 'polygon' or 'oneshelf' (DEFAULT = 'watershed')
% There are options related to how ice shelves are delineated. The
% 'watershed' options converts MUA into a structured grid and then uses
% image processing techniques to quickly define connected floating regions
% as individual ice shelves. The graph option constructs a graph network 
%from the triangulation to find floating nodes that are connected in order
% to seperate each ice shelf. The 'polygon' option creates polygons out of
% each individual GL segment and the MeshBoundaryCoordinates and defines
% individual ice shelves as floating nodes within these polygons. In
% general the 'watershed' option will be considerably faster while the
% 'polygon' option is much slower but generally more robust. The 'oneshelf'
% option should only be used when you are certain that you will only ever
% have one distinct ice shelf in your domain and this will make the code
% considerably faster. Note: In this case any checks on shelf area etc. are
% ignored.
%
% - PICO_opts.SmallShelfMelt (DEFAULT = 0)
% The melt rate applied to floating nodes that for whatever reason are not
% within an ice shelf as determined by this code.
%
% - PICO_opts.FloatingCriteria: 'GLthreshold' or 'StrictDownstream' (DEFAULT = 'StrictDownstream')
% Melt is only applied at nodes that meet this floating criteria which is
% either determined by CtrlVar.GLthreshold or nodes strictly downstream of
% the grounding line.
%
% - PICO_opts.nmax (DEFAULT = 5)
% Maximum number of boxes possible within an ice shelf.
% 
% - PICO_opts.minArea (DEFAULT = 2e9)
% Minimum area (m^2) for a floating region to be considered an ice shelf.
%
% - PICO_opts.minNumShelf (DEFAULT = 20)
% Minimum number of floating nodes (polygon + graph) or grid cells 
% (watershed) for a given region to be considered an ice shelf.
%
% - PICO_opts.InfoLevel (DEFAULT = 100)
% Same functionality as CtrlVar.InfoLevel, choose between 0,1,10,100 with 
% higher InfoLevels providing more information during the PICO run.
%
% -------------------------------------------------------------------------
% Options related to the ambient ocean properties: 
%
% - PICO_opts.BasinsInterpolant: Where the user wants to have multiple ocean
% basins with different properties, this should be a scatteredInterpolant 
% covering the entire domain and numbered from 1 to n, where n is the maximum 
% number of ocean basins. See Reese et al (2018) for an example of how to 
% define the ocean basins. If this is empty the model will assume only one
% basin in the entire domain (and hence use only one T0 and S0).
% 
% - PICO_opts.Tbasins: A vector of length n, where each element (i) of
% Tbasins is the ambient ocean temperature of basin number (i). If no
% BasinsFile is given Tbasins should be a scalar giving the ambient ocean
% temperature of the entire domain.
%
% - PICO_opts.Sbasins: A vector of length n, where each element (i) of
% Sbasins is the ambient ocean salinity of basin number (i). If no
% BasinsFile is given Sbasins should be a scalar giving the ambient ocean
% salinity of the entire domain.
%
% -------------------------------------------------------------------------
% Options specific to the melt rate model:
%
% - PICO_opts.C1: Overturning strength (DEFAULT = 1e6)
% - PICO_opts.gamTstar: turbulent temp. exch. coeff. (DEFAULT = 2e-5)
% 
% -------------------------------------------------------------------------
% Options related to PICO_opts.algorithm = 'watershed'
%
% - PICO_opts.PICOres (DEFAULT = 6e3/3e3/1e3 depending on mesh size)
% The resolution of the structured grid that the watershed option uses,
% this can potentially slow things down a lot and a higher number (lower
% resolution) typically seems to work fine.
%
% - PICO_opts.ContinentArea (DEFAULT = 5e10)
% This helps to identify the main grounding line of interest by setting an
% area threshold below which grounded areas are treated as islands and
% hence these grounding lines are ignored for the purposes of defining the
% boxes. If you want to include all grounding lines then set this to 0 but
% think carefully about your decision and check how this parameter affects
% your box geometry. The default value is set so that a simulation for the
% whole of Antarctica will include Alexander Island but ignore Berkner
% island to replicate the behaviour of the PICO implementation in PISM.
% Note that the graph algorithm also makes use of this number.
% 
% -------------------------------------------------------------------------
% Options specific to PICO_opts.algorithm = 'polygon'
%
% - PICO_opts.MeshBoundaryCoordinates
% This should be the same n x 2 matrix that is provided to
% Ua2D_InitialUserInput.m
%
% - PICO_opts.persistentBC (DEFAULT = 0)
% To avoid reordering mesh boundary coordinates every time PICO is called,
% set this parameter to 1

if nargin < 7
    error('Some PICO inputs appear to be undefined');
elseif nargin<8
    warning('PICO_opts undefined, using only default values... ARE YOU SURE YOU WANT TO DO THIS?');
    PICO_opts = struct;
else
    PICO_opts = varargin{1};
end

[Mk,ShelfID,T0,S0,Tkm,Skm,q,PBOX,Ak] = PICO_driver(UserVar,CtrlVar,MUA,GF,h,rhoi,rhow,PICO_opts);