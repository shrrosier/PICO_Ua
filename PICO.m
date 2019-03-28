
%% PICO melt rate parameterisation v0.5 ===================================
% An Ua implementation of the Potsdam Ice-shelf Cavity mOdel (PICO),
% details of the model can be found in Reese et al. 2018
% https://www.the-cryosphere.net/12/1969/2018/
%
%
% This is an Ua module intended to be called from within 'DefineMassBalance.m'
% The main output is a melt rate, given in metres per year
%
% -------------------------------------------------------------------------
% To run PICO, use the following syntax:
%
% ab = PICO_driver(CtrlVar,MUA,GF,h,rhoi,rhow,PICO_opts)
%
% where:
%
%   CtrlVar: Control Variables (Ua structure)
%   MUA: Mesh (Ua structure)
%   GF: Grounded/floating mask (Ua structure)
%   h: ice thickness (MUA.Nnodes x 1)
%   rhoi: ice density (scalar)
%   rhow: water density (scalar)
%   PICO_opts: structure containing various PICO options (details below)
%
%% PICO_opts ==============================================================
%
% This is a structure similar to CtrlVar that contains all the options
% available when running the PICO model. It is not strictly necessary to
% define this structure when calling PICO but it is highly recommended,
% particularly for problem specific fields e.g. related to the ambient
% ocean fields. 
% -------------------------------------------------------------------------
% General options:
%
% - PICO_opts.algorithm: 'watershed' or 'polygon' (DEFAULT = 'watershed')
% There are two options related to how ice shelves are delineated. The
% 'watershed' options converts MUA into a structured grid and then uses
% image processing techniques to quickly define connected floating regions
% as individual ice shelves. The 'polygon' option creates polygons out of
% each individual GL segment and the MeshBoundaryCoordinates and defines
% individual ice shelves as floating nodes within these polygons. In
% general the 'watershed' option will be considerably faster while the
% 'polygon' option is much slower but generally more robust.
%
% - PICO_opts.SmallShelfMelt (DEFAULT = 0)
% The melt rate applied to floating nodes that for whatever reason are not
% within an ice shelf as determined by this code.
%
% - PICO_opts.FloatingCriteria: 'GLthreshold' or 'StrictDownstream' (DEFAULT = 'GLthreshold')
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
% Minimum number of floating nodes (polygon) or grid cells (watershed) for
% a given region to be considered an ice shelf.
%
% -------------------------------------------------------------------------
% Options related to the ambient ocean properties: 
%
% - PICO_opts.BasinsFile: Where the user wants to have multiple ocean
% basins with different properties, this should be the name of a .mat file
% containing a scatteredInterpolant covering the entire domain and numbered
% from 1 to n, where n is the maximum number of ocean basins. See Reese et
% al (2018) for an example of how to define the ocean basins. If 
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
% Options specific to PICO_opts.algorithm = 'watershed'
%
% - PICO_opts.PICOres (DEFAULT = 6e3/3e3/1e3 depending on mesh size)
% The resolution of the structured grid that the watershed option uses,
% this can potentially slow things down a lot and a higher number (lower
% resolution) typically seems to work fine.
%
% -------------------------------------------------------------------------
% Options specific to PICO_opts.algorithm = 'polygon'
%
% - PICO_opts.MeshBoundaryCoordinates
% This should be the same n x 2 matrix that is provided to
% Ua2D_InitialUserInput.m
