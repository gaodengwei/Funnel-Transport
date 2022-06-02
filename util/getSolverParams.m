function options = getSolverParams(SDPsolver)
% GETSOLVERPARAMS Return sdp solver parameters.
%
% Returns solver parameters as a Yalmip sdpsettings structure
% Parameter selected only for MOSEK and SeDuMi, otherwise using YALMIP's
% default ones.
%
% INPUTS:
%
%   "SDPsolver" -
%
% OUTPUTS:
%
%   "" -
%
% NOTES:
%
% - This helper function is from: ROA - Matlab codes (including distributions of YALMIP and GloptiPoly 3, and using either MOSEK or SeDuMi which should be installed) for computing estimates of the region of attraction of a polynomial control system, following the developments described in D. Henrion, M. Korda, Convex computation of the region of attraction of polynomial control systems, LAAS-CNRS Research Report 12488, August 2012. Developed by Milan Korda.
% See Gloptipoly library and ROA Matlab codes in
% http://homepages.laas.fr/henrion/
% and references therein.
%

%==========================================================================

if(strcmpi(SDPsolver,'mosek'))
    if(~exist('mosekopt.m'))
        disp('***MOSEK not found, leaving the solver selection to Yalmip***')
        SDPsolver = '';
    end
end

if(strcmpi(SDPsolver,'sedumi'))
    if(~exist('sedumi.m'))
        disp('***SeDuMi not found, leaving the solver selection to Yalmip***')
        SDPsolver = '';
    end
end

switch SDPsolver
    case 'mosek'
        disp('SDP Solver: MOSEK')
        options = sdpsettings('solver','mosek-sdp',...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS',1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-12);
    case 'sedumi'
        disp('SDP Solver: SeDuMi')
        options = sdpsettings('solver','sedumi','sedumi.eps',1e-12);
    case ''
        options = [];
    otherwise
        disp(['SDP Solver: ' SDPsolver])
        options = sdpsettings('solver',SDPsolver);
end
