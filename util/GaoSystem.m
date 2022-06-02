% /*! @GaoSystem.m
% *************************************************************************
% <PRE>
% file.name       : GaoSystem.m
% related files   :
% function&ablity : build system
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/8/29     1.00        gaodengwei   build
% 2016/11/10    1.00        gaodengwei   modify: extend to highest class
% 2017/12/04    1.00        gaodengwei   modify: add time_invariant_flag
% any dynamic model
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
%
% output:
% *************************************************************************
classdef GaoSystem
    %% ==============dynmaic system======================
    properties
        mark;
        eps = 1e-10;        % precision of the system
        timenow = 0;        % simuliation start time
        end_step = 300;     % simulantion total step
        temp = 1;           % simulantion step
        timeNum = 100;      % funnel section between start and goal
        breaks;             % time step in calculaion
        getNumStates;       % dynamic dimension
        getNumInput;        % control dimension
        getNumoutput;       % output dimension
        getNumObstacles;    % number of obstacle
        A; B;   % linearlize coefficient
        K;      % feedback gain
        INPUTS;  % system set
        SetVariable;  % system control by gaosystem
        % plot part in simulition
        PlotObj;     % plot function
        parSet;
        RealMap;     % real map
        FunTraj;     % nominal trajetcory function
        FunContr;    % nominal control function
        Maxinterval; % time interval
        umin;umax;   % control sat
        p_t;p_u;p_x; % messploy coefs
        time_invariant_flag =false;
    end
    % ========================== start method=========================
    methods
        function obj = GaoSystem(dynamic_mark,numx,numu)
            obj.getNumStates = numx;      % Number of States
            obj.getNumoutput = numx;      % number of output
            obj.getNumInput = numu;       % Number of Inputs
            if dynamic_mark>5
                error('need to add new dynamic function');
            end
            switch (dynamic_mark)
                case 0 %LTI system
                    obj.mark = 0;
                    obj.time_invariant_flag = true;
                case 1 %'car'
                    obj.mark = 1;
                case 2 %'VideoRay2D'
                    obj.mark = 2;
                case 3 %'Sapcecraft'
                    obj.mark = 3;
                case 4  %VideoRay3D
                    obj.mark = 4;
                otherwise
                    obj.mark = 5;
            end
            obj.INPUTS = [];
            obj.getNumObstacles=[];
            obj.A = [];
            obj.B = [];
            obj.K = [];
            obj.SetVariable.flag_rho = 1;%1:rho calculate\0:rho load
            obj.SetVariable.ROAflag = 2;% 1:sample;2:SOS
            obj.PlotObj = [];
            obj.parSet = [];
            obj.FunTraj = [];
            obj.FunContr = [];
            % messploy coef
            obj.p_t = msspoly('t');
            obj.p_u = msspoly('u',numu);
            obj.p_x = msspoly('x',numx);
            
            % control sat
            obj.umin = [];
            obj.umax = [];
        end
        % =================================================================
        % plot figure
        function obj = PlotFigure(obj,varargin)
            if isempty(varargin)
                return
            else
                plotfigurecode(obj,varargin{1})
            end
        end
        
        function obj = BuildMap(obj,map,shrinkValue)
            % draw the map in 3D
            realSize = size(map)/shrinkValue;
            height = 10;
            accuracyControl = 1;
            I=imresize(map,realSize*accuracyControl);
            [x,y]=size(I);
            X=1:x;
            Y=1:y;
            [xx,yy]=meshgrid(X/accuracyControl,Y/accuracyControl);
            hh=height*flipud(I);
            x = reshape(xx,size(xx,1)*size(xx,2),1);
            y = reshape(yy,size(yy,1)*size(yy,2),1);
            z = reshape(hh,size(hh,1)*size(hh,2),1);
            hadleMap = figure('Visible', 'off');
            trisurf(delaunay(x, y), x, y, z);
            colormap copper;                 % Default color map.
            set(gca, 'Position', [0 0 1 1]); % Fill the figure window.
            axis equal vis3d;                % Set aspect ratio.
            shading interp;                  % Interpolate color across faces.
            camlight left;                   % Add a light over to the left somewhere.
            lighting gouraud;                % Use decent lighting.
            obj.RealMap = hadleMap;
        end
        
        function sys = updateNominal(sys,OL_time,OL_state,OL_control)
            %             OL_time : open loop time step
            %             OL_state: open loop state
            %             OL_control: open loop control
            
            if length(OL_state) == 1
                sys.FunTraj = OL_state;                     % goal points
                sys.FunContr = OL_control;                  % open-loop control
            else
                [m,n] = size(OL_state);
                [mc,nc] = size(OL_control);
                tm = length(OL_time);
                if m~=tm
                    if n~=tm
                        error('state is not find in time')
                    else
                        OL_state = OL_state';
                    end
                end
                if mc ~=tm||nc~=sys.getNumInput
                    error('wrong control')
                end
                statefun = spline(OL_time,OL_state');               % openloop state for the nominal trajectory
                controlfun = spline(OL_time,OL_control');           % openloop control for the nominal trajectory
                sys.Maxinterval = [0 max(OL_time)];                 % cost time
                sys.FunTraj = polyniminalTrajectory(statefun);      % nonminal trajectory
                sys.FunContr = polyniminalTrajectory(controlfun);   % nonminal control
            end
        end
        
        function obj = setInputLimits(obj,umin,umax)
            % Guards the input limits to make sure it stay consistent
            
            if (isscalar(umin)), umin=repmat(umin,obj.getNumInput,1); end
            if (isscalar(umax)), umax=repmat(umax,obj.getNumInput,1); end
            
            sizecheck(umin,[obj.getNumInput,1]);
            sizecheck(umax,[obj.getNumInput,1]);
            if (any(obj.umax<obj.umin)), error('umin must be less than umax'); end
            obj.umin = umin;
            obj.umax = umax;
        end
        
        function obj = tv_poly_linearize(obj)
            % linearize the system and get the coefs in every points
            x0 = obj.FunTraj;
            u0 = obj.FunContr;
            if size(x0.eval(0),2) ~= 1, error('x0 must be a column.'); end
            if size(u0.eval(0),2) ~= 1, error('u0 must be a column.'); end
            
            obj.A = @(t)dynamicA(obj,t,x0.eval(t),u0.eval(t));
            obj.B = @(t)dynamicB(obj,t,x0.eval(t),u0.eval(t));
        end
        
        function [A,B,C,D,x0dot,y0] = linearize(obj,t0,x0,u0)
            % Uses the geval engine to linearize the model around the nominal
            % point, at least for the simple case.
            
            nX = obj.getNumStates;
            nU = obj.getNumInput;
            [~,df] = geval(@obj.dynamics,t0,x0,u0);
            A = df(:,1+(1:nX));
            B = df(:,nX+1+(1:nU));
            
            if (nargout>2)
                [~,dy] = geval(@obj.output,t0,x0,u0);
                C = dy(:,1+(1:nX));
                D = dy(:,nX+1+(1:nU));
                if (nargout>4)
                    x0dot = dynamics(obj,t0,x0,u0);
                    if (nargout>5)
                        y0 = output(obj,t0,x0,u0);
                    end
                end
            end
        end
        
        function sys = timecalculation(sys,varargin)
            % set time breaks
            if size(varargin)>0
                Num = varargin{1};
            else
                Num = sys.timeNum;
            end
            if sys.SetVariable.flag_rho == 1
                taus = (linspace(sys.Maxinterval(1),sys.Maxinterval(2),Num))'; % set the funnel with some intervals, the num is time count
            else
                taus = (linspace(sys.Maxinterval(1),sys.Maxinterval(2),Num))';
                warning('Do not calculate the funnel')
            end
            sys.timeNum = Num;
            sys.breaks = taus;
        end
        
        function dynamicf = getDynamics(obj,t)
            % get dynamic handle
            t=max(min(t,obj.Maxinterval(end)),obj.Maxinterval(1));
            dynamicf = obj.handle(t);
        end
        
        % system output is the the state
        function y = output(obj,t,x,u)
            % default full state feedback
            y = x;
        end
        
        function obj = setfield(obj,field)
            %             add field to system
            obj.INPUTS.field = field;
        end
        
        function tf = isTI(obj)
            % Returns true if the system is time-invariant
            tf = obj.time_invariant_flag;
        end
        
        function f = TransSystem(obj,varargin)
            % trans sys to V.x0 or nominal trajectory
            x = obj.p_x; u =obj.p_u; t=obj.p_t;
            if isTI(obj) % time invarant
                V = varargin{1};
                if ~isa(V,'quadraticV_function'), error('V wrong'); end
                f = obj.dynamics(t,x,u);
                if ~isempty(V.x0)
                    f = TransFormSystem(f,x,V.x0);
                end
            else %time varying
                %                trans to nominal trajectory
                order = varargin{1};
                [f, ~]= TaylorExpansion(sys,order); % taylor expantion in order
            end
        end
    end
end






