classdef DAODE
    % @copyright dengwei version 1 at 2019
    % <Funnel Transport: Approximating Finite-Horizon Reachable Sets>
    properties
        vars            % number of variables
        uncertain = 0   % default 0
        u               % feedback control (optional)
        order           % polynomial maximum degree  
        LsFlag = 0;     % 0: direct remove; 1:QP remove; 2: lsqcurvefit remove
        R = 0.5;        % initial radius default
    end
    properties(Hidden)
        sys = []; 
        % Initialize method parameters. 
        A = [1/5, 3/10, 4/5, 8/9, 1, 1];
        B = [1/5         3/40    44/45   19372/6561      9017/3168       35/384
            0           9/40    -56/15  -25360/2187     -355/33         0
            0           0       32/9    64448/6561      46732/5247      500/1113
            0           0       0       -212/729        49/176          125/192
            0           0       0       0               -5103/18656     -2187/6784
            0           0       0       0               0               11/84];
    end
    
    methods
        function obj = DAODE(x,u,order,R,methods)
            obj.vars = x;
            obj.u = u;
            obj.order = order; 
            if nargin>3, obj.R = R; end
            if nargin > 4
                obj.LsFlag = methods;
            end
        end
        
        function zdot = dynamics(obj,t,x)
            % it is not a general way to add the uncertain, 
            % more general thing need to do here
            zdot = subs(obj.sys(t),obj.vars,x)+obj.uncertain;
        end
        
        function obj = addNoise(obj,w)
            obj.uncertain = w;
        end
         
        function xT = integrate(obj,tspan,nstep)
            xm = msspoly('x',length(obj.vars));  
            x = xm;
            xs = sym('x',[1,length(x)],'real')';  
            t = tspan(1);
            Hspan = abs(tspan(2)-tspan(1));
            tdir = sign(tspan(2) - tspan(1)); 
            h = tdir*Hspan/nstep;
            neq = length(obj.vars);
            f = msspoly(zeros(neq,6));
            xT = x;
            while true
                odeFcn = @obj.dynamics;
                try % find a ploy expression
                    f0 = odeFcn(t,x);
                    f(:,1) = cleanPow(f0,obj.order,obj.LsFlag,xs,obj.R);
                catch
                    keyboard
                    f0 = odeFcn(t,xs); 
                    f(:,1) = cleanPow(f0,obj.order,obj.LsFlag,xs,obj.R);
                    odeFcn = @(t,a)subs(f(:,1),xm,a); 
                end
                
                hA = h * obj.A;
                hB = h * obj.B;
                f(:,2) = cleanPow(odeFcn(t+hA(1),x+f*hB(:,1)),obj.order,obj.LsFlag,xs,obj.R);
                f(:,3) = cleanPow(odeFcn(t+hA(2),x+f*hB(:,2)),obj.order,obj.LsFlag,xs,obj.R);
                f(:,4) = cleanPow(odeFcn(t+hA(3),x+f*hB(:,3)),obj.order,obj.LsFlag,xs,obj.R);
                f(:,5) = cleanPow(odeFcn(t+hA(4),x+f*hB(:,4)),obj.order,obj.LsFlag,xs,obj.R);
                f(:,6) = cleanPow(odeFcn(t+hA(5),x+f*hB(:,5)),obj.order,obj.LsFlag,xs,obj.R);
                
                x = x + f*hB(:,6);
                t = t + hA(6); 
                x = cleanPow(x,obj.order,obj.LsFlag,xs,obj.R);
                xT = [xT,x];
                if tdir>0
                    if t >= tspan(2)
                        break;
                    end
                else
                    if t <= tspan(2)
                        break;
                    end
                end
            end
        end
         
        function obj = TaylorMap(obj,plant,xtraj,utraj,K) 
            % Taylor mapping the ploynomianl
            if isa(xtraj,'Trajectory')
                if nargin<5
                    K{1}= ConstantTrajectory(zeros(length(obj.u),length(obj.vars))); 
                    K{2} = ConstantTrajectory(zeros(length(obj.vars),1)); 
                end
                obj = TaylorMapTimevarying(obj,plant,xtraj,utraj,K);
            elseif isa(xtraj,'double')
                if nargin<5, K=ConstantTrajectory(zeros(length(obj.u),length(obj.vars))); end
                obj = TaylorMapFix(obj,plant,xtraj,utraj,K);
            else
                error('no such class')
            end
        end
        %% Taylor Map
        function obj = TaylorMapTimevarying(obj,plant,xtraj,utraj,K) 
            % time varying mapping
            p_xu = [obj.vars;obj.u]; 
            ctf = PolynominalSystem(plant,[],[],[],[],eye(plant.getNumStates),-xtraj,0); 
            dtf = PolynominalSystem(plant,[],[],[],[],eye(plant.getNumStates),xtraj,0); 
            ltvsys = PolynominalSystem(plant,[],[],[],[],K{1},K{2},1); 
            cInput = PolynominalSystem(plant,[],[],[],[],eye(plant.getNumInput),utraj,1); 
            feedbacksys = cascade(ctf,ltvsys); 
            feedbacksys = cascade(feedbacksys,cInput); 
            xdothat = @(t)obj.build_poly(@feedbacksys.dynamics,t,xtraj.eval(t),K{2}.eval(t),obj.order,p_xu); 
            zdothat = @(t)obj.newdyn(t,ctf,dtf,xdothat);
            obj.sys = zdothat;
        end
        
        function obj = TaylorMapFix(obj,plant,x0,u0,K)
            % time invariant mapping 
            contr = u0+K*(obj.vars-x0);
            try
                xdot = @(t)plant.dynamics(t,obj.vars,contr);
            catch
                p_xu = [obj.vars;obj.u]; 
                xdot = @(t)obj.build_poly(@plant.dynamics,0,x0,u0,obj.order,p_xu,contr);
            end 
             
            if any(x0~=0)
                zdothat = @(t)TransFormSystem(xdot(t),x,x0);
            else
                zdothat = xdot;
            end
            obj.sys = zdothat;
        end
        
        %% General Map by FT
        function obj = GeneralMap(obj,plant,xtraj,utraj,K)
              
            if isa(xtraj,'Trajectory')
                if nargin<5
                    K{1}= ConstantTrajectory(zeros(length(obj.u),length(obj.vars))); 
                    K{2} = ConstantTrajectory(zeros(length(obj.vars),1)); 
                end
                obj = GeneralMapTimevarying(obj,plant,xtraj,utraj,K);
            elseif isa(xtraj,'double')
                if nargin<5, K=ConstantTrajectory(zeros(length(obj.u),length(obj.vars))); end
                obj = GeneralMapFix(obj,plant,xtraj,utraj,K);
            else
                error('no such class')
            end 
        end
         
        function obj = GeneralMapFix(obj,plant,x0,u0,K)
            % time invariant mapping 
            contr = u0+double(K)*(obj.vars-x0); 
            try
                plant.dynamics(0,obj.vars,contr);
                xdot = @(t)plant.dynamics(t,obj.vars,contr);  
            catch
                p_xu = [obj.vars;obj.u]; 
                xdot = @(t)obj.build_poly(@plant.dynamics,0,x0,u0,obj.order,p_xu,contr); 
            end  
             
            if any(x0~=0)
                zdothat = @(t)TransFormSystem(xdot(t),x,x0);
            else
                zdothat = xdot;
            end
            obj.sys = zdothat; 
        end
         
        function obj = GeneralMapTimevarying(obj,plant,xtraj,utraj,K)
            % time varying mapping
            xs = sym('x',[1,length(obj.vars)],'real')';
%             us = sym('u',[1,length(obj.u)],'real')'; 
            ctf = PolynominalSystem(plant,[],[],[],[],eye(plant.getNumStates),-xtraj,0); 
            dtf = PolynominalSystem(plant,[],[],[],[],eye(plant.getNumStates),xtraj,0);  
            Control = @(t,x)(utraj.eval(t)+K{1}.eval(t)*(x-xtraj.eval(t))+K{2}.eval(t));
            xdot = @(t)obj.build_fun(@plant.dynamics,t,xs,Control(t,xs),obj.order,obj.R); 
            obj.vars = xs;
            obj.LsFlag = 3;   % enforce the truction
            zdothat = @(t)obj.newdyn(t,ctf,dtf,xdot);  
            obj.sys = zdothat; 
        end
    end
     
    methods(Access = private)
        function xdot = build_fun(obj,plant,t,x,u,order,R) 
            numx = length(x); numu = length(u);
            xdot = plant(t,x,u);
%             xdot = subs(p,obj.u,contr);
        end
        
        function xdot = build_poly(obj,plant,t,x0,u0,order,p_xu,contr)
            if nargin < 8, contr=zeros(size(obj.u)); end
            nX = length(x0); nU=length(u0);
            xu0 = [x0;u0];
            xu = TaylorVar.init(xu0,order);
            x0 = xu(1:nX);
            u0 = xu(nX+(1:nU));
            p = getmsspoly(plant(t,x0,u0),p_xu-xu0);
            xdot = subs(p,obj.u,contr);
        end 
        
        function g = newdyn(obj,t,ctf,dtf,xdothat)
            x = obj.vars;
            z = obj.vars;
            [~,dc]=ctf.output(t,[],x);
            d = dtf.output(t,[],z);
            dcdt = dc(:,1); dcdx = dc(:,2:end);
            g = dcdt + dcdx*xdothat(t);
            g = subs(g,x,d); 
        end
         
    end
end


