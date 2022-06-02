classdef feedbackSys<GaoSystem
    properties
        sys
        xG
        uG
        ubar = [];
        xbar = [];
        y0
    end
    methods
        function obj = feedbackSys(sys,K,B,u0,x0)
            obj = obj@GaoSystem(sys.mark,sys.getNumStates,sys.getNumInput);
            if iscell(K)
                obj.K = K{1};obj.y0 = K{2};
            else
                obj.K = K;obj.y0 = []; 
            end 
                
            obj.B = B;
            obj.sys = sys;
            if nargin >3,obj.ubar = u0; end
            if nargin >4,obj.xbar = x0; end
            try
                obj.xG = sys.xG;
                obj.uG = sys.uG;
            catch
                % default, be careful here 
                obj.xG = x0;
                obj.uG = u0;
            end
        end
        
        function [xdot,df]= dynamics(obj,t,x,~)
            if isa(obj.ubar,'Trajectory')
                if isa(obj.K,'function_handle')
                    u = obj.ubar.eval(t)-obj.K(x(obj.getNumStates));
                elseif isa(obj.K,'Trajectory')
                    u = obj.ubar.eval(t)+obj.K.eval(t)*(x-obj.xbar.eval(t));
                end
            else
                % odd  here need to be update
                if isa(obj.K,'function_handle')
                    u = -obj.K(x(obj.getNumStates))*x;
                elseif isa(obj.K,'Trajectory')
                    u = obj.K.eval(t)*x; 
                else
                    u = -obj.K*(x-obj.xG);
                end
            end
             
            xdot = obj.sys.dynamics(t,x,u);

            if (nargout>1)
                [~,df0] = obj.sys.dynamics(t,x,0);
                df = df0- obj.B.eval(t)*obj.K.eval(t);
            end
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
 
    end
end