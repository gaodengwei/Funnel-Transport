classdef Pendulumdyn<GaoSystem
    properties
        b = 0.1;
        g=9.8;
        xG = [pi;0];
        uG=0;
        Q;R;uppercontrol
    end
    
    methods
        function obj = Pendulumdyn()
            mark = 0;
            obj = obj@GaoSystem(mark,2,1);
        end
        function [xdot,df,d2f,d3f]= dynamics(obj,t,x,u)
            xdot = [x(2);
                u-obj.b*x(2)-obj.g*sin(x(1))];
            if (nargout>1)
                %                 df = [0 1; -obj.g*cos(x(1)) -obj.b];
                % Compute Gradients:
                df = sparse(2,4);
                df(1,3) = 1;
                df(2,2) = -(obj.g*cos(x(1)));
                df(2,3) = -obj.b;
                df(2,4) = 1;
            end
            if (nargout>2)
                d2f = [0 0; obj.g*sin(x(1)) 0];
            end
            if (nargout>3)
                d3f = [0 0; obj.g*cos(x(1)) 0];
            end
        end
        
%         function xdot = dynamics_par(obj,t,x,u,g)
%             xdot  = dynamics(obj,t,x+g,u); 
%         end
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function A = dynamicA(obj, t, x, u)
            [~,df,~,~] = obj.dynamics(t,x,u);
            nX = obj.getNumStates;
            A = full(df(:,1+(1:nX)));
        end
        
        function B = dynamicB(obj, t, x, u)
            [~,df,~,~] = obj.dynamics(t,x,u);
            nX = obj.getNumStates;
            nU = obj.getNumInput;
            B = full(df(:,nX+1+(1:nU)));
        end
        
    end
end