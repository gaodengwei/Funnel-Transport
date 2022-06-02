classdef vandpol<GaoSystem
    methods
        function obj = vandpol()
            mark = 0; 
            obj = obj@GaoSystem(mark,2,1);
        end
        function [xdot,df]= dynamics(obj,t,x,u)
            
            xdot = [-x(2);
                    -x(2)*(1-x(1)^2)+x(1)+u];
%                 xdot = [x(2);
%                    0.5*x(2)*(1-x(1)^2)-x(1)+u];
            if (nargout>1)
                df = [0 -1; 1+2*x(1)*x(2) -1+x(1)^2];
            end
%             quiver(x(1),x(2),0,u,0.5,'y')
        end
        
        % systme output
        function y = output(obj,t,x,u)
            y = x;
        end
 
    end
end