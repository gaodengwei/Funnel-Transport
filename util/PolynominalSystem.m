classdef PolynominalSystem<GaoSystem
    % Implements
    %   xcdot = Ac*x + Bc*u + xcdot0
    %   y = C*x + D*u + y0
    properties
        Gaosys
        sysind
        C;D;y0;xcdot0;
        num_x;num_u;num_y;
        
    end
    
    methods
        function obj = PolynominalSystem(sys,xcdot0,Ac,Bc,C,D,y0,flagSys)
            num_x = max([size(Ac,1),size(Bc,1),size(xcdot0,1)]);
            num_u = max([size(Bc,2),size(D,2)]);
            num_y = max([size(C,1),size(D,1),size(y0,1)]);
            obj = obj@GaoSystem(sys.mark, sys.getNumStates, sys.getNumInput);
            obj.num_x = num_x;
            obj.num_y = num_y;
            obj.num_u = num_u;
            if flagSys == 1;
                obj.Gaosys = sys;
            else
                obj.Gaosys = [];
            end
            obj.sysind = stateIndicesForCombination(obj);
            if (isempty(Ac))
                obj.A = sparse(num_x,num_x);
            else
                sizecheck(Ac,[num_x,num_x]);
                if isnumeric(Ac), Ac = ConstantTrajectory(Ac); end
                obj.A = Ac;
            end
            if (isempty(Bc))
                obj.B = sparse(num_x,num_u);
            else
                sizecheck(Bc,[num_x,num_u]);
                if isnumeric(Bc), Bc = ConstantTrajectory(Bc); end
                obj.B = Bc;
            end
            if (isempty(C))
                obj.C = sparse(num_y,num_x);
            else
                sizecheck(C,[num_y,num_x]);
                if isnumeric(C), C = ConstantTrajectory(C); end
                obj.C = C;
            end
            if (isempty(D))
                obj.D = sparse(num_y,num_u);
            else
                sizecheck(D,[num_y,num_u]);
                if isnumeric(D), D = ConstantTrajectory(D); end
                obj.D = D;
            end
            if (isempty(y0))
                obj.y0 = sparse(num_y,1);
            else
                sizecheck(y0,[num_y,1]);
                if isnumeric(y0), y0 = ConstantTrajectory(y0); end
                obj.y0 = y0;
            end
            if (isempty(xcdot0))
                obj.xcdot0 = sparse(num_x,1);
            else
                sizecheck(xcdot0,[num_x,1]);
                obj.xcdot0 = xcdot0;
            end
        end
        
        function xcdot = dynamics(obj,t,x,u)
            [x1,x2]=decodeX(obj,x);
            [y1,y2]=getOutputs(obj,t,x,u);
            xcdot=dynamics(obj.Gaosys,t,x1,sat1(obj,y2+u));    % sys with feedback control
        end
        
        function [y,dy]= output(obj,t,x,u)
            y=zeros(obj.num_y,1);
            if (obj.num_x)
                y=y+eval(obj.C,t)*x;
            end
            if (obj.num_u)
                y=y+eval(obj.D,t)*u;
            end
            if ~isempty(obj.y0)
                y=y+eval(obj.y0,t);
            end
            if (nargout>1)
                dy = zeros(obj.num_y,1+obj.num_x+obj.num_u);
                if (obj.num_x)
                    dy(:,1:1+obj.num_x)=dy(:,1:1+obj.num_x) + [deriv(obj.C,t)*x, eval(obj.C,t)];
                end
                if (obj.num_u)
                    dy(:,[1,2+obj.num_x:end])=dy(:,[1,2+obj.num_x:end])+[deriv(obj.D,t)*u,eval(obj.D,t)];
                end
                if ~isempty(obj.y0)
                    dy(:,1)=dy(:,1)+deriv(obj.y0,t);
                end
            end
        end
        
        function sys = cascade(sys1,sys2)
            % combine sys1 and sys2: output of sys1 is the input of sys2
            % x1dot = A1x1 + B1u + xcdot01
            % x2dot = A2x2 + B2(C1x1+D1u+y01) + xcdot02
            % y = C2x2 + D2(C1x1+D1u + y01) + y02
            newA=[]; newC=[];
            newB = [sys1.B; sys2.B*sys1.D];
            newxcdot0 = [sys1.xcdot0; sys2.B*sys1.y0 + sys2.xcdot0];
            
            newD = sys2.D*sys1.D;
            newy0 = sys2.D*sys1.y0 + sys2.y0;
            
            sys = PolynominalSystem(sys2.Gaosys,newxcdot0,newA,newB,newC,newD,newy0,1);
        end
    end
    
    % function inside this class
    methods (Access=private)
        function [x1,x2] = decodeX(obj,x)
            x1=x(obj.sysind);
            x2=x([]);
        end
        function u1=sat1(obj,u1)
            if ~isempty(obj.Gaosys.umin)||~isempty(obj.Gaosys.umax)
                u1=min(max(u1,obj.Gaosys.umin),obj.Gaosys.umax);
            end
        end
        function u2=sat2(obj,u2)
            if ~isempty(obj.Gaosys.umin)||~isempty(obj.Gaosys.umax)
                u2=min(max(u2,obj.Gaosys.umin),obj.Gaosys.umax);
            end
        end
        function sysind = stateIndicesForCombination(sys)
            ind=0;
            n=sys.getNumStates;
            sysind = ind+(1:n)';
        end
        function [y1,y2] = getOutputs(obj,t,x,u)
            [x1,x2]=decodeX(obj,x);
            y1=output(obj.Gaosys,t,x1,u);       % output shouldn't depend on u
            y2=output(obj,t,x2,sat2(obj,y1));   % u=u0+Kx
        end
        
    end
end