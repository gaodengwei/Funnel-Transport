classdef quadraticV_function
    % create by dengwei at 2017.12.4
    %  x'*S*x + x'*s1 + s2
    %  S can be doubles or trajectories (yielding a time-varying quadratic)
    
    properties
        S;
        s1;s2;
        p_t = msspoly('t',1);
        time_invariant_flag = true;
        dim;
        p_x;
        x0 = []; % the frame
        x  % mssploy
        invrho = 1;
    end
    
    methods
        function obj = quadraticV_function(S,s1,s2)
            n = length(s1);
            obj.dim=n; 
            obj.p_x = msspoly('x',n);obj.x = obj.p_x;
            % check the S and s1 and s2
            typecheck(S,{'double','Trajectory'});
            sizecheck(S,[n n]);
            obj.S = S; 
            if isa(S,'Trajectory') obj.time_invariant_flag = false; end
            
            if (nargin>1)
                typecheck(s1,{'double','Trajectory'});
                sizecheck(s1,[n 1]);
                obj.s1 = s1;
                if isa(s1,'Trajectory') obj.time_invariant_flag = false; end
            else
                obj.s1 = zeros(n,1);
            end
            
            if (nargin>2)
                typecheck(s2,{'double','Trajectory'});
                sizecheck(s2,1);
                obj.s2 = s2;
                if isa(s2,'Trajectory') obj.time_invariant_flag = false; end
            else
                obj.s2 = 0;
            end 
        end
        
        function ts = getBreaks(obj)
            ts = obj.S.getBreaks; 
        end
        
        function Vpoly = getPoly(obj,t)
            x = obj.p_x;
            if isTI(obj)
                Vpoly = x'*obj.S*x + x'*obj.s1 + obj.s2;
            else
                Vpoly = x'*obj.S.eval(t)*x + x'*obj.s1.eval(t) + obj.s2.eval(t);
            end
            Vpoly = Vpoly*obj.invrho;
        end
        
        function pVpt = getPolyTimeDeriv(obj,t)
            if isTI(obj)
                pVpt = 0;
            else
                x = obj.p_x;
                pVpt = x'*obj.S.deriv(t)*x + x'*obj.s1.deriv(t) + obj.s2.deriv(t);
            end
        end
        
        function V = inFrame(obj,x0)
                if  isa(x0,'polyniminalTrajectory')
                    D = -eye(obj.dim);  %feedback system -xbar
                    c = x0;
                    V = quadraticV_function(D'*obj.S*D,D'*(obj.S+obj.S')*c + D'*obj.s1,c'*obj.S*c + c'*obj.s1 + obj.s2); 
                    V.x0 = x0;
                elseif isa(x0,'double')
                    obj.x0 = x0; 
                    obj.time_invariant_flag = true;
                    V = obj;
                elseif isa(x0,'msspoly')
                    obj.x0 = x0; 
                    obj.time_invariant_flag = true;
                    V=obj;
                else
                    error('cannot build such frame')
                end  
        end
        
        function V = mtimes(a,b)
            % support simple scaling of Lyapunov functions via multiplication by
            % a (scalar) double
            if ~isa(b,'quadraticV_function')
                % then a must be the lyapunov function.  swap them.
                tmp=a; a=b; b=tmp;
            end
            typecheck(a,{'numeric','Trajectory'});
            try
                V = quadraticV_function(a*b.S, a*b.s1, a*b.s2);
            catch
                V = quadraticV_function(a*b.S, b.s1, b.s2);
            end
        end 
        
        function obj = eval(obj,t,x)
            if isTI(obj) 
                t=0;
            elseif nargin<2 || isempty(t)
                error('you must specify a time'); 
            end
            obj = double(subs(obj.getPoly(t),[obj.t;obj.x],[t;x]));
        end 
        
        function V = extractquadraticV_function(obj)
            if (~isTI(obj)) error('not implemented yet'); end
            
            Vpoly = obj.getPoly();
            x = obj.p_x;
            
            if (deg(Vpoly,x)>2) error('not quadratic'); end
            
            Ss=double(.5*subs(diff(diff(Vpoly,x)',x),x,0*x));
            s1s=double(subs(diff(Vpoly,x),x,0*x))';
            s2s=double(subs(Vpoly,x,0*x));
            
            V = quadraticV_function(Ss,s1s,s2s);
        end
        function obj = mrdivide(obj,b) 
            typecheck(b,'double');
            obj.invrho = 1/b;
        end
    end
    
    methods 
        function tf = isTI(obj)
            % Returns true if the system is time-invariant
            tf = obj.time_invariant_flag;
        end
        
    end
end


