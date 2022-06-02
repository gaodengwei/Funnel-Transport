% /*! @V_function.m
% *************************************************************************
% <PRE>
% file.name       : V_function.m
% related files   :
% function&ablity : Calculate Funnel
% author          : gaodengwei
% version         : 1.00
% --------------------------------------------------------------------------------
% remarks         :
% --------------------------------------------------------------------------------
% record of modify :
% date          version     name         content
% 2016/7/28     1.00                     build
% </PRE>
% ********************************************************************************
%
% * right(c)
%
% *************************************************************************
% input:
% Pmess: the invariant set in every nonminal points
% K: feedback gain function
% t: msspoly time variable
% taus: nominal trajectory sample - based points
% x: msspoly state variable
% xmess: nominal trajectory function
% u: msspoly control variable
% ui: nominal trajectory control function
%
% output:
% *************************************************************************
classdef V_function
    %  x'*S*x
    %  S can be doubles or trajectories (yielding a time-varying quadratic)
    
    properties
        Spp;
        getbreak;
        x;t;
        x0;
        w;
        Vpoly; %default is []
    end
    
    methods
        function obj = V_function(t,x,S,taus,varargin)
            obj.x = x;                      % variable in lyapunov function
            obj.t=t;                        % variable in lyapunov function
            obj.getbreak = taus;            % time break
            obj.Spp = S;                    % spline form matrix
            obj.x0 = [];
            if nargin>4
                obj.w = varargin{1};
            else
                obj.w = [];
            end
        end
        
        function obj = eval(obj,t,x)
            obj = double(subs(obj.getPoly(t),[obj.t;obj.x],[t;x]));
        end
        
        function Vpoly = getPoly(obj,varargin)
            %             if nargin<2
            %                 if ~isempty(obj.x0)
            %                     Vpoly = (obj.x-obj.x0)'*obj.Spp*(obj.x-obj.x0);
            %                 else % default is origin
            %                     Vpoly = obj.x'*obj.Spp*obj.x;
            %                 end
            %             else
            %                 tt=varargin{1};
            %                 if ~isempty(obj.x0)
            %                     Vpoly = (obj.x-obj.x0.eval(tt))'*ppval(obj.Spp,tt)*(obj.x-obj.x0.eval(tt));
            %                 else
            %                     Vpoly = obj.x'*ppval(obj.Spp,tt)*obj.x;
            %                 end
            %             end
            if isempty(obj.Vpoly)
                if nargin<2
                    Vpoly = obj.x'*obj.Spp*obj.x;
                else
                    tt=varargin{1};
                    Vpoly = obj.x'*ppval(obj.Spp,tt)*obj.x;
                end
            else
                if nargin<2
                    Vpoly = obj.Vpoly;
                else
                    tt=varargin{1};
                    Vpoly = obj.Vpoly(tt);
                end
            end
        end
        
        function pVpt = getPolyTimeDeriv(obj,t)
            % calculate the time derivative for a poly
            [b,c,~,k,d] = unmkpp(obj.Spp);
            if (k==1)  % handle the case of too-high order
                ppform = mkpp(b,0*c(:,1),d);
            else
                for i=1:k-1
                    cnew(:,i) = (k-i)*c(:,i);
                end
                ppform = mkpp(b,cnew,d);
            end
            df = ppval(ppform,t);
            pVpt = obj.x'*df*obj.x;
        end
        
        function obj = updateV(obj,a)
            a = polyniminalTrajectory(a);
            b = polyniminalTrajectory(obj.Spp);
            newS = a*b;
            xp = obj.x0;  % add frame if exist
            obj = V_function(obj.t, obj.x, newS.pp, obj.getbreak);
            obj.x0 = xp;
        end
        
        function obj = inFrame(obj,x0)
            obj.x0 = x0;
            if isa(x0,'polyniminalTrajectory')
                obj.Vpoly = @(t)subss(obj.getPoly(t),obj.x,obj.x-x0.eval(t));
            else
                obj.Vpoly = subss(obj.getPoly,obj.x,obj.x-x0);
            end
        end
        
        function output = PdiffVE(obj,E2,t,dim)
            n = length(obj.x);
            m = dimension(E2);
            x00 = obj.x0.eval(t);
            f = obj.getPoly(t);
            xx = obj.x;
            %             if m ~= n
            %                 error('Pontryagin difference must be of the same dimension.');
            %             end
            % make the PP as symmetric
            no_dims=1:length(x00);  no_dims(dim)=[];
            f = subs(f,xx(no_dims),x00(no_dims));
            xx = xx(dim);
            H = double(0.5*diff(diff(f,xx)',xx));
            PP = H^-1;
            E1 = ellipsoid(x00(dim),PP);
            
            % compete difference
            if  length(dim)>3
                bound = Pontryagin_diff(E1,E2);
            else
                bound = minkdiff(E1,move2origin(E2));
                output = bound;
                return
            end
            if dim == m
                output = bound;
            else
                % convex hull
                k = convhulln(bound');
                outputdim = bound(1:dim,k');
                m = convhull(outputdim');
                output = outputdim(:,m);
            end
            output = output+repmat((x00(1:dim,:)),1,size(output,2));
        end
        
        function obj = updatebreaks(obj,ts)
            obj.getbreak = ts;
            [b,c,~,~,d] = unmkpp(obj.Spp);
            obj.Spp = mkpp(b+ts(1),c,d);
            if (b(end)-b(1))~=(ts(end)-ts(1))
                error('check the funnel shift')
            end
            % obj.Spp = spline(obj.Spp.breaks+ts(1),ppval(obj.Spp,obj.Spp.breaks));
        end
        
        function b = mtimes(a,b)
            % support simple scaling of Lyapunov functions via multiplication by
            % a (scalar) double
            if ~isa(b,'V_function')
                % then a must be the lyapunov function.  swap them.
                tmp=a; a=b; b=tmp;
            end
            typecheck(a,'numeric');
            sizecheck(a,1);
            
            b.Vpoly = a*b.Vpoly;
        end
        
        function S = extractVFunction(obj,t)
            Vpoly = obj.getPoly(t);
            
            if (deg(Vpoly,obj.x)>2) error('not quadratic'); end
            
            S=double(.5*subs(diff(diff(Vpoly,obj.x)',obj.x),obj.x,0*obj.x));
        end
        
    end
end


