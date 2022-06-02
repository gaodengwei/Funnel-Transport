classdef (InferiorClasses = {?ConstantTrajectory}) polyniminalTrajectory < Trajectory
% this function creat by dengwei at 2017.2.4 
% this function updates by dengwei at 2017.12.4
% methods include:
% loadobj:load the exist polyniminalTrajectory
% eval:value of polyniminalTrajectory
% fnder: Derivative for input polyniminalTrajectory
% deriv: Derivative at particular time
% shiftTime: shift Time of polyniminalTrajectory
% cutTime: cut timespan [t1 t2] as new polyniminalTrajectory
% scaleTime: compress the polyniminalTrajectory with time sclae
% approxNearestPoint: find the nearest point of polyniminalTrajectory
% append: Append a polyniminalTrajectory to this one, creating a new trajectory that starts where this object starts and ends where the given trajectory ends.
% refine: insert the break of trajectory    
% rebreaks: redefine the break of trajectory    
properties
        pp
    end
    
    methods
        function obj = polyniminalTrajectory(ppform)
            if isnumeric(ppform)
                ppform = ConstantTrajectory(ppform);
            end
            if isa(ppform,'ConstantTrajectory')
                ppform = mkpp([-inf,inf],ppform.eval(0),ppform.dim);
            end
            obj = obj@Trajectory(ppform.dim);
            obj.pp = polyniminalTrajectory.minimalOrder(ppform);
            obj.tspan = [min(obj.pp.breaks) max(obj.pp.breaks)];
            
        end
    end
    
    methods (Static)
        function obj = loadobj(a)
            obj = polyniminalTrajectory(a.pp);
        end
    end
    
    methods 
        function y = eval(obj,t)
            t=max(min(t,obj.tspan(end)),obj.tspan(1));
            y = ppvalSafe(obj.pp,t);  % still benefits from being safe (e.g. for supporting TaylorVar)
        end
        
        function dtraj = fnder(obj,order)
            if nargin<2
                order = 1;
            end
            % this requires the curve-fitting toolbox, so i'm implementing it myself below
            %      dtraj = polyniminalTrajectory(fnder(obj.pp,order));
            
            % first handle order=1
            [b,c,l,k,d] = unmkpp(obj.pp);
            if (order>k-1)  % handle the case of too-high order
                dtraj = polyniminalTrajectory(mkpp(b,0*c(:,1),d));
                return;
            end
            
            for i=1:k-1
                cnew(:,i) = (k-i)*c(:,i);
            end
            dtraj = polyniminalTrajectory(mkpp(b,cnew,d));
            
            if (order>1)
                dtraj = fnder(dtraj,order-1);
            end
        end
        
        function df = deriv(obj,t)
            [b,c,~,k,d] = unmkpp(obj.pp);
            if (k==1)  % handle the case of too-high order
                ppform = mkpp(b,0*c(:,1),d);
            else
                for i=1:k-1
                    cnew(:,i) = (k-i)*c(:,i);
                end
                ppform = mkpp(b,cnew,d);
            end
            df = ppval(ppform,t);
        end  
        
        function nobj = shiftTime(obj,offset)
            typecheck(offset,'double');
            sizecheck(offset,[1 1]);
            obj.tspan = obj.tspan + offset;
            obj.pp.breaks = obj.pp.breaks + offset;
            nobj = polyniminalTrajectory(obj.pp);
        end
        
        function nobj = cutTime(obj,Tspan)
            typecheck(Tspan,'double');
            sizecheck(Tspan,[1 2]);
            if Tspan(1)<obj.tspan(1)||Tspan(2)>obj.tspan(2)
                error('wrong time interval')
            end
            obj.tspan = Tspan;
            breaks = [Tspan(1):1:Tspan(2)];
            state = obj.eval(breaks);
            fun = spline(breaks,state); 
            nobj = polyniminalTrajectory(fun);
        end
        
        function nobj = scaleTime(obj,scale)
            % Scale the time of the polyniminalTrajectory uniformly. The resulting trajectory has
            % the duration of each segment multiplied by [scale]. The coefficients are
            % also scaled in order to preserve the shape of the trajectory.
            [breaks, flat_coefs, l, k, d] = unmkpp(obj.pp);
            
            breaks = breaks * scale;
            if d == 1
                coefs = flat_coefs;
                for i = 1:size(coefs, 1)
                    coefs(i,:) = coefs(i,:) .* bsxfun(@power, (1/scale), (size(coefs, 2)-1):-1:0);
                end
            else
                % Even though they give us l, k, d in that order, coefs should actually be a d-by-l-by-k array.
                coefs = reshape(flat_coefs, [prod(d), l, k]);
                for j = 1:(prod(d))
                    for i = 1:size(coefs, 2)
                        coefs(j,i,:) = coefs(j,i,:) .* reshape(bsxfun(@power, (1/scale), (size(coefs, 3)-1):-1:0), 1, 1, []);
                    end
                end
                coefs = reshape(coefs, [d, l, k]);
            end
            if d == 1
                nobj = polyniminalTrajectory(mkpp(breaks, coefs));
            else
                nobj = polyniminalTrajectory(mkpp(breaks, coefs, d));
            end
        end
        
        function nobj = uminus(obj)
            obj.pp.coefs = -obj.pp.coefs;
            nobj = polyniminalTrajectory(obj.pp);
        end
        
        function tf = eq(a,b)
            % only implement the trivial case of pptraj=scalar const
            % (which is what I need right now)
            if isscalar(a)
                tmp=b;b=a;a=tmp;
            end
            if isscalar(b)
                % first check if it's a constant
                if any(a.pp.coefs(:,1:end-1)~=0), tf=false; return; end
                tf = all(a.pp.coefs(:,end)==b);
            else
                error('not implemented yet');
            end
        end
        
        function tf = ne(a,b)
            tf = ~eq(a,b);
        end
        
        function t = getBreaks(obj)
            t = obj.pp.breaks;
        end
        
        function [dmin,tmin]=approxNearestPoint(obj,xc)
            %compute the closest point to xc on the trajectory
            if(length(obj.getBreaks())==1)
                dmin=norm(xc-obj.eval(obj.getBreaks()));
                tmin=obj.getBreaks();
            else
                if(length(xc)~=obj.dim) error('Incorrect state dimension'); end
                tBreaks=obj.getBreaks;
                xdiff=repmat(xc,1,length(tBreaks))-obj.eval(tBreaks());
                [dminBreak,tminBreakIndex]=min(sum(xdiff.*xdiff,1));
                if(tminBreakIndex==1)
                    a=1;b=2;
                elseif(tminBreakIndex==length(tBreaks))
                    a=tminBreakIndex-1;
                    b=tminBreakIndex;
                else
                    a=tminBreakIndex-1;
                    b=tminBreakIndex;
                    c=tminBreakIndex+1;
                end
                coefs=obj.pp.coefs((a-1)*obj.dim+(1:obj.dim),:);
                coefs(:,end)=coefs(:,end)-xc;
                dtcoefs=repmat(fliplr(0:(obj.pp.order-1)),obj.dim,1).*coefs;
                dtcoefs=dtcoefs(:,1:(end-1));
                dtPoly=zeros(1,2*obj.pp.order-2);
                for i=1:obj.dim
                    dtPoly=conv(dtcoefs(i,:),coefs(i,:))+dtPoly;
                end
                dtPolyRoots=roots(dtPoly);
                dtPolyRoots=dtPolyRoots(isreal(dtPolyRoots)&dtPolyRoots>=0&dtPolyRoots<=(tBreaks(b)-tBreaks(a)));
                if(isempty(dtPolyRoots))
                    tminCandidate=[tBreaks(a) tBreaks(b)];
                else
                    tminCandidate=[tBreaks(a) tBreaks(a)+dtPolyRoots' tBreaks(b)];
                end
                if(tminBreakIndex~=1&&tminBreakIndex~=length(tBreaks))
                    coefs=obj.pp.coefs((b-1)*obj.dim+(1:obj.dim),:);
                    coefs(:,end)=coefs(:,end)-xc;
                    dtcoefs=repmat(fliplr(0:(obj.pp.order-1)),obj.dim,1).*coefs;
                    dtcoefs=dtcoefs(:,1:(end-1));
                    dtPoly=zeros(1,2*obj.pp.order-2);
                    for i=1:obj.dim
                        dtPoly=conv(dtcoefs(i,:),coefs(i,:))+dtPoly;
                    end
                    dtPolyRoots=roots(dtPoly);
                    dtPolyRoots=dtPolyRoots(isreal(dtPolyRoots)&dtPolyRoots>=0&dtPolyRoots<=(tBreaks(c)-tBreaks(b)));
                    if(isempty(dtPolyRoots))
                        tminCandidate=[tminCandidate tBreaks(c)];
                    else
                        tminCandidate=[tminCandidate tBreaks(b)+dtPolyRoots' tBreaks(c)];
                    end
                end
                xdiffCandidate=obj.eval(tminCandidate)-repmat(xc,1,length(tminCandidate));
                [dmin,dminIndex]=min(sum(xdiffCandidate.*xdiffCandidate,1));
                dmin=sqrt(dmin);
                tmin=tminCandidate(dminIndex);
            end
        end
        
        function traj = ctranspose(traj)
            [breaks,coefs,l,k,d] = unmkpp(traj.pp);
            if length(d)<2
                d = [1 d];
            elseif length(d)>2
                error('ctranspose is not defined for ND arrays');
            else
                coefs = reshape(coefs,[d,l,k]);
                coefs = permute(coefs,[2 1 3 4]);
                d=[d(end),d(1:end-1)];
            end
            traj = polyniminalTrajectory(mkpp(breaks,coefs,d));
        end 
        
        function c = plus(a,b)
            if ~isequal(size(a),size(b))
                error('must be the same size');  % should support scalars, too (but don't yet)
            end
            if any(size(a)==0)  % handle the empty case
                c = ConstantTrajectory(zeros(size(a)));
                return;
            end
            if isa(a,'ConstantTrajectory') a=double(a); end
            if isa(b,'ConstantTrajectory') b=double(b); end
            
            if isnumeric(a)  % then only b is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(b.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('plus is not defined for ND arrays'); end
                coefs = reshape(coefs,[d,l,k]);
                for i=1:l,
                    coefs(:,:,i,end)=a+coefs(:,:,i,end);
                end
                c=polyniminalTrajectory(mkpp(breaks,coefs,[size(a,1) d(2)]));
                return;
            elseif isnumeric(b) % then only a is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('plus is not defined for ND arrays'); end
                coefs = reshape(coefs,[d,l,k]);
                for i=1:l,
                    coefs(:,:,i,end)=coefs(:,:,i,end)+b;
                end
                c=polyniminalTrajectory(mkpp(breaks,coefs,[d(1) size(b,2)]));
                return;
            end
            
            if ~isa(a,'polyniminalTrajectory') || ~isa(b,'polyniminalTrajectory')
                % kick out to general case if they're not both pp trajectories
                c = plus@Trajectory(a,b);
                return;
            end
            
            [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
            [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp);
            
            if ~isequal(abreaks,bbreaks)
                breaks = unique([abreaks,bbreaks]);
                a = refine(a,breaks);
                b = refine(b,breaks);
                
                [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
                [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp);
            end
            if (ak>=bk)
                coefs=acoefs; coefs(:,end-bk+1:end)=coefs(:,end-bk+1:end)+bcoefs;
            else
                coefs=bcoefs; coefs(:,end-ak+1:end)=coefs(:,end-ak+1:end)+acoefs;
            end
            
            c = polyniminalTrajectory(mkpp(abreaks,coefs,ad));
        end
        
        function a = inv(a)
            if a.pp.order==1
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                coefs = reshape(coefs,[d,l]);  % k should be 1
                for i=1:l
                    coefs(:,:,i) = inv(coefs(:,:,i));
                end
                a = polyniminalTrajectory(mkpp(breaks,coefs,d));
            else
                a = inv@Trajectory(a);
            end
        end
        
        function c = minus(a,b)
            c = plus(a,uminus(b));
        end
        
        function c = mtimes(a,b)
            if any([size(a,1) size(b,2)]==0)  % handle the empty case
                c = ConstantTrajectory(zeros(size(a,1),size(b,2)));
                return;
            end
            if isa(a,'ConstantTrajectory') a=double(a); end
            if isa(b,'ConstantTrajectory') b=double(b); end
            
            if isnumeric(a)  % then only b is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(b.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('mtimes is not defined for ND arrays'); end
                if isscalar(a), cd = d; elseif isscalar(b), cd = size(a); else cd = [size(a,1),d(2)]; end
                coefs = full(reshape(coefs,[d,l,k])); a=full(a);
                for i=1:l, for j=1:k,
                        c(:,:,i,j)=a*coefs(:,:,i,j);
                    end, end
                c=polyniminalTrajectory(mkpp(breaks,c,cd));
                return;
            elseif isnumeric(b) % then only a is a polyniminalTrajectory
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                if length(d)<2, d=[d 1]; elseif length(d)>2, error('mtimes is not defined for ND arrays'); end
                if isscalar(a), cd = d; elseif isscalar(b), cd = size(a); else cd = [size(a,1),d(2)]; end
                coefs = full(reshape(coefs,[d,l,k])); b=full(b);
                for i=1:l, for j=1:k,
                        c(:,:,i,j)=coefs(:,:,i,j)*b;
                    end, end
                c=polyniminalTrajectory(mkpp(breaks,c,[d(1) size(b,2)]));
                return;
            end
            
            
            if ~isa(a,'polyniminalTrajectory') || ~isa(b,'polyniminalTrajectory')
                % kick out to general case if they're not both pp trajectories
                c = mtimes@Trajectory(a,b);
                return;
            end
            
            %      c = polyniminalTrajectory(fncmb(a.pp,'*',b.pp));  % this seems to fail on simple test cases??
            
            [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
            [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp);
            
            if ~isequal(abreaks,bbreaks)
                breaks = unique([abreaks,bbreaks]);
                a = refine(a,breaks);
                b = refine(b,breaks);
                
                [abreaks,acoefs,al,ak,ad] = unmkpp(a.pp);
                [bbreaks,bcoefs,bl,bk,bd] = unmkpp(b.pp); 
            end
            
            if (length(ad)<2) ad=[ad 1];
            elseif (length(ad)>2) error('mtimes not defined for ND arrays'); end
            if (length(bd)<2) bd=[bd 1];
            elseif (length(bd)>2) error('mtimes not defined for ND arrays'); end
            
            acoefs = reshape(acoefs,[ad,al,ak]);
            bcoefs = reshape(bcoefs,[bd,bl,bk]);
            
            %       ( sum a(:,:,j)(t-t0)^(k-j) ) ( sum b(:,:,j)(t-t0)^(k-j) )
            
            cbreaks = abreaks; % also bbreaks, by our assumption above
            if isscalar(a), cd = bd; elseif isscalar(b) cd = ad; else cd = [ad(1) bd(2)]; end
            cl = al;  % also bl, by our assumption that abreaks==bbreaks
            ck = ak+bk-1;
            
            ccoefs = zeros([cd,cl,ck]);
            for l=1:cl
                for j=1:ak  % note: could probably vectorize at least the inner loops
                    for k=1:bk
                        %            order_a = ak-j; order_b = bk-k;  order_c = order_a+order*b;
                        ccoefs(:,:,l,ck-(ak-j)-(bk-k))=ccoefs(:,:,l,ck-(ak-j)-(bk-k)) + acoefs(:,:,l,j)*bcoefs(:,:,l,k);
                    end
                end
            end
            c = polyniminalTrajectory(mkpp(cbreaks,ccoefs,cd));
        end
        
        function c = vertcat(a,varargin)
            % be careful for this by dengwei at 2017.12.4
            if isnumeric(a)
                a = ConstantTrajectory(a);
            end
            if isa(a,'ConstantTrajectory')
                breaks=[];
                coefs = a.pt;
                l=1; k=1; d=size(coefs);
            elseif isa(a,'polyniminalTrajectory')
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                d = size(a);  % handles case where d is a scalar, instead of [d,1]
                coefs = reshape(coefs,[d,l,k]);
            else
                c = vertcat@Trajectory(a,varargin{:});
                return;
            end
            for i=1:length(varargin)
                if isnumeric(varargin{i})
                    varargin{i}=ConstantTrajectory(varargin{i});
                end
                if isa(varargin{i},'ConstantTrajectory')
                    breaks2=[];
                    coefs2 = varargin{i}.pt;
                    l2=1; k2=1; d2=size(coefs2);
                elseif isa(varargin{i},'polyniminalTrajectory')
                    [breaks2,coefs2,l2,k2,d2]=unmkpp(varargin{i}.pp);
                    d2 = size(varargin{i}); % handles case where d2 is a scalar, instead of [d,1]
                else
                    c = vertcat@Trajectory(a,varargin{:});
                    return;
                end
                if ~isequal(d(2:end),d2(2:end))
                    error('incompatible dimensions');
                end
                if isempty(coefs2)
                    continue;
                elseif isempty(coefs)
                    breaks = breaks;
                    coefs = coefs2;
                    d=d2;k=k2;l=l2;
                    continue;
                end
                if isempty(breaks)&&isempty(breaks2)
                    % all constanttrajectories so far
                    coefs = vertcat(coefs,coefs2);
                    d=size(coefs);
                    continue;
                elseif isempty(breaks)
                    [breaks,coefs,l,k,d] = unmkpp(zoh(breaks2,repmat(coefs,[d*0+1,length(breaks2)])));
                    coefs = reshape(coefs,[d,l,k]);
                elseif isempty(breaks2)
                    [breaks2,coefs2,l2,k2,d2] = unmkpp(zoh(breaks,repmat(coefs2,[d2*0+1,length(breaks)])));
                    coefs2 = reshape(coefs2,[d2,l2,k2]);
                elseif ~isequal(breaks,breaks2)
                    warning('polyniminalTrajectory:DifferentBreaks','vertcat for pptrajectories with different breaks is not supported (yet).  kicking out to function handle version');
                    c = vertcat@Trajectory(a,varargin{:});
                    return;
                end
                if (k<k2)
                    coefs = reshape(coefs,prod(d)*l,k);
                    coefs = [zeros(prod(d)*l,k2-k),coefs];  % pad with zeros
                    k=k2;
                    coefs = reshape(coefs,[d,l,k]);
                elseif (k2<k)
                    coefs2 = reshape(coefs2,prod(d2)*l2,k2);
                    coefs2 = [zeros(prod(d2)*l2,k-k2),coefs2];  % pad with zeros
                    k2=k;
                    coefs2 = reshape(coefs2,[d2,l2,k2]);
                end
                coefs2=reshape(coefs2,[d2,l2,k2]);
                d = [d(1)+d2(1),d(2:end)];
                coefs = [coefs;coefs2];
            end
            if numel(d)==2 && d(2)==1
                d = d(1); % column vectors are a special case that's handled differently by the spline class
                coefs = reshape(coefs, [d, l, k]);
            end
            c = polyniminalTrajectory(mkpp(breaks,coefs,d)); 
        end
        
        function varargout = subsref(a,s)
            if (length(s)==1 && strcmp(s(1).type,'()'))
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                coefs = reshape(coefs,[d,l,k]);
                if (length(s.subs)==1 && length(d)>1)
                    subs = cell(1,length(d));
                    [subs{:}] = ind2sub(d,s.subs{:});
                    s.subs = subs;
                end
                s.subs = {s.subs{:},':',':'};
                coefs = subsref(coefs,s);
                d=size(subsref(a.eval(a.tspan(1)),s));
                if numel(d)==2 && d(2)==1, d = d(1); end  % column vectors are a special case that's handled differently by the spline class
                varargout{1} = polyniminalTrajectory(mkpp(breaks,coefs,d));
            else % use builtin
                varargout=cell(1,max(nargout,1));
                [varargout{:}] = builtin('subsref',a,s);
            end
        end
        
        function a = subsasgn(a,s,b)
            if (length(s)==1 && strcmp(s(1).type,'()'))
                if isempty(a) % handle the special case
                    [breaks,coefs,l,k,d] = unmkpp(b.pp);
                    e=[];
                    d_extended = [d,l,k];
                    coefs = reshape(coefs,d_extended);
                    s.subs = {s.subs{:},':',':'};
                    e = subsasgn(e,s,coefs);
                    a = polyniminalTrajectory(mkpp(breaks,e,d_extended(1:end-2)));
                    return;
                end
                if isnumeric(a) % then b must be a polyniminalTrajectory
                    breaks = b.getBreaks();
                    a = polyniminalTrajectory(zoh(breaks,repmat(a,[1+0*size(a),length(breaks)])));
                elseif isa(a,'ConstantTrajectory')
                    breaks = b.getBreaks();
                    a = polyniminalTrajectory(zoh(breaks,repmat(a.pt,[1+0*size(a),length(breaks)])));
                end
                typecheck(a,'polyniminalTrajectory');  % i believe this is the only way this method would get called
                [breaks,coefs,l,k,d] = unmkpp(a.pp);
                if isnumeric(b)
                    b = polyniminalTrajectory(zoh(breaks,repmat(b,[1+d*0,length(breaks)])));
                elseif isa(b,'ConstantTrajectory')
                    b = polyniminalTrajectory(zoh(breaks,repmat(b.pt,[1+d*0,length(breaks)])));
                end
                typecheck(b,'polyniminalTrajectory');
                [breaks2,coefs2,l2,k2,d2] = unmkpp(b.pp);
                if ~isequal(breaks,breaks2)
                    a = subsasgn@Trajectory(a,s,b);
                end
                if (k<k2)
                    coefs =  [zeros(prod(d)*l,k2-k),coefs];
                    k=k2;
                elseif (k2<k)
                    coefs2 = [zeros(prod(d2)*l2,k-k2),coefs2];  % pad with zeros
                    k2=k;
                end
                coefs = reshape(coefs,[d,l,k]);
                coefs2 = reshape(coefs2,[d2,l2,k2]);
                s.subs = {s.subs{:},':',':'};
                coefs = subsasgn(coefs,s,coefs2);
                d = size(coefs); d=d(1:end-2);
                a = polyniminalTrajectory(mkpp(breaks,coefs,d));
            else
                a = subsasgn@Trajectory(a,s,b);
            end
        end
        
        function newtraj = append(obj, trajAtEnd)
            % Append a polyniminalTrajectory to this one, creating a new trajectory that
            % starts where this object starts and ends where the given trajectory
            % ends.
            %
            % This will throw an error if the trajectory to append does not start
            % where the first trajectory ends.  This is useful if you did a bunch
            % of peicewise simulation and now want to combine them into one
            % object.
            %
            % @param trajAtEnd trajectory to append
            % @retval newtraj new polyniminalTrajectory object that is the combination of
            % both trajectories
            
            if ~isa(obj,'polyniminalTrajectory')
                obj = polyniminalTrajectory(obj);
                obj.tspan = [-10^8,trajAtEnd.tspan(1)];
                obj.pp.breaks = obj.tspan;
            end
            if ~isa(trajAtEnd,'polyniminalTrajectory')
                trajAtEnd = polyniminalTrajectory(trajAtEnd);
                trajAtEnd.tspan = [obj.tspan(2),10^8];
                obj.pp.breaks = obj.tspan;
            end
            
            % check for time condition
            firstEnd = obj.pp.breaks(end);
            secondStart = trajAtEnd.pp.breaks(1);
            
            if (firstEnd ~= secondStart)
                error(strcat('Cannot append trajectories that do not start/end at the same time.', ...
                    'First trajectory ends at t = ',num2str(firstEnd), ...
                    ' but the second trajectory starts at t = ', num2str(secondStart)));
            end
            
            
            % check for the same dimensions
            if (obj.pp.dim ~= trajAtEnd.pp.dim)
                error(strcat('Cannot append trajectories with different dimensionality.', ...
                    'First trajectory has pp.dim = ', num2str(obj.pp.dim), ...
                    ' but the second trajectory has pp.dim = ', num2str(trajAtEnd.pp.dim)));
            end
            
            % check for the same dimensions
            if (obj.dim ~= trajAtEnd.dim)
                error(strcat('Cannot append trajectories with different dimensionality.', ...
                    'First trajectory has dim = ', num2str(obj.pp.dim), ...
                    ' but the second trajectory has dim = ', num2str(trajAtEnd.pp.dim)));
            end
            
            % check for the same order
            if (obj.pp.order < trajAtEnd.pp.order)
                [breaks,coefs,l,k,d] = unmkpp(obj.pp);
                coefs = reshape(coefs,prod(d)*l,k);
                coefs = [zeros(prod(d)*l,trajAtEnd.pp.order-k),coefs];  % pad with zeros
                k=trajAtEnd.pp.order;
                obj.pp = mkpp(obj.pp.breaks,reshape(coefs,[d,l,k]),d);
            elseif (obj.pp.order > trajAtEnd.pp.order)
                [breaks,coefs,l,k,d] = unmkpp(trajAtEnd.pp);
                coefs = reshape(coefs,prod(d)*l,k);
                coefs = [zeros(prod(d)*l,obj.pp.order-k),coefs];  % pad with zeros
                k=obj.pp.order;
                trajAtEnd.pp = mkpp(trajAtEnd.pp.breaks,reshape(coefs,[d,l,k]),d);
            end
            
            newtraj = obj;
            
            newtraj.pp.dim = obj.pp.dim;
            newtraj.pp.order = obj.pp.order; 
            newtraj.pp.pieces = obj.pp.pieces + trajAtEnd.pp.pieces; 
            newtraj.pp.breaks = [obj.pp.breaks trajAtEnd.pp.breaks(2:end)];
            newtraj.pp.coefs = [obj.pp.coefs; trajAtEnd.pp.coefs]; 
            
            newtraj.dim = obj.dim;
            
            newtraj.tspan = [min(obj.pp.breaks) max(trajAtEnd.pp.breaks)];
        end % append
        
        function obj = refine(obj,newbreaks)
            obj = polyniminalTrajectory(pprfn(obj.pp,newbreaks));
        end
        
        function obj = rebreaks(obj,newbreaks)
            if newbreaks==obj.pp.breaks
                return;
            end
            newtraj = ppvalSafe(obj.pp,newbreaks);
            newpp = spline(newbreaks,newtraj);
            obj = polyniminalTrajectory(newpp);
        end
        
        function h = plot(obj,varargin)
%             ts = linspace(obj.tspan(1),obj.tspan(end),100);
            % modeify by dengwei at 2018.1.31
            % modeify by dengwei at 2018.5.17
            ts = obj.pp.breaks;
            xtraj = obj.eval(ts);
            if nargin>1
                clr = varargin{1};
            else
                clr = 'k';
            end
            gca;
            if nargin>2
                plotdim = varargin{2}; 
                if length(plotdim)==3
                    h = plot3(xtraj(plotdim(1),:),xtraj(plotdim(2),:),xtraj(plotdim(3),:),clr);
                elseif length(plotdim)==2
                    h = plot(xtraj(plotdim(1),:),xtraj(plotdim(2),:),clr);
                elseif length(plotdim)==1
                    h = plot(ts,xtraj(plotdim(1),:),clr);
                else
                    error('no such a fucking dim')
                end
            else 
                if obj.dim==3
                    h = plot3(xtraj(1,:),xtraj(2,:),xtraj(3,:),clr);
                elseif obj.dim==2
                    h = plot(xtraj(1,:),xtraj(2,:),clr);
                elseif obj.dim==1
                    h = plot(ts,xtraj,clr);
                else
                    error('no such a fucking dim')
                end
            end 
        end
        
    end
    
    methods (Static=true)
        function ppform = minimalOrder(ppform)
            if ppform.order>1 && all(ppform.coefs(:,1)==0)
                ppform.coefs = ppform.coefs(:,2:end);
                ppform.order = ppform.order-1;
                ppform = polyniminalTrajectory.minimalOrder(ppform); % recurse
            end
        end
    end
end
