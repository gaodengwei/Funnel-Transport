classdef Trajectory
  
  properties
    tspan
    ts=[];    % default sample times of the model
  end
  
  properties (SetAccess=protected,GetAccess=protected)
    dim   
  end
  
  methods (Abstract=true)
    y = eval(obj,t);
    t = getBreaks(obj); 
  end
 
  methods 
    function obj = Trajectory(dim) 
      obj.dim = dim;
    end
 
    function y = output(obj,t,x,u)
      if (isnumeric(obj.dim))
        y = reshape(obj.eval(t),[],1);
      else 
        y=nan;
      end
    end
 
    function ydot = deriv(obj,t)
      ydot = eval(fnder(obj),t);
    end
    
    function yddot = dderiv(obj,t)
      yddot = eval(fnder(fnder(obj)),t);
    end
 
    function mobj = uminus(obj)
      mobj = FunctionHandleTrajectory(@(t)-obj.eval(t),obj.dim,obj.getBreaks,@(t)-obj.deriv(t));
    end
 
    function [a,b,breaks] = setupTrajectoryPair(a,b)
        if (isnumeric(a)) a = ConstantTrajectory(a); end
        typecheck(a,'Trajectory'); 
        if (isnumeric(b)) b = ConstantTrajectory(b); end
        typecheck(b,'Trajectory');
        tspan = [max(a.tspan(1),b.tspan(1)),min(a.tspan(2),b.tspan(2))];
        if (tspan(2)<tspan(1))
          error('Trajectory:IncompatibleTimesForConcatenation','concatenated trajectories do not overlap in time');
        end
        
        breaks = unique([reshape(a.getBreaks,1,[]),reshape(b.getBreaks,1,[])]);
        breaks = breaks(breaks>=tspan(1) & breaks<=tspan(2));
        breaks = breaks';
    end
    
    function b = ctranspose(a)
      s = size(a);
      b = FunctionHandleTrajectory(@(t) ctranspose(a.eval(t)),s([end,1:end-1]),a.getBreaks);
    end
    
    function b = inv(a)
      s = size(a);
      b = FunctionHandleTrajectory(@(t) inv(a.eval(t)),s([end,1:end-1]),a.getBreaks);
    end
    
    function c = vertcat(a,varargin)
      c=a;
      for i=1:length(varargin)
        [c,b,breaks]=setupTrajectoryPair(c,varargin{i});
        sc = size(c); sb = size(b);
        if (length(sc)~=length(sb) || any(sc(2:end)~=sb(2:end)))
          error('dimensions 2:end must match');
        end 
        c = FunctionHandleTrajectory(@(t) vertcat(c.eval(t),b.eval(t)),[c.dim(1)+b.dim(1),c.dim(2:end)],breaks);
      end
    end
    
    function c = horzcat(a,varargin)
      c=a;
      for i=1:length(varargin)
        [c,b,breaks]=setupTrajectoryPair(c,varargin{i});
        cdim=size(c); bdim=size(b);
        if (length(cdim)~=length(bdim) || any(cdim([1,3:end])~=bdim([1,3:end])))
          error('dimensions 1 and 3:end must match');
        end
        c = FunctionHandleTrajectory(@(t) horzcat(c.eval(t),b.eval(t)),[cdim(1),cdim(2)+bdim(2),cdim(3:end)],breaks);
 
      end
    end
    
    function c = power(a,b)
      [a,b,breaks]=setupTrajectoryPair(a,b);
      c = FunctionHandleTrajectory(@(t) a.eval(t).^b.eval(t),max(size(a),size(b)),breaks);
    end
    
    function a = prod(x,dim)
      if (nargin<2)
        dim=find(size(x)>1,'first');
      end
      a = FunctionHandleTrajectory(@(t) prod(x.eval(t),dim), [x.dim(1:dim-1),x.dim(dim+1:end)],x.getBreaks);
    end
    
    function c = times(a,b)
      [a,b,breaks]=setupTrajectoryPair(a,b);
      c = FunctionHandleTrajectory(@(t) a.eval(t).*b.eval(t),size(a),breaks);
    end
    
    function c = mtimes(a,b)
      if (ndims(a)>2 || ndims(b)>2) error('only defined for two-d matrices'); end
      if (~isscalar(a) && ~isscalar(b) && size(a,2) ~= size(b,1)) error('dimension mismatch'); end

      [a,b,breaks]=setupTrajectoryPair(a,b);
      c = FunctionHandleTrajectory(@(t) a.eval(t)*b.eval(t),size(zeros(size(a))*zeros(size(b))),breaks);
    end
    
    function c = plus(a,b)
      if (ndims(a) ~= ndims(b)) error('dimension mismatch'); end
      if (~isscalar(a) && ~isscalar(b) && any(size(a)~=size(b))) error('dimension mismatch'); end
      
      [a,b,breaks]=setupTrajectoryPair(a,b);
      c = FunctionHandleTrajectory(@(t) a.eval(t)+b.eval(t),max(size(a),size(b)),breaks);
    end
    
    function c = minus(a,b)
      if (ndims(a) ~= ndims(b)) error('dimension mismatch'); end
      if (~isscalar(a) && ~isscalar(b) && any(size(a)~=size(b))) error('dimension mismatch'); end
      
      [a,b,breaks]=setupTrajectoryPair(a,b);
      c = FunctionHandleTrajectory(@(t) a.eval(t)-b.eval(t),max(size(a),size(b)),breaks);
    end
    
    function a = subsasgn(a,s,b)
      [a,b,breaks]=setupTrajectoryPair(a,b);      
      a = FunctionHandleTrajectory(@(t) subsasgn(a.eval(t),s,b.eval(t)),size(subsasgn(a.eval(breaks(1)),s,b.eval(breaks(1)))),breaks);
    end

    function n=numel(varargin)
      % have to set numel to 1, because calls to subsref 
      % with . and {} automatically ask for numel outputs.
      % note that I *hate* that numel ~= prod(size(obj))
      n=1;
    end
    
    function varargout = subsref(a,s)
      if (length(s)==1 && strcmp(s(1).type,'()'))
        breaks = a.getBreaks();
        varargout=cell(1,max(nargout,1));
        [varargout{:}] = FunctionHandleTrajectory(@(t) subsref(a.eval(t),s),size(subsref(a.eval(breaks(1)),s)),breaks);
      else  % use builtin
        varargout=cell(1,max(nargout,1));
        [varargout{:}] = builtin('subsref',a,s);
      end
    end
 
    function s = size(obj,dim)
      s=obj.dim;
      if (length(s)==1) s=[s,1]; end
      if (nargin>1) s=s(dim); end
    end
    
    function l = length(obj)
      s = size(obj);
      l = max(s);
    end
 
    function ts = getTimeSpan(obj)
      ts = obj.tspan;
    end
 
    
    function traj = trim(obj,newtspan)
      sizecheck(newtspan,[1 2]);
      assert(newtspan(1)>=obj.tspan(1));
      assert(newtspan(2)<=obj.tspan(2));
      breaks = getBreaks(obj);
      breaks(breaks<newtspan(1))=[];
      breaks(breaks>newtspan(2))=[];
      breaks = unique([newtspan(1),breaks,newtspan(2)]);
      traj = FunctionHandleTrajectory(@(t) eval(obj,t), obj.dim, breaks);
    end
    
     
    
    function h=fnplt(obj,plotdims)
      if (nargin>1 && ~isempty(plotdims) && any(plotdims>prod(obj.dim) | plotdims<1)) error('plotdims out of range'); end
      breaks=obj.getBreaks();
      if iscolumn(breaks) breaks = breaks'; end 
      m=5; t=linspace(0,1,m)'; n=length(breaks)-1;
      ts = repmat(1-t,1,n).*repmat(breaks(1:end-1),m,1) + repmat(t,1,n).*repmat(breaks(2:end),m,1);
      ts = ts(:);
      pts = obj.eval(ts);
      if (prod(obj.dim)==1)
        h=plot(ts,squeeze(pts),'b.-','LineWidth',1,'MarkerSize',5);
        xlabel('t');
        ylabel('x');
      elseif (nargin>1 && ~isempty(plotdims) && length(plotdims)==1)
        h=plot(ts,squeeze(pts(plotdims,:)),'b.-','LineWidth',1,'MarkerSize',5);
        xlabel('t'); 
        ylabel('x');
      else
        if (nargin<2 || isempty(plotdims)) plotdims=[1,2]; end
        h=plot(pts(plotdims(1),:),pts(plotdims(2),:),'b-',pts(plotdims(1),[1:m:end,end]),pts(plotdims(2),[1:m:end,end]),'b.','LineWidth',1);% ,'MarkerSize',5);
 
        xlabel('t');
        ylabel('x');
      end
    end
    
    function [xnear,tnear,dxnear,dtnear] = closestPoint(obj,x,knot_ind)
      % returns the closest point on the trajectory to the sample point x
      %
      % @param x the sample point
      % @param knot_ind (optional) when available, this indicates the index
      % of the closest knot point.  
      % 
      % @retval xnear the closest point on the trajectory
      % @retval tnear the time on the trajectory associated with xnear
      
      t=obj.getBreaks();
      if (nargin<3)
        xt = obj.eval(t);
        xbar = x(:,ones(1,length(t))) - xt;
        d = sum(xbar.^2,1);
        [dmin,knot_ind] = min(d);
      else
        xt = obj.eval(t(knot_ind));
        xbar = x-xt;
        dmin = sum(xbar.^2,1);
      end
       
      tspan = [t(max(knot_ind-1,1)),t(min(knot_ind+1,length(t)))];
      options = optimset('fminbnd');
      options.TolX = 1e-15;
      [tnear,dmin,exitflag]=fminbnd(@(t) sum((x-obj.eval(t)).^2,1),tspan(1),tspan(end),options);
      
      if (exitflag<1) 
        warning('fminbnd exited with exitflag %d',exitflag);
      end
      xnear = obj.eval(tnear);
      
      if (nargout>2)
        xdmin = obj.deriv(tnear);
        
        if (min(abs(tnear-t([1,end])))<=1e-14) % then I'm at one of the rails (value should be larger than options.TolX above)
          dtnear = zeros(1,obj.num_y);
        else
          xddmin = obj.dderiv(tnear);
          dtnear = xdmin'/(xdmin'*xdmin + (xnear-x)'*xddmin);
        end
        dxnear = xdmin*dtnear;
      end
    end
    
    function d = distance(obj,x)
      xnear = closestPoint(obj,x);
      d = norm(xnear-x);
    end
    
    function [d,dd] = distanceSq(obj,x)
      if (nargout>1)
        [xnear,~,dxnear,~] = closestPoint(obj,x);
        d = (xnear-x)'*(xnear-x);
        dd = 2*(xnear-x)'*(dxnear - eye(obj.num_y));
      else
        xnear = closestPoint(obj,x);
        d = (xnear-x)'*(xnear-x);
      end
    end
    
    % added at 2018 by dengwei 
    function ts = getSampleTime(obj)
      % As described at http://www.mathworks.com/help/toolbox/simulink/sfg/f6-58760.html
      % to set multiple sample times, specify one *column* for each sample
      % time/offset pair.
      % The default behavior is continuous time for systems with only continuous
      % states, and discrete time (with sample period 1s) for systems only
      % discrete states, and inherited for systems with no states.  For
      % systems with both discrete and continuous states, an error is
      % thrown saying that this function should be overloaded to set the
      % desired behavior.
      if ~isfield(obj,'ts')
            ts = [0;0];  % continuous time// default
      elseif ~isempty(obj.ts)
        ts = obj.ts; 
%       elseif (obj.num_xc>0 && obj.num_xd==0)
%         ts = [0;0];  % continuous time, no offset
%       elseif (obj.num_xc==0 && obj.num_xd>0)
%         ts = [1;0];  % discrete time, with period 1s.
%       elseif (obj.num_xc==0 && obj.num_xd==0)
%         ts = [-1;0]; % inherited sample time
      else
        error('systems with both discrete and continuous states must implement the getSampleTime method or call setSampleTime to specify the desired behavior');
      end
    end
    
    function obj = setSampleTime(obj,ts)
      % robust method for setting default sample time
      %
      % @param ts a 2-by-n matrix with each column containing a sample time
      %    redundant colums are eliminated automatically.
      % only a few possibilities are allowed/supported
      %   inherited, single continuous, single discrete, single continuous+single
      %   discrete (note: disabled single continuous + single discrete
      %   because it wasn't obviously the right thing... e.g. in the
      %   visualizer who asked for the output to be at fixed dt, but after
      %   combination, the output gets called continuously).

      ts = unique(ts','rows')';

      if size(ts,2)>1  % if only one ts, then all is well
        if any(ts(1,:)==-1)  % zap superfluous inherited
          ts=ts(:,ts(1,:)~=-1);
        end
        if sum(ts(1,:)>0)>1 % then multiple discrete
          error('cannot define trajectory have different discrete sample times');
        end
        if sum(ts(1,:)==0)>1 % then multiple continuous
          error( 'continuous time'' and ''continuous time, fixed in minor offset'' sample times');
        end
        if sum(ts(1,:)>=0)>1 % then both continuous and discrete
          error( 'cannot define trajectory that have both continuous and discrete sample times');
        end
      end
      obj.ts = ts;
    end
  end
end
