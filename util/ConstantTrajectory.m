classdef ConstantTrajectory<Trajectory
    % Trivial instance of a trajectory as a constant.
    % build a trajectory with a constant form
    % build by dengwei 2016
    properties
        pt  
    end
    
    methods
        function obj = ConstantTrajectory(pt)
            obj = obj@Trajectory(size(pt));
            obj.pt=pt;
            obj.tspan = [-inf,inf];
            obj.dim = size(pt);
        end
        
        function a = double(obj)
            a=double(obj.pt);
        end
        
        function dtraj = fnder(obj)
            % Implements the (trivial) derivative.
            dtraj = ConstantTrajectory(0*obj.pt);
        end
        
        function y = eval(obj,t)
            % Return the fixed point for all t.
            if isscalar(t)
                y=obj.pt;
            else
                d = size(obj.pt);
                dd = 0*d+1;  % 1 for all values
                if (d(end)==1),
                    dd(end) = length(t);
                else
                    dd(end+1) = length(t);
                end
                y = repmat(obj.pt,dd);
            end
        end
        
        function t = getBreaks(obj)
            % Return a single break (at t=0).
            t = 0;
        end
        
        function a = subsasgn(a,s,b)
            if isa(b,'polyniminalTrajectory')
                breaks = getBreaks(b);
                a = polyniminalTrajectory(zoh(breaks,repmat(a.pt,1,length(breaks))));
                a = subsasgn(a,s,b);
                return;
            end
            if isa(a,'ConstantTrajectory')
                if isnumeric(b)
                    a = subsasgn(a.pt,s,b);
                elseif isa(b,'ConstantTrajectory')
                    a = subsasgn(a.pt,s,b.pt);
                else
                    a = subsasgn@Trajectory(a,s,b);
                end
            else % b must be a ConstantTrajectory
                a = subsasgn(a,s,b.pt);
            end
        end
        
        function obj = uminus(obj)
            obj.pt = -obj.pt;
        end
        
        function a = ctranspose(a)
            a = ConstantTrajectory(a.pt');
        end
        
        function c = mtimes(a,b)
            if any([size(a,1) size(b,2)]==0)  % handle the empty case
                c = sparse(zeros(size(a,1),size(b,2)));
                return;
            end
            if isa(a,'ConstantTrajectory')
                a=a.pt;
            elseif isa(a,'polyniminalTrajectory')
                c = mtimes@polyniminalTrajectory(a,b);
                return;
            end
            if isa(b,'ConstantTrajectory')
                b=b.pt;
            elseif isa(b,'polyniminalTrajectory')
                c = mtimes@polyniminalTrajectory(a,b);
                return;
            end
            if isnumeric(a)&&isnumeric(b)
                c=ConstantTrajectory(a*b);
            else
                error('do not define num times ConstantTrajectory');
            end
        end
        
        function c = plus(a,b)
            if isa(a,'ConstantTrajectory') a=a.pt;
            elseif ~isa(a,'numeric') c = plus@Trajectory(a,b); return; end
            if isa(b,'ConstantTrajectory') b=b.pt;
            elseif ~isa(b,'numeric') c = plus@Trajectory(a,b); return; end
            c=ConstantTrajectory(a+b);
        end
        
        function c = inv(a)
            c = ConstantTrajectory(inv(a.pt));
        end
        function ydot = deriv(obj,t)
            ydot = eval(fnder(obj),t);
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

    end
end
