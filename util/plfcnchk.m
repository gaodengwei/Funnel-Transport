function [outx,xexpr,inargs,emsg] = plfcnchk(inx,textnum,defarg)
% build by dengwei gao 2017.12.4
% trans char inline function to x. which can using feval
if (nargout>2), inargs = {}; end
if (nargout>3)
    emsg.message = '';
    emsg.identifier = '';
    emsg = emsg(zeros(0,1)); % make sure emsg is the right dimension
end

% Remove array notation from inline functions
xexpr = inx;
if (nargout>1) & isa(inx, 'inline')
    xexpr = char(inx);
    if (length(xexpr) > 2)
        xexpr(findstr(xexpr, '.*')) = [];
        xexpr(findstr(xexpr, '.^')) = [];
        xexpr(findstr(xexpr, './')) = [];
    end
end

% Check for function vs. identity expression such as 't',
% and convert to a consistent form
outx = inx;
numonly = 0;             % indicates a numeric expression, no variables
e = 0;
if (ischar(inx))
    numonly = isempty(symvar(inx));
    if (numonly & nargin>1 & textnum==1), return; end
    
    e = exist(inx);
    dovector = (e~=2 & e~=5) | (e==5 & numonly);
    
    % Always vectorize an inline rather than a character string,
    % because we don't want to add 'ones(size(x))' to the expression.
    if (dovector), outx = inline(inx); end
else
    dovector = isa(inx,'inline');
end
if (dovector)
    [outx,emsg] = fcnchk(outx, 'vectorized');
else
    [outx,emsg] = fcnchk(outx);
end

% Display error unless caller is going to deal with it
if (nargout<4), error(emsg); end

% Get function expression for label, if not already done
if (isa(outx,'function_handle')), xexpr = func2str(outx); end

% Sometimes we can check the number of inputs expected
if (nargout>2)
    if (isa(outx, 'inline'))
        inargs = argnames(outx);
        
        % Some inputs such as 'sin' will become inline functions with
        % an argument 'x' that does not really appear in the input.
        % In that case, give the argument an "anonymous" name.
        if (ischar(inx) & isequal(inargs,{'x'}) & ...
                ~ismember({'x'},symvar(inx)))
            inargs = {''};
        end
    elseif isa(outx,'function_handle')
        if ismember('(',xexpr) && ismember(')',xexpr)
            % For function handles that are anonymous functions, get the
            % names of the args and define a reasonable expression
            A = find(xexpr=='(',1);
            B = find(xexpr==')',1);
            inargs = strtrim(strread(xexpr(A+1:B-1),'%s','delimiter',','));
            xexpr = strtrim(xexpr(B+1:end));
        end
    elseif (ischar(outx) & ((e==2) | (e==5 & ~isempty(symvar(outx)))))
        if (e==2)
            k = nargin(outx);
        else
            k = 1;  % assume one argument for builtin function
        end
        if (k>=1), inargs{1} = ''; end
        if (k>=2), inargs{2} = ''; end
        if (nargin>2)
            xexpr = [xexpr '(' defarg ')'];
        elseif (k==2)
            xexpr = [xexpr '(x,y)'];
        else
            xexpr = [xexpr '(x)'];
        end
    end
end
end