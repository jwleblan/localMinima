% sfigure - Silent figure command
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   see figure.m
%
% PURPOSE:
%   This function implements the interface of figure.m but does not capture
%   focus.  The "s" is for silent.
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
% Copyright (C) 2011 Michigan Technological University
%-------------------------------------------------------------------------------
%}
function varargout = sfigure(varargin)
if(nargin)
        if ishandle(varargin{1})
            h = varargin{1};
            set(0, 'CurrentFigure', h);
            if(nargin>1)
                set(h,varargin{2:end});
            end
        elseif isnumeric(varargin{1})
            % We don't escape the focus here yet
            if ~mod(varargin{1},1)
                h = figure(varargin{1});
                if(nargin>1)
                    set(h,varargin{2:end});
                end
            else
                % figure not able to allocate non-integer handles yet
                h = figure(varargin{2:end});
            end
            
        else
            % We don't escape the focus here yet
            h = figure(varargin{:});

        end
        
else
    % We don't escape the focus here yet
    h = figure();
end

% Assign the handle
if(nargout>0)
    varargout{1} = h;
end

end