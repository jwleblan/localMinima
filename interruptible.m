% interruptible - Allow interruption of executing code
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   status = interruptible()
%
% PURPOSE:
%   This function generates a window that allows for the interruption of matlab.
%   The best use of interruptible is to explicilty design interrupts into your
%   program, however, a break-point option "break" is also provided.
%  
% INPUT:
% 
% OUTPUT:
%   status  - True if an interrupt was requested and false otherwise
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
% Copyright (C) 2015 Joel W. LeBlanc
%-------------------------------------------------------------------------------
%}
function status = interruptible()
    
    % Default Values
    figData = figDataFun;
    status = false;
    
    % Create the figure if necessary
    if isempty(figData.h) || ~ishandle(figData.h)
        figData = createFigure(figData);
    end
    
    if figData.callPending
        % Reset the flag
        figData.callPending = false;
        figDataFun(figData);
        
        % Address the request
        if(nargout)
            status = true;
        else
            [db,workInd] = dbstack;
            dbstop('in',...
            db(workInd+1).file,...
            'at',...
            sprintf('%d',db(workInd+1).line+1));
        end
        
        % Reset visual notification
        set(figData.h,'Name','Interrupt...');
        drawnow;
    end
    
    % Flush the event queue
    % (This is the entry point for the callbacks below)
    drawnow;

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   figData = createFigure(figData)
%
% PURPOSE:
%   This function generates the figure based on figData, and updates that data
%   to be consistent with the newly created figure.
%  
% INPUT:
% 
% OUTPUT:
%   figData  - Structure describing the figure's state (see figDataFun)
%-------------------------------------------------------------------------------
%}
function figData = createFigure(figData)

tagStr = 'Interruptible Main Window';
figureOptions = {...
        'BusyAction','cancel',...
        'ButtonDownFcn',@cb_interruptFun,...
        'CloseRequestFcn',@cb_CloseRequestFcn,...
        'Color',[1 0 0],...
        'HandleVisibility','off',...
        'HitTest','off',...
        'IntegerHandle','off',...
        'Interruptible','off',...   % Ironic
        'MenuBar','none',...
        'Name','Interrupt...',...
        'NumberTitle','off',...
        'OuterPosition',figData.OuterPosition,...
        'Tag',tagStr};

% Make sure we're cleaning up after ourselves (good diagnostic too)
showHiddenHandlesState = get(0,'ShowHiddenHandles');
set(0,'showHiddenHandles','on');
obj = findobj('Type','figure','Tag',tagStr);
set(0,'showHiddenHandles',showHiddenHandlesState);
if numel(obj)
    warning('Unexpected Hanging Interruptible Window');
    delete(obj);
    figData.h = [];
end


% Generate a new figure
if ~isempty(figData.h) && ishandle(figData.h)
    sfigure(figData.h,figureOptions{:});
else
    figData.h = sfigure(figureOptions{:});
end

figDataFun(figData);
    
end

function cb_interruptFun(h,data)

[db,workInd] = dbstack;

% Make sure we don't remain the current object
set(0,'currentFigure',[]);

% If we are not being called by a function
if(numel(db)<3)
    return;
end

% If we got called at someone elses event queue flush
if(~strcmp(db(workInd+1).name,mfilename))
    figData = figDataFun();
    
    if(figData.callPending)  % User's freaking out a bit
        disp('Multiple Interrupts Recieved... Breaking Out.');
        
        % Reset visual notification
        set(figData.h,'Name','Interrupt...');
        drawnow;
        
        % Setup a break point if possible
        dbstop('in',...
            db(workInd+1).file,...
            'at',...
            sprintf('%d',db(workInd+1).line+1));
        
        % Clear the previous notification
        figData.callPending = false;
        
    else
        % Set a flag so we deal with this interrupt next time we're called
        figData.callPending = true;
        
        % Give the user some visual confirmation
        set(h,'Name','Pending...');
        drawnow;
    end
    figDataFun(figData);
    
else
    % If the user wants the status
    if ~evalin('caller','figData.callPending')
        evalin('caller','status=true;');
        
    else    % Generate the breakpoint
        dbstop('in',...
            db(workInd+2).file,...
            'at',...
            sprintf('%d',db(workInd+2).line+1));
    end
end

end

function cb_CloseRequestFcn(h,data)

% Record the last location
figData = figDataFun;
figData.OuterPosition = get(h,'OuterPosition');
figDataFun(figData)

% Delete the figure
delete(h);
end

function varargout = figDataFun(varargin)
persistent DATA
%mlock;              % This is needed because people get clear happy
if nargin
    DATA = varargin{1};
    
elseif isempty(DATA)
    % Initialization
    DATA = struct('h',[],'OuterPosition',[25 25 175 100],'callPending',false);
end

if nargout
	varargout{1} = DATA;
end
end
