% processPrecond - Process preconditioners
%{   
% SYNTAX:
%   [A,AInv,b,<x0>] = processPrecond(A,AInv,b,<maxSVDSize>,<x0,precondX0>)
%
% DESCRIPTION: 
%   Breaks down the preconditioner.  This function will strip ill-conditioned
%   dimensions, and properly handles diagonal (vector) preconditioners.  If x0
%   is passed, it is moved into the preconditioned space.
%
%   If AInv is passed but A is not, then A will not be formed.  When AInv is a
%   vector it gets processed, and otherwise it is left alone.  This is because
%   it is faster to solve the linear systems than perform the decomposition in
%   most cases, and because we're not going to perform a O(n^3) calculation just
%   to check your math.  The upshot is that if you form both matricies we'll
%   assume you know best.  If you give us vectors we'll check everything because
%   we can quickly do so, and for almost all A's (all but those violating
%   maxSVDDim) we need to break them down, so we'll check your math while we're
%   at it.
%
% INPUTS:
%   A           - Preconditioner as a vector or matrix.
%   AInv        - Inverse preconditioner as a vector or matrix.
%   b           - Preconditioner offset or []
%   maxSVDSize  - When a decomposition is needed and the minimum dimension size
%                 exceeds this number it will not be performed and direct system
%                 solving will be used instead.
%                 Default: 3000
%   x0          - Initial vector
%   precondX0   - True if x0 is in the preconditioned space and false otherwise
%
% OUTPUT:
%   Numerically stable preconditioners and x0 in the preconditioned space
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order.
%
% Copyright (C) 2013-2019 Joel W. LeBlanc
%}
function [A,AInv,b,x0] = processPrecond(A,AInv,b,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parse Input Arguments %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxSVDSize = 3000;
% We want to bias K to contain AInv
if(isempty(AInv))
    AInvEmpty = true;
    if(isempty(A))  % Just being nice and checking this
        error('One of A or AInv must be passed');
    end
    K = A;
else
    AInvEmpty = false;
    K = AInv;
end
KSize = size(K);

if(nargin < 5)
    classK = class(K);
    numX0 = 0;
    x0 = [];
    
    if(nargin==4) % Got maxSVDDim
        maxSVDSize = varargin{1};
    end
else
    if(nargin==6) % Got maxSVDDim
        maxSVDSize = varargin{1};
        x0 = varargin{2};
        precondX0 = varargin{3};
    else
        x0 = varargin{1};
        precondX0 = varargin{2};
    end
    classK = class(x0);
    numX0 = numel(x0);
end
tol = max(KSize) * eps(cast(max(K(:)),classK));

if(KSize(2)==1)     % isvector(K)
    %ind = find(abs(K)<tol);
    ind = find(K==0);
    clear('K');     % Avoid a lazy copy down the road
    numInd = numel(ind);
   
    if(numInd)
        ASize = [KSize(1)-numInd,KSize(1)];
        temp = 1:ASize(2); temp(ind) = [];
        
        if(AInvEmpty)
            % A exists and needs to be expanded into a sparse matrix.  Ind
            % points to dimensions that form a null-space.  Reforming AInv
            % is computationally efficient so we'll just do it.  This also
            % ensures a matched A, AInv pair
            A(ind) = [];
            AInv = sparse(temp,1:ASize(1),1./A,ASize(2),ASize(1),ASize(1));
            A = sparse(1:ASize(1),temp,A,ASize(1),ASize(2),ASize(1));
        else
            % AInv exists and needs to be expanded into a sparse matrix. Ind
            % points to dimensions that form a null-space
            AInv(ind) = [];
            A = sparse(1:ASize(1),temp,1./AInv,ASize(1),ASize(2),ASize(1));
            AInv = sparse(temp,1:ASize(1),AInv,ASize(2),ASize(1),ASize(1));
        end
        
        % At this point A and AInv are populated and ASize points to the full
        % matrix size of these matricies even if they remain vectors
        if(numX0)
            if(precondX0)
                if(numX0==ASize(2))
                    % We neeed to eliminate a dimension provided by x0
                    if(isempty(b))
                        b = AInv*x0(temp);
                        b(ind) = x0(ind);
                        x0 = zeros(ASize(1),1);
                    end
                    warning('%d dimension(s) of x0 in a preconditioner null-space have been removed',...
                        numInd);

                elseif(numX0==ASize(1))
                    % It is REALLY unlikely this was intended... Warn strongly
                    warning('x0 seemed to account for preconditioner null-space.  This is very likely a mistake.');
                
                else
                    error('Unexpected length of x0');
                end
                
            else % X0 is not in the preconditioned space
                if(numX0==ASize(2))
                    % We neeed to eliminate a dimension provided by x0
                    if(isempty(b))
                        b = x0;
                        x0 = zeros(ASize(1),1);
                    else
                        b(ind) = x0(ind);   % Lock down the dimensions to be eliminated
                        x0 = A*(x0-b);
                    end
                    warning('processPrecond:dimRemoved',...
                        '%d dimension(s) of x0 in a preconditioner null-space have been removed',...
                        numInd);
                else
                    error('Unexpected length of x0');
                end
            end
        end
        
    else
        if(AInvEmpty)
            AInv = 1./A;
        else
            A = 1./AInv;
        end
        if(numX0)
            if(~precondX0 && numX0==KSize(1))
                if(isempty(b))
                    x0 = A.*x0;
                else
                    x0 = A.*(x0-b);
                end
            elseif(numX0~=KSize(1))
                error('Unexpected length of x0');
            end
        end
    end
else
    if(AInvEmpty)
        if(min(KSize)<=maxSVDSize)  % We're OK to SVD
            if issparse(A)
                [U,S,V] = svds(A,min(KSize));
            else
                [U,S,V] = svd(A,'econ');
            end
            s = diag(S);
            ind = find(s<tol,1,'first');
            
            if(isempty(ind))
                S = sparse(1:KSize(1),1:KSize(1),s);
                A = S*V';
                AInv = V*inv(S);
                if(numX0)
                    if(precondX0 && numX0==KSize(1))
                        x0 = U'*x0;
                    elseif(~precondX0 && numX0==KSize(2))
                        % Some part of x0 may lay in the null-space
                        if(isempty(b))
                            b = x0;
                            x0 = zeros(KSize(1),1);
                        else
                            if(norm((x0-b)-V*(V'*(x0-b))) > tol)
                                warning('%d dimension(s) of x0 in a preconditioner null-space have been removed',...
                                diff(KSize));
                            end
                            x0 = A*(x0-b);
                        end
                    else
                        error('Unexpected length of x0');
                    end
                end
            elseif(ind>1)
                s(ind:end) = [];
                V(:,ind:end) = [];
                ind = ind-1;
                %ASize = [ind,KSize(2)];
                temp = 1:ind;
                S = sparse(temp,temp,s);
                
                A = S*V';
                AInv = V*inv(S);    % This is now well defined
                
                if(numX0)
                    if(precondX0)
                        if(numX0==KSize(1))
                            % We neeed to eliminate a dimension provided by x0
                            % while rotating to the new space
                            x0 = (x0.'*U(:,1:ind)).';
                            warning('%d dimension(s) of x0 in a preconditioner null-space have been removed',...
                                KSize(1)-ind);
                            
                        elseif(numX0==ind)
                            % It is REALLY unlikely this was intended... Warn strongly
                            warning('x0 seemed to account for preconditioner null-space.  This is very likely a mistake.');
                            
                        else
                            error('Unexpected length of x0');
                        end
                        
                    else % X0 is not in the preconditioned space
                        if(numX0==KSize(2))
                            % We neeed to eliminate a dimension provided by x0
                            if(isempty(b))
                                b = x0;
                            else
                                b = b-x0;
                            end
                            x0 = zeros(ind,1);
                            warning('%d dimension(s) of x0 in a preconditioner null-space have been removed',...
                                KSize(1)-ind);
                        else
                            error('Unexpected length of x0');
                        end
                    end
                    
                end
                
            else
                error('Preconditioner has no null-space complement');
            end
            
        else                      
            % AInv will not be populated.  Trust A is correct and solve systems
            if(numX0)
                if(precondX0 && numX0~=KSize(1))
                    error('Unexpected length of x0');
                elseif(numX0==KSize(2))
                    if(isempty(b))
                        x0 = A*x0;
                    else
                        x0 = A*(x0-b);
                    end
                elseif(~precondX0)
                    error('Unexpected length of x0');
                end
            end
            
        end
    else    % We already have AInv
        % We are not going to form A if it doesn't already exist
        if(numX0)
            if(precondX0)   % NEW IF BLOCK %K is AInv down here
                if(numX0==KSize(2))
                    % Everything looks good
                else
                    error('Unexpected length of x0');
                end
                
            else % X0 is not in the preconditioned space
                if(numX0==KSize(1))
                    if(isempty(b))
                        % Use b to preserve x0
                        b = x0;
                        x0 = zeros(KSize(2),1);
                    else
                        % We neeed to project into the preconditioned space
                        x0 = AInv\(x0-b);
                        
                        % WARN???
                    end
                else
                    error('Unexpected length of x0');
                end
            end
        end
        
    end
    
end

end