% Multiple subscripts from linear index
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   varargout = ind2sub(siz,ind,<vec>)
%
% PURPOSE:
%   This function is a rewrite of Matlab's builtin ind2sub.  It performs the
%   same functions as the built-in, faster, and with a more sane interface.
%   Matlab allows the last index to be linear regardless of the number of
%   remaining dimensions.  Thus, if you ask for fewer outputs than the
%   dimensionality indicated by size, the last dimension is linearized.  If,
%   however, you ask for only a single output, this code also allows a vector
%   to be returned (this is a common use case).
%  
% INPUT:
%   siz         - [1 x n] Size of the array being linearly indexed
%   ind         - Vector of k linear indicies
%   vec         - [Logical indicating whether a vector output should be returned
%                 Default: false (Matlab's behavior)
%
% OUTPUT:
%   varargout   - The requested sub-indicies.  If ind is an array, these will
%                 also be arrays.  If the vec flag is passed this will be [k n]
%
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
% SEE ALSO
%   elmat/ind2sub
%-------------------------------------------------------------------------------
%
%   Copyright (C) 2009-2019 Joel W. LeBlanc
%   This program comes with ABSOLUTELY NO WARRANTY;
%   This is free software, and you are welcome to redistribute it under certain
%   conditions; see gpl-3.0.txt for more details.
%
%-------------------------------------------------------------------------------
%}

% Copyright (C) 2009-2017  Joel W. LeBlanc <jwleblan@gmail.com>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License version 3, as published
% by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
function varargout = ind2subVec(siz,ind,varargin)

n = numel(siz);

% We will produce multiple output vectors
if(nargout>1 || ~(nargin>2 && varargin{1}))
    if(nargout<=n)
        n = nargout;
        k = [1 cumprod(siz(1:n-1))];
    else
        k = [1 cumprod(siz(1:n-1))];
        k(n) = k(n)*prod(siz(n+1:end));
        varargout(n+1:nargout) = {ones(size(ind))};
    end
    
    for i=n:-1:2
        t = ceil(ind/k(i));
        ind = ind-(t-1)*k(i);
        varargout{i} = t;
    end
    varargout{1} = ind;
     
else                            % We will produce 1 vector output
    % Escape empty outputs... needed because matlab will not properly broadcast
    % dimensions from empty vector.
    % i.e. t = zeros([0 4]); t(:,3) = 2;
    if isempty(ind)
        varargout{1} = zeros([0 n]);
        return
    end
    k = [1 cumprod(siz(1:n-1))];

    for i=n:-1:2 
        t = ceil(ind/k(i));
        ind = ind-(t-1)*k(i);
        varargout{1}(:,i) = t;
    end
    varargout{1}(:,1) = ind;
end
     
end
