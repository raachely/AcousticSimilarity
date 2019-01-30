function result = usev0dialog(varargin)

%   Copyright 2007-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $
if (isempty(varargin))
    result = false;
else
    if (ischar(varargin{1}) && (strcmpi(varargin{1}, 'v0')))
        result = true;
    else
        result = false;
    end
end
end