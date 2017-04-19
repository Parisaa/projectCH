function W = genRandW(nantx,nanty,nw, unif)
% genRandW:  Generates random weights aligned in random angles

% Parameters
dsep = 0.5;     % antenna spacing in wavelengths
if (nargin < 4)
    unif = false;
end

% Get displacements per antenna
dx = repmat( (0:nantx-1)',1,nanty);
dx = dsep*dx(:);
dy = repmat( (0:nanty-1),nantx,1);
dy = dsep*dy(:);

% Generate random angles
if (unif)
    nwx = round(sqrt(nw));
    nwy = ceil(nw/nwx);
    angx = repmat(linspace(0,pi*((nwx-1)+0.5)/nwx,nwx)',1,nwy);
    angx = angx(:);
    angy = repmat(linspace(0,pi*((nwy-1)+0.5)/nwy,nwy),nwx,1)-pi/2;
    angy = angy(:);
else
    angx = 2*pi*rand(nw,1);
    angy = 2*pi*rand(nw,1);
end

% Compute spatial signature of each angle
theta = dx*(cos(angx).*sin(angy))' + dy*(sin(angx).*sin(angy))';
W = exp(2*pi*1i*theta);

% Puncture to reduce back down to nw
I = randperm(nw);
W = W(:,I);
