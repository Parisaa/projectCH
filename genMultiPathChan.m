function [Stx,Srx,angtx,angrx] = genMultiPathChan(nantvec,param)
% genMultiPathChan:  Generates a random multipath MIMO channel
%
% Returns matrices Srx = npath x nrx and Stx npath x ntx
% of the RX and TX spatial response (including gain) for each path.
% Given TX and RX BF vectors wrx and wtx, the gain is given by
%
%   h  = sum( (Srx*wrx).*(Srx*wtx) )
%
% and the average power is
%
%   hpow = sum( abs( (Srx*wrx).*(Stx*wtx) ).^2 )

% Get parameters.  
adly = param.adly;
angStdMean = param.angStdMean;
ashad = param.ashad;
nclam = param.nclam;

% Parameters
dsep = 0.5;     % antenna spacing in wavelengths
nsub = 20;      % num subpaths per cluster
ndir = 2;       % num of directions (TX and RX)

% Generate random number of clusters
nct = max(1, poissrnd(nclam, 1,1));

% Compute random powers on the clusters
% and the real gains on each subpath
powc = 10.^(0.1*ashad*randn(nct,1) + adly*log10(rand(nct,1)));
powc = powc / sum(powc);
gsub = repmat(sqrt(powc/nsub)',nsub,1);
gsub = gsub(:);

for idir = 1:ndir
 
    % Compute offset for antenna structures.
    if (idir == 1)
        ioff = 0;
    else
        ioff = 2;
    end

    % Get antenna dimnesions
    nantx = nantvec(ioff + 1);
    nanty = nantvec(ioff + 2);
    nant = nantx*nanty;
    
    % Compute displacement per antenna
    dx = repmat( (0:nantx-1)',1,nanty);
    dx = dsep*dx(:);
    dy = repmat( (0:nanty-1),nantx,1);
    dy = dsep*dy(:);           
    
    % Generate random angles for the cluster
    % The elevation angle is the same for all cluster
    angcx = 2*pi*rand(nct,1);
    angcy = 2*pi*rand(1)*ones(nct,1);
    angStdx = -log(rand(nct,1))*2*pi/360*angStdMean(ioff + 1);
    angStdy = -log(rand(nct,1))*2*pi/360*angStdMean(ioff + 2);
    
    % Generate random angles for each subpath
    angx = repmat(angcx',nsub,1) + repmat(angStdx',nsub,1).*randn(nsub,nct);    
    angy = repmat(angcy',nsub,1) + repmat(angStdy',nsub,1).*randn(nsub,nct); 
    angx = angx(:);
    angy = angy(:);
    
    % Compute spatial signature of each path
    theta = dx*(cos(angx).*sin(angy))' + dy*(sin(angx).*sin(angy))';
    S = exp(2*pi*1i*theta).*repmat(sqrt(gsub)',nant,1);
    
    % Store the spatial matrix
    if (idir==1)
        Stx = S';
        angtx = [angx angy];
    else
        Srx = S';
        angrx = [angx angy];
    end
    
end