clear all;
disp('Generating channel trace...');
chan = ChanMod();
%divOrder = 4;
[nantxTx,nantyTx,nantxRx,nantyRx] = ...
    deal(chan.nantvec(1),chan.nantvec(2),chan.nantvec(3),chan.nantvec(4));
nantRx = nantxRx*nantyRx;
nw = 40;

blk = 1; 

unifMeas = 1;
% generate random uniformly spaced beamforming vectors
% for measuring the power in different direcions.
WRx = genRandW(nantxRx,nantyRx,nw,unifMeas);

% Generate path parameters
chan.genPathParams();

% Compute UE velocity angle.  To look at the worst case,
% the UE is directed in the centroid of the angles of arrival
np = size(chan.Srx,1);
grx = sum(abs(chan.Srx).^2,2);
nang = 2;
avgAng = sum(repmat(grx,1,nang).*exp(1i*chan.angrx))./repmat(sum(grx),1,nang);
cosDir = real(repmat(conj(avgAng),np,1).*exp(1i*chan.angrx));
cosDir = prod(cosDir,2);
wi = genRandW(4,4,1,1); 
% Compute the Doppler shift of each path
ueVelKmh = 3;
ueVelMs = ueVelKmh/3.6; % UE velocity in m/s
vc = 3e8;               % speed of light in m/s
chan.tsMs = 1e-3*4*32;
deltheta = cosDir*chan.tsMs*chan.fcGHz*1e6*ueVelMs/vc;

% Compute the beamformed channel:  H(t) = urx'*Srx'
nt = size(chan.blkdB,1);
H = zeros(nt,chan.nfreq);

lamRx = zeros(nw,1);
pRx = zeros(nw,1);
nfreq = 4;
np = size(chan.Srx,1);
nt = 1000;
Srxt = zeros(np,nantRx,nt*nfreq);
Qrxt = zeros(nantRx,nantRx,nt*nfreq);

for ifreq = 1:chan.nfreq
    % Generate phase rotations for each path
    theta0 = rand(1,np)*2*pi;
    P = exp(1i*(repmat(theta0,nt,1) + (0:nt-1)'*deltheta'));

    % Compute the power in the path
%     H(:,ifreq) = sum(repmat(chan.srx.*chan.stx,nt,1).*P,2);


    for it = 1 : nt
        p = P(it,:);
        ite =(ifreq-1)*nt+it;
        for inantRx = 1 : 16
            SrxtTemp = chan.Srx(:,inantRx).*transpose(p);
            Srxt(:,inantRx,ite) = SrxtTemp;
        end
        Qrxt(:,:,ite) = Srxt(:,:,ite)'*Srxt(:,:,ite);
        lamRx(it,ifreq) = (norm(wi'*Qrxt(:,:,ite)*wi).^2)/(wi'*wi);
    end
end

% Normalize to an average gain of 1
chan.wvar = 10.^(-0.1*chan.snrSig);
HpowAvg = mean(abs(H(:)).^2);
H = H / sqrt(chan.nfreq*HpowAvg);

% Compute true wideband SNR
Hpow = sum(abs(H).^2,2);
chan.snrTrueNB = 10*log10(Hpow/chan.wvar);
chan.snrTrue = chan.snrTrueNB + chan.blkdB;

% Compute estimated SNR
if (chan.blk)
    H = H.*repmat(10.^(0.05*chan.blkdB),1,chan.nfreq);
end
Hn = H + sqrt(chan.wvar/2)*(randn(nt,chan.nfreq)+1i*randn(nt,chan.nfreq));
Hpown = sum(abs(Hn).^2,2);
chan.snrRaw = 10*log10( max(chan.snrMin, Hpown/chan.wvar-chan.nfreq));

% Save reference symbols
chan.Href = H;
chan.Hrefn = Hn;