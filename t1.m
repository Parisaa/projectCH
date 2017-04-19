disp('Generating channel trace...');
chan = ChanMod();
%divOrder = 4;
[nantxTx,nantyTx,nantxRx,nantyRx] = ...
    deal(chan.nantvec(1),chan.nantvec(2),chan.nantvec(3),chan.nantvec(4));
nantRx = nantxRx*nantyRx;
nw = 40;
nt = 1000;
nfreq =1;
unifMeas = 1;
% generate random uniformly spaced beamforming vectors
% for measuring the power in different direcions.
WRx = genRandW(nantxRx,nantyRx,nw,unifMeas); 

% Generate path parameters
chan.genPathParams();
np = size(chan.Srx,1);

chan.compBFDir();


% Generate SNR trace
chan.genSNRTrace();

% Get power values
lamRx = zeros(nw,nt,nfreq);
lamRxNT = zeros(nw,nfreq);
pRx = zeros(nw,1);

Srxt = zeros(np,nantRx,nt*nfreq);
SrxtC = zeros(np,nantRx,nt*nfreq);
Qrxt = zeros(nantRx,nantRx,nt*nfreq);
for iw = 1 : nw
    wi = WRx(:,iw);
    for ifreq = 1 : nfreq
        for it = 1 : nt
            p = chan.P(it,:,ifreq);
            ite =(ifreq-1)*nt+it;
            for inantRx = 1 : nantRx
                SrxtTemp = chan.Srx(:,inantRx).*transpose(sqrt(p));
                Srxt(:,inantRx,ite) = SrxtTemp;
                SrxtTempC = chan.Srx(:,inantRx)'.*sqrt(p);
                SrxtC(:,inantRx,ite) = transpose(SrxtTempC);
            end
            Qrxt(:,:,ite) = transpose(SrxtC(:,:,ite))*SrxtC(:,:,ite);
%             Qrxt(:,:,ite) = Srxt(:,:,ite)'*Srxt(:,:,ite);
            lamRx(iw,it,ifreq) = real(wi'*Qrxt(:,:,ite)*wi)/(wi'*wi);
            lamRxNT(iw,ifreq) = real(wi'*(chan.Srx'*chan.Srx)*wi)/(wi'*wi);
        end
          
    end   
end
meanAvg = mean(lamRx,2);
figure;
for i = 1 : nw
    hold on
    plot([1:nt],lamRx(i,:,1)-meanAvg(i));
end

figure
plot([1:nt],real(chan.P(:,1,1)));



return;

% nt = size(chan.snrTrueNB,1);
nt = 100;
t = (0:nt-1)*chan.tsMs*1e-3;
figure;
plot(t, [chan.snrTrue chan.snrRaw]);
legStr1 = {'snrTrue','snrRaw'};
% plot(t, [chan.snrTrueNB chan.snrTrue chan.snrRaw]);
% legStr1 = {'snrTrueNB','snrTrue','snrRaw'};
legend(legStr1);

% Get true channel and scale to the SNR
h0 = chan.Href;
[nt0,nsig] = size(h0);
hpow = mean(abs(h0(:)).^2);
h0 = sqrt(10^(0.1*snrSigAvg)/hpow)*h0;

% Get true wideband spec efficiency
snr = abs(h0).^2*snrDataScaleLin;
y0 = mean(log2(1 + snr), 2);

%y0 = 10*log10(sum(abs(h0).^2,2));
t = (0:nt-1)*chan.tsMs*1e-3;
figure;
plot(t, y0);

% Downsample data
nper = round(TperMs/chan.tsMs);
h0p = h0(1:nper:end,:);
y0p = y0(1:nper:end,:);
ntp = size(h0p,1);

% Add noise and make noisy measurements
w = sqrt(wvar/2)*(randn(ntp,nsig) + 1i*randn(ntp,nsig));
z = (abs(h0p+w).^2)/wvar;
snrp = max(z-1,0)*snrDataScaleLin;
yhat0p = mean(log2(1 + snrp), 2);

% Plot the raw estimated spectral efficiency
tp = t(1:nper:end);
figure;
plot(tp, [y0p yhat0p]);



return
% Compute raw estimate
[seHat, seVar] = seEst(z, snrDataScaleLin);
yhat1p = mean(seHat,2);
yvar1p = 1/nsig*mean(seVar,2);

methTest = {'jlest', 'ar'};
nmeth = length(methTest);
yhatMeth = zeros(ntp,nmeth);
for imeth = 1:nmeth
    
    meth = methTest{imeth};
    if strcmp(meth, 'jlest')
        % Jump linear estimator
        est = JLChanEst();
        est.set('dSEStd',0.2, 'pseg', 0.03, 'slopeStd', 0.5, 'len', 5, 'measStd', 0.0 );
        est.genSSModel();
        [xhatm, yhatm, plike] = est.filter(yhat1p,yvar1p);
    elseif strcmp(meth, 'ar')
        a = 0.5;
        yhatm = filter(a, [1 a-1], yhat0p, yhat0p(1));
    end
    
    % Store result
    yhatMeth(:,imeth) = yhatm;                

end

% Plot results
legStr = {'true', methTest{:}};
figure;
plot(tp, [y0p yhatMeth]);
grid on;
set(gca,'FontSize',16);
xlabel('Time (sec)');
ylabel('Spec eff (bps/Hz)');
legend(legStr);

