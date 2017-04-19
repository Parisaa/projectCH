classdef ChanMod < hgsetget
    % Millimeter channel model with blocking
    
    properties
        
        % Path parameters
        Stx, Srx;       % Spatial signatures of each path
        angtx, angrx;   % Angles of each path
        stx, srx;       % Coefficients after BF vectors are applied
        P;
        
        % SNR parameters
        snrDat = -10;    % nominal data SNR in best direction
        WtotMHz = 500;   % total bandwidth
        TperMs = 10;     % sync signal period in ms between
                         % directional measurements
        TsigMs;          % sync signal duration in ms
        snrSig;          % nominal sync signal SNR in dB   
        
%         % Blocking data file
%         blk = true;  % 0=no blocking
%         fileId = 3;  % Blocking file id 
%         fnames = {'22032','634625', '430088'};
        
        % Blocking data
        blkdB;  % Blocking value in dB
        tsMs;   % sample time in ms.
        
        % SNR measurement parameters
        nfreq = 1;       % number of frequency measurements.
        fcGHz = 28;      % frequency for selecting covariance parameters
        ueVelKmh = 3;   % UE velocity in km/h
        nantvec = [8 8 4 4];        % Antenna vector
        sigOver = 0.05;  % Fraction overhead for synchronization signals
        
        
        % SNR values
        snrMin = 1e-2;  % minimum SNR in linear scale
        snrTrueNB;      % True wideband SNR without blocking
        snrTrue;        % True wideband SNR with blocking
        snrRaw;         % Noisy SNR measurement
        Href, Hrefn;    % Reference symbols raw and noisy values
        wvar;           % noise power
        
    end
    
    methods
        
        % Constructor
        function obj = ChanMod()
        end
        
        
        
        % Read blocking data
        function readBlkData(obj)
            % Read the file
            fn = sprintf('data/%s.tmp',obj.fnames{obj.fileId});
            fileID = fopen(fn,'r');
            x = fread(fileID);
            fclose(fileID);
            ns0=128; % Throw away the first 128 numbers, which are useless anyway.
            
            % Reshape the data to 8 words and convert to U64
            nx = size(x,1);
            nb = 8;
            ns = floor(nx/nb)-ns0;
            x = reshape(x(ns0*nb+1:(ns0+ns)*nb),nb,ns);
            I = [0 1 2 3 4 5 6 7]; % This is called "Little Endian" byte order.
            y = 2.^(8*I)*x;
            
            % Reshape the data to groups of NFFT samples.
            ny = length(y);
            nfft = 128;
            nsym = floor(ny/nfft);
            y = reshape(y(1:nsym*nfft),nfft,nsym);
            
            % Find the noise level in each symbol by taking the lowest
            % pnoise percent of each symbol
            pnoise = 0.5;
            ysort = sort(y);
            nnoise = round(pnoise*nfft);
            yn = mean(ysort(1:nnoise,:));
            
            % Find power by capturing signal level above noise
            thresh = 1;
            yt = repmat(thresh*yn,nfft,1);
            y1 = y.*(y > yt);
            ypow = 10*log10(sum(y1));
            
            % Measure power relative to peaks and store in blocking data
            p = 0.3;
            ys = sort(ypow,'descend');
            yp = mean(ys(1:round(p*nsym)));
            obj.blkdB = ypow'-yp;
            obj.tsMs = 1e-3*4*32;   % Sample time in ms
        end
        
        % Generate path parameters
        function genPathParams(obj)
            % Load the channel parameters
            load data\chanParam28;
            
            % Get parameters.
            param.adly = adly;
            param.angStdMean = angStdMean;
            param.ashad = ashad;
            param.nclam = nclam;
            [obj.Stx,obj.Srx,obj.angtx,obj.angrx] = ...
                genMultiPathChan(obj.nantvec,param);
        end
        
        % Compute the BF directions and SNR
        function compBFDir(obj)
            
            % Compute the optimal beamforming vectors in the RX direction
            [Urx,~] = svd(obj.Srx',0);
            obj.srx = Urx(:,1)'*obj.Srx';
            
            % The TX sends omni directionally.  So, we transmit out of the
            % one of the antennas.
            obj.stx = obj.Stx(:,1)';
            
            % Compute the BF gain in the TX direction.  This will not be
            % available on the sync channel
            s = svd(obj.Stx,0);
            ntx = size(obj.Stx,2);
            txBFGain = ntx*max(s)/sum(s);
            
            % Sync signal duration.
            nrx = size(obj.Srx,2);
            obj.TsigMs = obj.TperMs*obj.sigOver/nrx;
            
            % Compute the signal SNR
            %  snrSig0 = snrDat*Wtot*Tsig/txBFGain;
            obj.snrSig = obj.snrDat + 10*log10(obj.WtotMHz*obj.TsigMs*1e3/txBFGain);
            
        end         
        
        % Find the best TX and RX directions
        function genSNRTrace(obj)
            
            obj.tsMs = 1e-3*4*32;
            
            % Compute UE velocity angle.  To look at the worst case,
            % the UE is directed in the centroid of the angles of arrival
            np = size(obj.Srx,1);
            grx = sum(abs(obj.Srx).^2,2);
            nang = 2;
            avgAng = sum(repmat(grx,1,nang).*exp(1i*obj.angrx))./repmat(sum(grx),1,nang);
            cosDir = real(repmat(conj(avgAng),np,1).*exp(1i*obj.angrx));
            cosDir = prod(cosDir,2);
            
            % Compute the Doppler shift of each path
            ueVelMs = obj.ueVelKmh/3.6; % UE velocity in m/s
            vc = 3e8;               % speed of light in m/s
            deltheta = cosDir*obj.tsMs*obj.fcGHz*1e6*ueVelMs/vc;
                        
            % Compute the beamformed channel:  H(t) = urx'*Srx'
%             nt = size(obj.blkdB,1);
            nt = 1000;
            HP = zeros(nt,obj.nfreq);
            theta0 = zeros(obj.nfreq,np);
            

            for ifreq = 1:obj.nfreq
                % Generate phase rotations for each path
                theta0(ifreq,:) = rand(1,np)*2*pi;
                obj.P(:,:,ifreq) = exp(1i*(repmat(theta0(ifreq,:),nt,1) + (0:nt-1)'*deltheta'));
                
                % Compute the power in the path
                HP(:,ifreq) = sum(repmat(obj.srx.*obj.stx,nt,1).*obj.P(:,:,ifreq),2);
                
            end
            
            
            
            % Normalize to an average gain of 1
            obj.wvar = 10.^(-0.1*obj.snrSig);
            HpowAvg = mean(abs(HP(:)).^2);
            HP = HP / sqrt(obj.nfreq*HpowAvg);
            
            % Compute true wideband SNR
            Hpow = sum(abs(HP).^2,2);
            obj.snrTrue = 10*log10(Hpow/obj.wvar);
%             obj.snrTrue = 10*log10(Hpow/obj.wvar);
%             obj.snrTrue = obj.snrTrueNB + obj.blkdB;
            
            % Compute estimated SNR
%             if (obj.blk)
%                 H = H.*repmat(10.^(0.05*obj.blkdB),1,obj.nfreq);
%             end
            Hn = HP + sqrt(obj.wvar/2)*(randn(nt,obj.nfreq)+1i*randn(nt,obj.nfreq));
            Hpown = sum(abs(Hn).^2,2);
            obj.snrRaw = 10*log10( max(obj.snrMin, Hpown/obj.wvar-obj.nfreq));
            
            % Save reference symbols
            obj.Href = HP;
            obj.Hrefn = Hn;
        end
        
    end
    
end