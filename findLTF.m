function [LTFpeakPos] = findLTF(rxData)

global USESIM RSSthresh;

sigin = rxData;
ns = 16;

smsEng_windowSize = 16;
smsEng = 0; %smoothed energy level
smsEngQ = zeros(64,1);

SNRQsize = 64;
SNRQ = zeros(SNRQsize, 1);
SNRthresh = 4;
SNRQthresh = 0.6;

selfcorrQsize = 16;
selfcorrQ = zeros(selfcorrQsize, 1);
selfcorrQthresh = 0.8; 

prevPktPos = 0;
corr2Engthresh = 0.9;
avgCrosscorr = zeros(16, 1);

LTFpeakPos = 1;

while ns < length(sigin) - 161
	selfcorr = sum(sigin(ns+1:ns+16).*conj(sigin(ns-15:ns)));
    energy = sum(abs(sigin(ns+1:ns+16)).^2);

    % --- Update smoothed energy level and SNR queue -----
    smsEng = smsEng*(1-1/smsEng_windowSize) + ...
        1/smsEng_windowSize*abs(sigin(ns+1))^2;
    
    if (USESIM && abs(smsEng) <= 0)
        smsEng = 1e-8;
    end 
    
    %sliding window
    smsEngQ(end+1) = smsEng;
    smsEngQ(1) = [];
    
    if (smsEngQ(1) <= 0)
        SNR = 0;
    else
        SNR = 10*log10(smsEng/smsEngQ(1));
    end
    
    if (SNR > SNRthresh && energy > RSSthresh)
        SNRQ(end+1) = 1;
    else
        SNRQ(end+1) = 0;
    end
    
    SNRQ(1) = [];
    %{
    if DEBUG_OUT
        smsEngOut(ns) = smsEng;
        selfcorrOut(ns) = selfcorr;
        energyOut(ns) = energy;
    end
    %}
    maxSTFcrossCorrPeak = 0;
    
    % ----- Update self-correlation queue -------
    if ((ns - prevPktPos) > 2*SNRQsize && ...
        (abs(selfcorr)/energy>corr2Engthresh && ...
        abs(selfcorr)/energy<1/corr2Engthresh))
        selfcorrQ(end+1) = 1;
    else
        selfcorrQ(end+1) = 0;
    end
    selfcorrQ(1) = [];

    % ---- Decision: self-correlation and energy detection ----
    if (sum(SNRQ(end-16:end))/16< SNRQthresh || ...
        sum(selfcorrQ)/length(selfcorrQ) < selfcorrQthresh)
        %keep sliding
        ns = ns + 1;
        continue;
    end

    % ---- Self-corr peak detected ---- 
	if (DEBUG_OUT)
    	fprintf(1, '!!! Self-corr peak detected at sample %d level %f\n'...
            , ns, energy);
	end
    firstSTFpos = ns; %STF located
	selfcorrPeak = energy;

    % ---- Search for the first peak via cross-corr with STF ----
    af = 1; 
    peakPos = 0; 
    
    for k = 1:16
        crosscorrSTF = sum(sigin(ns+k:ns+k+15).*conj(STF_t(2:17)));
        
        if k==1
            avgCrosscorr(1) = abs(crosscorrSTF);
            continue;
        end
        avgCrosscorr(k) = af*abs(crosscorrSTF) + (1-af)*avgCrosscorr(k-1);
    end
    
    for k=2:15
        % detect a peak in the avgCrosscorr curve
        if (avgCrosscorr(k) > maxSTFcrossCorrPeak &&...
               avgCrosscorr(k) > avgCrosscorr(k-1) &&...
               (k<4||avgCrosscorr(k) > avgCrosscorr(k-3)) &&...
               (k<6||avgCrosscorr(k) > avgCrosscorr(k-5)) &&...
               avgCrosscorr(k) > avgCrosscorr(k+1) &&...
               (k>16-3 || avgCrosscorr(k) > avgCrosscorr(k+3)) &&...
               (k>16-5 || avgCrosscorr(k) > avgCrosscorr(k+5)))
                
            maxSTFcrossCorrPeak = avgCrosscorr(k);
            peakPos = k;
            %{
            if DEBUG_OUT
                %fprintf(1, 'peakPos %d, maxSTFcrossCorrPeak %f\n', ...
                %        peakPos, maxSTFcrossCorrPeak);
            end
            %}
        end
    end
    
    ns = ns + peakPos;


    % ------- Detect other peaks ---------
	endSTF = 0;
    peakPos = 0;
    % Calculate avgCrosscorr of all the following sample positions
    numPeak = 1;
    psamples = 16*9;
    
    for k = 1:psamples %max number of peaks is 9
        crosscorrSTF = sum(sigin(ns+k:ns+k+15).*conj(STF_t(2:17)));       
        if k==1
            avgCrosscorr(1) = abs(crosscorrSTF);
            continue;
        end
        avgCrosscorr(k) = af*abs(crosscorrSTF) + (1-af)*avgCrosscorr(k-1);
    end 
    
    for k=2:psamples 
        if (avgCrosscorr(k) > maxSTFcrossCorrPeak*0.8 &&...
           avgCrosscorr(k) > avgCrosscorr(k-1) && ... 
           (k<3||avgCrosscorr(k) > avgCrosscorr(k-2)) && ... 
           (k<4||avgCrosscorr(k) > avgCrosscorr(k-3)) &&...
           (k<6||avgCrosscorr(k) > avgCrosscorr(k-5)) &&...
           (k>psamples-1 || avgCrosscorr(k) > avgCrosscorr(k+1)) &&...
           (k>psamples-2 || avgCrosscorr(k) > avgCrosscorr(k+2)) &&...
           (k>psamples-3 || avgCrosscorr(k) > avgCrosscorr(k+3)) &&...
           (k>psamples-5 || avgCrosscorr(k) > avgCrosscorr(k+5)))
   
            numPeak = numPeak + 1;
            peakPos = k;
            if DEBUG_OUT
            %fprintf(1, 'peak at %d, avgCrosscorr=%g\n', ...
            %        ns+k, avgCrosscorr(k));
            end
        end

        if k-peakPos > 18 % no more peaks
            if DEBUG_OUT
                %fprintf(1, 'No more peaks at %d\n', ns+k);
            end
            break;
        end
         
		
		% End of STF: self-corr ends (if there is a peak nearby,
        % the distance will not be greater than 16)
        tempselfcorr = sum(sigin(ns+k+1:ns+k+16).*conj(sigin(ns+k-15:ns+k)));
        tempselfcorr = abs(tempselfcorr);
        tempEng = sum(abs(sigin(ns+k+1:ns+k+16)).^2);
        if (tempselfcorr/tempEng>corr2Engthresh && ...
            tempselfcorr/tempEng<1/corr2Engthresh )
            selfcorrQ(end+1) = 1;
        else
            selfcorrQ(end+1) = 0;
        end
        selfcorrQ(1) = [];

        if (sum(selfcorrQ(end-15:end))/16 < selfcorrQthresh)
            if DEBUG_OUT
               fprintf(1, 'case 1 End of STF at %d\n', ns+k);
            end
			endSTF = ns + k + 16;
			break; 
        end		
        
        if (ns+k-firstSTFpos > 160)
            if DEBUG_OUT
                fprintf(1, 'case 2 End of STF at %d\n', ns+k);
            end
            endSTF = ns+k;
            break;
        end
        
    end %end for k=2:psamples

    if (numPeak < 2) 
            SNRQ = zeros(SNRQsize,1);
            selfcorrQ = zeros(selfcorrQsize,1); 
            fprintf(1, 'pkt lost, num of detected STF peaks is %d\n',...
                numPeak);
            continue; 
    end

    %STF ends, start point for LTF detection
    ns = endSTF;   
    
    % Get the exact position of LTF via crosscorr    
    for t = ns+10:ns+30        
        crosscorrLTF(t-ns-9) = abs(sum(sigin(t+1:t+16).*conj(LTF_t(1:16))));
    end   
    
	[~, peakPos] = max(crosscorrLTF(1:end));

    LTFpeakPos = ns+10-1 + peakPos;
    
end

end