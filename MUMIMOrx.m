% OFDM MU-MIMO packet decoder
function [finalSINR,avgThru] = MUMIMOrx(rxID, rxCompare)

global VLTF_shift LTF_shift LTF_t STF_t numTxAntenna numRxAntenna numOFDMsymb newCHout rxData_board DECODE_DATA USESIM osamp DEBUG_OUT signalPower noisePower sinr_subc useShannon usePilotCompen pktType modulationM useNewSINR modeID; 
avgThru = 0;
finalSINR = 0;
otherTo48=zeros(1,64);
RSSthresh = 0.0005; %sum(abs(sigin(1, 90:106)).^2);
rxaID = 1;  %ID of Rx antenna
skipsamples = 50;


if USESIM 
    osamp = 1;
    maxRead = 10000;
    for rxa = 1:numRxAntenna
        if (modeID == 1)||(modeID == 2)
            fname = sprintf('tmp_sim_modulated_Rx%d.txt', rxID);
        elseif modeID == 3 % open loop mimo
            fname = sprintf('tmp_sim_modulated_Rx%d.txt', rxa);
        end
        fin = fopen(fname, 'rb');
        mt = fread(fin, [2, maxRead], 'float32');
        fclose(fin); 
        sigin(rxa,:) = mt(1,:) + 1i*mt(2,:);
    end
else 
    % Process the received samples
    for rxa = 1:numRxAntenna
        rxData_this = rxData_board(:,rxa,rxID);   
        sigin(rxa,:) = decimate(rxData_this, osamp);         
        % Remove DC Offset (DCO) from RxData
        %[sigin(rxa,:)] = warplab_correctDCO(sigin(rxa,:),16);
    end



end %end if USESIM



% =========== Initialize parameters ===========
numPkt = 0;
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

% =========== Loop through all samples ============
ns = skipsamples;

if (USESIM)
    for k = 1:length(smsEngQ)
        smsEngQ(k) = 1e-8;
    end
end

while ns < length(sigin(rxaID,:)) - 161
    
%% start: locate LTF
   % --- Update self-corr and energy level
    selfcorr = sum(sigin(rxaID, ns+1:ns+16).*conj(sigin(rxaID, ns-15:ns)));
    energy = sum(abs(sigin(rxaID, ns+1:ns+16)).^2);

    % --- Update smoothed energy level and SNR queue -----
    smsEng = smsEng*(1-1/smsEng_windowSize) + 1/smsEng_windowSize*abs(sigin(rxaID, ns+1))^2;
    
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
    	fprintf(1, '!!! Self-corr peak detected at sample %d level %f\n', ns, energy);
	end
    firstSTFpos = ns; %STF located
	selfcorrPeak = energy;

    % ---- Search for the first peak via cross-corr with STF ----
    af = 1; 
    peakPos = 0; 
    
    for k = 1:16
        crosscorrSTF = sum(sigin(rxaID, ns+k:ns+k+15).*conj(STF_t(2:17)));
        %{
        if DEBUG_OUT 
            crosscorrOut(ns+k) = crosscorrSTF; 
            selfcorrOut(ns+k) = sum(sigin(rxaID, ns+1:ns+16).*conj(sigin(rxaID, ns-15:ns)));
            energyOut(ns+k) = sum(abs(sigin(rxaID, ns+1:ns+16)).^2);
        end
        %}
        
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
        crosscorrSTF = sum(sigin(rxaID, ns+k:ns+k+15).*conj(STF_t(2:17)));       
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
         
		
		% End of STF: self-corr ends (if there is a peak nearby, the distance will not be greater than 16)
        tempselfcorr = sum(sigin(rxaID, ns+k+1:ns+k+16).*conj(sigin(rxaID, ns+k-15:ns+k)));
        tempselfcorr = abs(tempselfcorr);
        tempEng = sum(abs(sigin(rxaID, ns+k+1:ns+k+16)).^2);
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
            fprintf(1, 'pkt lost, num of detected STF peaks is %d\n', numPeak);
            continue; 
    end

    %STF ends, start point for LTF detection
    ns = endSTF;   
    
    % Get the exact position of LTF via crosscorr    
    for t = ns+10:ns+30        
        crosscorrLTF(t-ns-9) = abs(sum(sigin(rxaID, t+1:t+16).*conj(LTF_t(1:16))));
    end   
    
	[~, peakPos] = max(crosscorrLTF(1:end));

    LTFpeakPos = ns+10-1 + peakPos;
    
    

%% end: locate LTF


    %one packet is detected
    numPkt = numPkt + 1;   
    
    
    st = LTFpeakPos - 32;
	if DEBUG_OUT
    	%fprintf(1, 'LTFpeakPos = %d ns = %d\n', LTFpeakPos, ns);
    end 
    
    %Estimate freq offset    
    fftoffset = 0;
    op = 32-fftoffset;
    co = sum(sigin(rxaID, (st+op+1):(st+op+64)).*conj(sigin(rxaID, (st+op+65):(st+op+65+63))));
    iFreq = angle(co)/64;

    % --- Estimate the channel from each tx antenna to each rx antenna --- 
    for atnt = 1:numTxAntenna
        
        for atn = 1:numRxAntenna
            % --- freq compensation ---
            delTheta = iFreq;
            for k = 1:64
                e0(k) = sigin(atn, st + op + (atnt-1)*160 + k);
                e0(k) = e0(k) * exp(1i*delTheta*((atnt-1)*160 + k-1));
            end


                % --- estimate channel of each subcarrier ---
                LTF_FFT = fft(e0, 64);
                for k = 2:64 %zero freq 0+1 does not carry data
                    if (k > 26+1 && k < 64-26+1)
                        continue;                
                    end
                    % --- Rx antenna atn, Tx antenna atnt, subcarrier k ---
                    LTF_CH(atn, atnt, k) = LTF_FFT(k)/LTF_shift(k);
                end
                
                
                
            % estimate signal power and noise floor
            % Use RX antenna 1    
            %if (atn == 1)

                for k = 1:64*2            
                    e1(k) = sigin(atn, st + op -7 + (atnt-1)*160 + k);
                    e1(k) = e1(k) * exp(1i*delTheta*((atnt-1)*160 + k-1)); 
                end
                freqLTF1 = fft(e1(1:64));
                freqLTF2 = fft(e1(65:128));                

                noisePower(rxID,atnt) = sum(abs((freqLTF1-freqLTF2)*(freqLTF1-freqLTF2)'));
                signalPower(rxID,atnt) = sum(abs(freqLTF1).^2) - noisePower(rxID,atnt);
                noisePower(rxID,atnt) = noisePower(rxID,atnt)/64*52;
                
                LTF_SINR(atnt) = 10*log10(signalPower(rxID,atnt)/noisePower(rxID,atnt));
                fprintf(1, '****** RX %d antenna %d, LTF from Tx ant %d, signal power %g, noise power %g, Estimated SINR %g dB\n', rxID, atn, atnt, signalPower(rxID,atnt), noisePower(rxID,atnt), LTF_SINR(atnt));

            %end
            

            
        end %end for atn = 1:numRxAntenna
        
    end %end for atnt = 1:numTxAntenna 
    
    joint_SINR(rxID) = 10*log10(sum(signalPower(rxID,:))/sum(noisePower(rxID,:)));            
    %fprintf(1, '****** RX %d, Joint SINR %g dB\n', rxID, joint_SINR(rxID));
     
    % ================ Decode VHT-LTF ===================
        % --- VHT-LTF estimates the combined effects of precoding and channel gain
        st = st + numTxAntenna*160;
        if (DEBUG_OUT)
            %fprintf(1, 'VLTF starts %d\n', st);
            for (m = st:st+79)
                energyOut(m) = sum(abs(sigin(rxaID, m+1:m+16)).^2);
            end 
        end
      
        % --- freq compensation ---
        delTheta = iFreq; %(ft1(1)+ft2(1))/2; %iFreq;
        op = 16-fftoffset;
        clear e0;
        for k = 1:64    
            %get data from Rx antenna 1
            e0(k) = sigin(1, st + op + k);
            e0(k) = e0(k) * exp(1i*delTheta*(numTxAntenna*160 + k-1));
        end
        
        % --- estimate channel of each subcarrier ---
        VLIF_FFT = fft(e0, 64);
        for k = 2:64
            if (k > 27 && k < 64-25)
                continue;
            end
            % --- Rx antenna atn, subcarrier k
            VLTF_CH(k) = VLIF_FFT(k)/VLTF_shift(k,rxID);  
        end 
	

     
	if (DEBUG_OUT)
       
        if rxID ==1
            figure(211);
            subplot(2,2,1);
            plot(abs([squeeze(LTF_CH(1,1,33:64)); 0; squeeze(LTF_CH(1,1,2:32))]), ...
                '-+r', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('LTF channel gain');
            title(stt); 
            subplot(2,2,2);
            plot(unwrap(angle([squeeze(LTF_CH(1,1,33:64)); 0; squeeze(LTF_CH(1,1,2:32))])), ...
                '-or', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('LTF phase');
            title(stt);  
            subplot(2,2,3);
            plot(abs([VLTF_CH(33:64) 0 VLTF_CH(2:32)]), ...
                '-xb', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('VLTF channel gain');
            title(stt); 
            subplot(2,2,4);
            plot(unwrap(angle([VLTF_CH(33:64) 0 VLTF_CH(2:32)])), ...
                '-xb', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('VLTF phase');
            title(stt);
        else  
       
            figure(211);
            subplot(2,2,1);
            plot(abs([squeeze(LTF_CH(1,2,33:64)); 0; squeeze(LTF_CH(1,2,2:32))]), ...
                '-+r', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('LTF channel gain');
            title(stt); 
            subplot(2,2,2);
            plot(unwrap(angle([squeeze(LTF_CH(1,2,33:64)); 0; squeeze(LTF_CH(1,2,2:32))])), ...
                '-or', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('LTF phase');
            title(stt);  
            subplot(2,2,3);
            plot(abs([VLTF_CH(33:64) 0 VLTF_CH(2:32)]), ...
                '-xb', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('VLTF channel gain');
            title(stt); 
            subplot(2,2,4);
            plot(unwrap(angle([VLTF_CH(33:64) 0 VLTF_CH(2:32)])), ...
                '-xb', 'LineWidth', 2, 'MarkerSize', 10);
            stt = sprintf('VLTF phase');
            title(stt); 
        end
	end 
	 
	if DEBUG_OUT
        
        for t = 1:160*numTxAntenna
            if (ns+64 > length(sigin(rxaID,:))) break; end
            
            selfcorrOut(ns) = sum(sigin(rxaID, ns+1:ns+16).*conj(sigin(rxaID, ns-15:ns)));
            energyOut(ns) = sum(abs(sigin(rxaID, ns+1:ns+16)).^2);            
            crosscorrOut(ns) = sum(sigin(rxaID, ns+1:ns+64).*conj(LTF_t(1:64)));
            ns = ns + 1;
        end
	else 
		ns = ns + 160*numTxAntenna; 
	end 
    
    %% ================ Decode OFDM data symbols ===================
    if DECODE_DATA
	inb = load('databits.dat');
    %totalSymb = length(inb(:,1))/48; %total number of OFDM symbols 
    nsymbol = 0;
    
	if (DEBUG_OUT)
		%fprintf(1, 'First sample of data at %d\n', st+80+1);
    end
    
    symbcount = ones(1,numRxAntenna);% bit symbol
	numErr = zeros(1,numRxAntenna);
    %while nsymbol < totalSymb % decode 10 OFDM symbols 

while nsymbol < numOFDMsymb % decode 10 OFDM symbols     
    nsymbol = nsymbol + 1; 
    
	st = st + 80;
	ns = ns + 80;
	 
	for (m = ns:ns+79)
    	energyOut(m) = sum(abs(sigin(rxaID, m+1:m+16)).^2);
	end
	 
    op = 16-fftoffset;

	clear e0;
    for k = 1:64 % freq compen
        for rxa = 1:numRxAntenna
            e0(rxa,k) = sigin(rxa, st + op + k);
            e0(rxa,k) = e0(rxa,k) * exp(1i*delTheta*(numTxAntenna*160+80+(nsymbol-1)*80+k-1));
        end
    end

    %FFT to recover frequency domain data    
    for rxa = 1:numRxAntenna
        datafft(rxa,:) = fft(e0(rxa,1:64), 64);
    end

      
    % k=1 is 0 freq
    for k = 2:64        
        if (k > 27 && k < 64-26+1)
            continue;
        end        
        if numTxAntenna>1
            if modeID <= 2
                datasymb(k) = datafft(k)/VLTF_CH(k);
            elseif modeID == 3 %open loop mimo
                datasymb(:,k) = inv(LTF_CH(:,:,k))*datafft(:,k);
            end
        else
            datasymb(k) = datafft(k)/LTF_CH(k);
        end
    end  

    
    % ------------- pilot compensation
    if usePilotCompen == 1
     for rxa = 1:numRxAntenna
        th1 = angle(datasymb(rxa,64-21+1));
        th2 = angle(datasymb(rxa,64-7+1));
        th3 = angle(datasymb(rxa,7+1));

        th4 = angle(datasymb(rxa,21+1));

        % phase drift between two adjacent subcarriers         
        dTheta = ((th3-th1)+(th4-th2))/2/(21+7);
        
        %avgTheta is the phase corresponding to the middle subcarrier 0
        %-21 + (14+28+42)/4 = 0 
        avgTheta = th1 + (((th2-th1) + (th3-th1) + (th4-th1)) / 4);

        %phase of leftmost subcarrier
        th = avgTheta-26*dTheta; 

        for k = (64-26+1):64
            datasymb(rxa,k) =  datasymb(rxa,k) * exp(-1i*th);             
            th = th + dTheta; 
        end 
        
        for k = 1:27
            datasymb(rxa,k) =  datasymb(rxa,k) * exp(-1i*th);            
            th = th + dTheta; 
        end
     end
        
    end
    % ----- end pilot compensation -----
     

%% ====== decoding ======
for rxa = 1:numRxAntenna
        for k = 1:64
          
            %other way to map subc, which is used in precoding
            if k<=31
                otherCount(k)=k+32;
            else
                otherCount(k)=k-32;
            end
            
            
            %ignore pilot, unused subc, pilot 
            if (k==1 || (k > 27 && k < 64-25) || k == (7+1) || k == (21+1) ...
                || k == (64-21+1) || k == (64-7+1))
                continue;
            end

            datasymb1(rxa,symbcount(rxa)) = datasymb(rxa,k);
            
            if symbcount(rxa) <=48
                otherTo48(otherCount(k))=symbcount(rxa);
            end
            
            %demodulation
         
            
            % the ID of the dot on the constellation
            
            %anti-normalization just for the wierd built in qamdemod deisgn
            if modulationM == 2
                scalingFactor = 1; 
            elseif modulationM == 4
                scalingFactor = norm(sqrt(modulationM)-1+1i*floor(sqrt(modulationM)-1)); 
            elseif modulationM == 8            
                scalingFactor = norm(3+1i);            
            elseif modulationM == 16              
                scalingFactor = norm(sqrt(modulationM)-1+1i*floor(sqrt(modulationM)-1));
            elseif modulationM == 32              
                scalingFactor = norm(5+7i);
            elseif modulationM == 64              
                scalingFactor = norm(sqrt(modulationM)-1+1i*floor(sqrt(modulationM)-1));                
            end
           
            
            %{
            if sqrt(modulationM)== round(sqrt(modulationM))
                scalingFactor = norm(sqrt(modulationM)-1+1i*floor(sqrt(modulationM)-1));
            else                
                scalingFactor = norm(sqrt(modulationM/2)-1+1i*floor(sqrt(modulationM/2)+1));
            end
            %}
            
            mapBitsID = qamdemod(datasymb(rxa,k)*scalingFactor,modulationM);
            
            receivedIDs(symbcount(rxa)) = mapBitsID;
            
            symb2BitNum = log2(modulationM);
            
            startIndex = (symbcount(rxa)-1)*symb2BitNum+1;
            endIndex = startIndex + (symb2BitNum-1);
            
            
            decodeBits(startIndex:endIndex) = de2bi(mapBitsID,log2(modulationM));
            
            receivedBitGroup = decodeBits(startIndex:endIndex);

            if modeID == 3 %for open loop mimo, mutiple streams for one node
                rxCompare = rxa;
            end

            transmitBitGroup = inb(startIndex:endIndex, rxCompare).';     
            numErr(rxa) = numErr(rxa) + sum(receivedBitGroup~=transmitBitGroup);
            symbcount(rxa) = symbcount(rxa) + 1;            
        end
        

    
        if (DEBUG_OUT)
            figure(230+rxa);
            
            plot(0, 0, 'ob', 'LineWidth', 2, 'MarkerSize', 10);
            hold on
            plot(real(datasymb(rxa,:)), imag(datasymb(rxa,:)), '+r', ...
                'LineWidth', 2, 'MarkerSize', 10);  
            hold off
            if ((modeID == 1) || (modeID ==2))
                stt = sprintf('Rx %d, Data bits, error %d, total %d\n', rxID, numErr(rxa), (symbcount(rxa)-1)*symb2BitNum);
            elseif modeID == 3
                stt = sprintf('Rx antenna %d, Data bits, error %d, total %d\n', rxa, numErr(rxa), (symbcount(rxa)-1)*symb2BitNum);
            end
            title(stt);
            pause(0.5);
        end
end % end for rxa=1:numRxAntenna
        
%% ====== end decoding ====== 
	end % end while (nsymbol < )

    
      

for rxa = 1:numRxAntenna
    %calculate BER
    totalBitNum = (symbcount(rxa)-1)*symb2BitNum;
    mBER(rxa) = numErr(rxa)/totalBitNum;
  
    %store the decoded symbols along with raw data
            for kk = 1:(symbcount(rxa)-1)
                
                startIndex = (kk-1)*symb2BitNum+1;
                endIndex = startIndex + (symb2BitNum-1);
                if modeID == 3 %for open loop mimo, mutiple streams for one node
                    rxCompare = rxa;
                end
                correspondingInput(kk) = qammod(bi2de(inb(startIndex:endIndex, rxCompare).'),modulationM);
                                    
                               
              
            end
            
    
    % --- calculate SINR --- 
  	pos_ct = zeros(48,1);
	neg_ct = zeros(48,1);
    
    SINR_subc_tmp = zeros(10,48);
    storeID = zeros(1,48);
    SINR_symbol = zeros(1,numOFDMsymb*48);
    
        for k = 1:numOFDMsymb*48
            
            subc = mod(k-1,48)+1;
            
            if useNewSINR            
                
                storeID(subc)=storeID(subc)+1;

                errorVector = correspondingInput(k)/scalingFactor-datasymb1(rxa,k);
                NPower = norm(errorVector)^2;                
                SPower = norm(correspondingInput(k)/scalingFactor).^2; %we normalize it to be this way

                SINR_symbol(k)=10*log10(SPower/NPower);
                SINR_subc_tmp(storeID(subc),subc) = SINR_symbol(k);
            
            else
            
                if (real(correspondingInput(k)) > 0)
                    %one more positive value
                    pos_ct(subc) = pos_ct(subc) + 1;
                    pos(pos_ct(subc), subc) = real(datasymb1(rxa,k));                
                else
                    %one more negative value
                    neg_ct(subc) = neg_ct(subc) + 1;
                    neg(neg_ct(subc), subc) = real(datasymb1(rxa,k));
                end
            
            end
        end
        
        if useNewSINR
            sinr_subc(rxa,:) = sum(SINR_subc_tmp)/numOFDMsymb;
        else
            for subc = 1:48
                    [~, ~, p] = find(pos(:,subc));
                    [~, ~, n] = find(-neg(:,subc));               
                    allres = [p; n];     %just putting vertical vectors together 
                    SP_subc=mean(allres.^2);
                    NP_subc=var(allres);
                    SINR = SP_subc/NP_subc;
                    sinr_subc(rxa,subc) = 10*log10(SINR);
            end
        end
      
        
        for subc = 1:48  
            thru_subc(subc) = mapSINRtoRate(sinr_subc(subc));
        end
        
       
     
         
        
        %when esitimating performance, exclude unused subcarrier
        finalSINR(rxa) = mean(sinr_subc(rxa,:));
   
        
        
        if useShannon
            avgThru(rxa) = 40e6/osamp*sum(log2(1+SINR));
        else
            avgThru(rxa) = mean(thru_subc);
        end 
    
    
        if modeID <= 2
            if pktType == 2
                    fprintf(1, '***** RXN %d, check with direction %d, error %d, total %d, BER %g, decoding SINR is %g\n', rxID, rxCompare, numErr, totalBitNum, mBER,finalSINR); 
            else
                fprintf(1, '***** RXN %d, error %d, total %d, BER %g, decoding SINR is %g\n', rxID, numErr, totalBitNum, mBER, finalSINR);         
            end
        else
                fprintf(1, '***** RX antenna %d, error %d, total %d, BER %g, decoding SINR is %g\n', rxa, numErr(rxa), totalBitNum, mBER(rxa), finalSINR(rxa));         
        end
 
end % end for rxa=1:numRxAntenna
    
    else
        fprintf(1, '***** RXN %d, polling packet, do not decode\n', rxID);
    
    end %end if DECODE_DATA
    
%% ------- end decoding OFDM data symbols ---------


%% ============= Visualize the results ==============
if (DEBUG_OUT) 
figure(200);
%plot(1:length(crosscorrOut), abs(crosscorrOut), '-^');
plot(1:length(selfcorrOut), abs(selfcorrOut), '-rx', ...
    1:length(energyOut), energyOut, '-bo', ...
    1:length(crosscorrOut), ...
    abs(crosscorrOut)/max(abs(crosscorrOut))*max(energyOut), '-g^', ...
    'LineWidth', 2);
end

%% ============= the channel matrix to feedback ==============
for txa = 1:numTxAntenna
    for rxa = 1:numRxAntenna   
        newCHout(rxID,txa,:) = [squeeze(LTF_CH(rxa,txa,33:64)).' squeeze(LTF_CH(rxa,txa,1:32)).'];
    end
end

if (numPkt == 1) break; end

end % end while ns < siglen








