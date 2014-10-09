%clc

clear all;



global numTxAntenna ...
       numRxAntenna LTF_up2 STF_up2 VLTF_t_up VLTF_freq newCHout rxData_board ...
       DECODE_DATA USESIM osamp numOFDMsymb sinr_subc otherTo48 useShannon usePilotCompen numRxNode pktType modulationM useNewSINR DEBUG_OUT modeID numRxAntenna rxData;
 
%% ====== Basic mode control 
USESIM = 1; % whether to use simulation
useSISO = 0; % use SISO transmission
modeID = 1; % 1: SU-MIMO; 2: MU-MIMO 3: open-loop MIMO
modulationM = 2; % QAM modulation mode, maximum supported 64

totalDataPkts = 3; % total data packets in experiment
DEBUG_OUT = 1; % show constellation or not
usePilotCompen = 1; % current way of pilot do not work for open loop MIMO

useFake = 1; % use a simple fake channel to debug

%% ====== PHY parameter configuration
TxDelay = 1000;
TxLength = 2^14-TxDelay; % Length of transmission. In [0:2^14-1-TxDelay]
CarrierChannel = 13; % Channel in the 2.4 GHz band. In [1:14] (avoid 1 to 11); 5GHz in [15:37]
TxGain_RF = 35;% Tx RF Gain. In [0:63] 
TxGain_BB = 1; % Tx Baseband Gain. In [0:3]
RxGain_BB = 10; % Rx Baseband Gain. In [0:31]
RxGain_RF = 1;%2; % Rx RF Gain. In [1:3]

%% ====== user selection related control

% select the user selection mode
useOSel = 0;
useRandom = 0;
useFix = 1;

% OPUS variations
randBKOF = 0; % use random backoff
useRBF = 0; % use random beamforming
useFair = 1; % use proportional fairness
pktInterval = 1; % how many packet for each CSI
usePartial = 1;  % user number selection              

%% ====== implicit controls


% the way to calculate SINR
if modulationM == 2
    useNewSINR = 0;
else
    useNewSINR = 1;
end

% determine receiver and antenna number
if useSISO
    numRxNode = 1; %the total number of receivers in the network
    numTxAntenna = 1; % number of tx antennas      
    numRxAntenna = 1; % number of rx antennas
else
    numTxAntenna = 2; % number of tx antennas
    if modeID == 1
        numRxNode = 1;
        numRxAntenna = 1; % number of rx antennas
    elseif modeID == 2
        numRxNode = numTxAntenna;
        numRxAntenna = 1; % number of rx antennas
    elseif modeID ==3
        numRxNode = 1;
        numRxAntenna = 2;
    end
end

% fixed user group
if useFix
    selectedUser = [1,2];
end

if numRxNode == 1
    rxID = 1;  
    selectedUser = 1;
end

%  determine the number of user to beamform to
if useOSel
    bfRxNode = 1; %the max number of Rx Nodes used to beamforming to
else
    if modeID==2
        bfRxNode = numTxAntenna; %the max number of Rx Nodes used to beamforming to
    else
        bfRxNode = 1;
    end
end


% the fake channel for test
% for MU-MIMO, you'd better add some crosstalk, or rx2 will not work
% because the LTF is lost if crosstalk is 0
if useFake 
   
    fakeCH = zeros(2, 2,64);
    size(fakeCH)
    for subc = 1:64
        fakeCH(:,:, subc) = [1 0.1;
                             0.1 1;];
        %}
        %0 entry will cause negative LTF estimated SINR, because the LTF
        %itself is not properly received
        %{
                fakeCH(:,:, subc) = [0 1;
                                     1 0;];
        %}
    end
end



% construction user distortion (to make the 100 user case work)
for rxn = 1:numRxNode
	for txa = 1:numTxAntenna
		for subc = 1:64
			tmpDistortion = rand + 1i*rand;
			tmpDistortion = tmpDistortion/norm(tmpDistortion);
			userDistortion(rxn,txa,subc)=tmpDistortion;
        end
    end
end

if USESIM osamp = 1; end

% ====== initialize
lastIsPolling = 0;
mDataAmount= zeros(1,numRxNode);
mTimeAmount = 0;
numPollingPkt = 0;
numDataPkt = 0;
colFlag = 0;
usePowerAllocation = 0;
pktBetweenSel = 0;
bkof_amount=zeros(1, totalDataPkts);
colLoc = []; %record the collision location
haveCSILast = 0; % used to support the peak reduction part

% ====== Read global variables 
MUMIMOglobalVars;

% ====== Less important control signals
numTxNode = 1;
txID = 1; % currently we only allow one transmitter with ID txID
pScale = 1; %USE this to control SINR 0.03 for low sinr
useFakeCol=0; % use the fake col to evaluate large number of usrs
probAllBeam = 0; % this do not work because pinv(a) is not orthogonal to null(a)
showAlignment = 0;
SU_focus = 1; %0: round robin 1: rx1 2:rx2
useShannon = 0;
distortForMany = 0; % test 100 usr case, distort every 20 traces
skipTrace = 30; %the first 30 traces are not used in trace driven simulation
FeedbackThreshold = 7; %the threshold for OPUS user to feedback 

% ====== implement QAM
% -------- Generate OFDM data symbol -----------
numSymbTotal = 48*numOFDMsymb;

symb2BitNum = log2(modulationM);
numBitsTotal = numSymbTotal*symb2BitNum;

    if modeID <= 2
        streamNum = numRxNode;
    elseif modeID == 3
        streamNum = numRxAntenna;
    end

fg = fopen('databits.dat', 'w');
for k = 1:numBitsTotal 
			for kk = 1:streamNum
					if (rand < 0.5)						 
						fprintf(fg, '1 ');   
                    else 
						fprintf(fg, '0 ');    
					end 
			end
	
	fprintf(fg, '\n');
end
fclose(fg);


inb = load('databits.dat');
numBitsInPkt = length(inb(:,1));
symbcount = 1;
posInSymb = 0;
posOut = 0;

while symbcount <= numSymbTotal

    posOut = posOut + 1;
    
    posInSymb = posInSymb + 1;
    if (posInSymb == 65)
        posInSymb = 1;
    end
     
    % Insert pilot
    if (posInSymb == (7+1) || posInSymb == (21+1) ...
        || posInSymb == (64-21+1) || posInSymb == (64-7+1))
        for rxn = 1:streamNum
            symbIn(posOut, rxn) = 1;
        end
        continue;
    end
    
     
    % unused subc
    % zero-freq in index1, minus-freq near half-sequence
    if (posInSymb==0+1 || (posInSymb > 26+1 && posInSymb < 64-26+1))
        
        for rxn = 1:streamNum
            symbIn(posOut, rxn) = 0;
        end 
      
        continue;
    end 

    for rxn = 1:streamNum
        
        %map bits to ID here
            startIndex = (symbcount-1)*symb2BitNum+1;
            endIndex = startIndex + (symb2BitNum-1);
            
        idOnConstellation = bi2de(inb(startIndex:endIndex,rxn).');       
        %do not normalize for higher level modulation!
        symbIn(posOut, rxn) = qammod(idOnConstellation,modulationM); 
        
    end     
 symbcount=symbcount+1;
end
 

% ========== Initializaton of WARP radio parameters ========
if USESIM == 0

USE_AGC = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the WARPLab experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create a vector of node objects
nodes = wl_initNodes(numRxNode + numTxNode);
node_tx = nodes(1);
node_rx = nodes(2:length(nodes));
%{
nodes = wl_initNodes(3);
node_tx = nodes(2);
node_rx = nodes(3);
%}

%Create a UDP broadcast trigger and tell each node to be ready for it
eth_trig = wl_trigger_eth_udp_broadcast;
wl_triggerManagerCmd(nodes,'add_ethernet_trigger',[eth_trig]);

%Get IDs for the interfaces on the boards. Since this example assumes each
%board has the same interface capabilities, we only need to get the IDs
%from one of the boards
[RFA,RFB] = wl_getInterfaceIDs(nodes(1));
RF_vector = [RFA,RFB];

%Set up the interface for the experiment
wl_interfaceCmd(nodes,sum(RF_vector),'tx_gains',TxGain_BB,TxGain_RF);
wl_interfaceCmd(nodes,sum(RF_vector),'channel',2.4,CarrierChannel);

if(USE_AGC)
    wl_interfaceCmd(nodes,sum(RF_vector),'rx_gain_mode','automatic');
    wl_basebandCmd(nodes,'agc_target',-10);
    wl_basebandCmd(nodes,'agc_trig_delay', 500);
    wl_basebandCmd(nodes,'agc_dco', true);
else
    wl_interfaceCmd(nodes,sum(RF_vector),'rx_gain_mode','manual');
    RxGainRF = 1; %Rx RF Gain in [1:3]
    RxGainBB = 15; %Rx Baseband Gain in [0:31]
    wl_interfaceCmd(nodes,sum(RF_vector),'rx_gains',RxGain_RF,RxGain_BB);
end


%We'll use the transmitter's I/Q buffer size to determine how long our
%transmission can be
%txLength = nodes(1).baseband.txIQLen;

%Set up the baseband for the experiment
wl_basebandCmd(nodes,'tx_delay',TxDelay);
wl_basebandCmd(nodes,'tx_length',TxLength);


end % end if USESIM==0

% ================= Start transmission rounds =====================

%initialize before the while loop
numPrecodeRounds = 0;
afterData = 1;

 while numDataPkt < totalDataPkts

     
     clear precodeM;
     clear iPrecodeM;
     
    
    
 numPrecodeRounds = numPrecodeRounds + 1;    
fprintf(1, '==== precodeRound %d \n', numPrecodeRounds);

%%  determine the type of the packet
% type1: polling
% type2: polling_prob
% type3: data

% polling packet
if modeID <=2
if afterData && (pktBetweenSel==0)
    pktType = 1;
    
    lastIsPolling = 1;
    
    fprintf(1, '==== precodeRound %d transmits polling packet\n', numPrecodeRounds);
    
    
    DECODE_DATA = 0;
    

    
    numPollingPkt = numPollingPkt + 1; %update the number of polling pkt
  
    
    afterData = 0;
    
    % reset the number of selected users who do not need to feedback CSI
    % because their channel aligment are good enough
    noCSINum = 0;

   
    if useOSel
        
        bfRxNode = 1; 
        
        %select the core

	 if useFair && numDataPkt>=1
	     %calculate average throughput for proportional fairness
	     alpha = 0.1;
	     for rxn = 1:numRxNode
             if numDataPkt == pktInterval %just after the first round
                PF_Thru(rxn) = 10e6;
             else
                PF_Thru(rxn) = (1 - alpha) * PF_Thru(rxn) + alpha*thruRecord(numDataPkt-1,rxn);
             end
	     end
	     selMtric = 1./(PF_Thru/10e6);
	     % select the core user, the one with higher selMtric has a higher probability to be selected
		 selectedUser = randsample(1:numRxNode, 1, true, selMtric)
	 else
		selectedUser = randi([1, numRxNode])
	 end


        userToAdd = selectedUser;

               
        numUserLeft = numTxAntenna - 1;
    else
        numUserLeft = 0;
    end

    
else


    % polling_prob packet   
    if numUserLeft > 0 

        pktType = 2;
        
        fprintf(1, '==== precodeRound %d transmits polling_prob packet\n', numPrecodeRounds);
        
        
        DECODE_DATA = 1;
        
        %after one polling_prob packet, we add one more user to the group
        numUserLeft = numUserLeft - 1;

    else %data packet
        
        pktType = 3;
        
        fprintf(1, '==== precodeRound %d transmits DATA packet\n', numPrecodeRounds);
        
        %the first one
        if afterData == 0
            pktBetweenSel = pktInterval;
        end


        pktBetweenSel = pktBetweenSel - 1;

        
        DECODE_DATA = 1;

        numDataPkt = numDataPkt + 1; %update the number of data pkt


        afterData = 1;

    end

end % end if afterData && (pktBetweenSel==0)

elseif modeID == 3 % open loop MIMO
        pktType = 3;
        fprintf(1, '==== precodeRound %d transmits DATA packet\n', numPrecodeRounds);
        DECODE_DATA = 1;
        numDataPkt = numDataPkt + 1; %update the number of data pkt
end % end if modeID <=2

% randomly select user only for data pkt
if useRandom && pktType == 3
    selectedUser = randsample(numRxNode,bfRxNode)
end

if pktType == 3
    if length(selectedUser)==numTxAntenna
        selectRecord(numDataPkt,:)=selectedUser;
    else
        selectRecord(numDataPkt,:)=[selectedUser zeros(1,numTxAntenna-length(selectedUser))];
    end
end



    %CSI cutting unreasonable points
    if USESIM && ((modeID==1)||(modeID == 2))
        if useOSel
            if pktType == 3
                CSI_last = CSI;
                haveCSILast = 1;
            end
        end

        if useFix || useRandom
            if numDataPkt>1
                CSI_last = CSI;
                haveCSILast = 1;
            end
        end
    end


if (pktType ~= 1) && lastIsPolling
    
    lastIsPolling = 0;
    

    % update channel
     disp('CSI is updated');
    CSI = newCHout;
    
         % compute CQI
         for rxn = 1:numRxNode        
             for subc = 7:59
                 if subc == 33 continue; end
                CQI_tmp(rxn,subc)= CSI(rxn,:,subc)*CSI(rxn,:,subc)';
             end
             CQI(rxn)=sum(CQI_tmp(rxn,:))/52;
             %realCQI(rxn)=sum(realCQI_tmp(rxn,:))/52;
            %FIXME: cut off the unreasonable peaks
             if CQI(rxn)>100 && USESIM && haveCSILast
                 for subc = 7:59
                     if subc == 33 continue; end
                     CSI(rxn,:,subc)=CSI_last(rxn,:,subc);
                    CQI_tmp(rxn,subc)= CSI(rxn,:,subc)*CSI(rxn,:,subc)';
                 end
                 CQI(rxn)=sum(CQI_tmp(rxn,:))/52;
             end

         end

         if pktType == 3
             CQI_store(numDataPkt,:) = CQI;
         end
         
     
     % compute CDI
      for rxn = 1:numRxNode        
             for subc = 7:59
                 if subc == 33 continue; end
                    CDI(rxn,:,subc)= CSI(rxn,:,subc)/sqrt(CQI_tmp(rxn,subc));
             end            
      end
      %{
      % test CDI mag
         for rxn = 1:numRxNode        
             for subc = 7:59
                 if subc == 33 continue; end
                mCDI_tmp(rxn,subc)= CDI(rxn,:,subc)*CDI(rxn,:,subc)';
             end
             mCDI(rxn)=sum(mCDI_tmp(rxn,:))/52;
         end
         
         mCDI
      %}
      
      %plot channel phase difference between ant 1 and 2
      %phaseDiff(numDataPkt) = angle(CDI(1,1,10))-angle(CDI(2,1,10));
      %phaseDiff(numDataPkt) = angle(CDI(1,1,10)/CDI(1,2,10));      
      %phaseDiff(numDataPkt) = sum(angle(CDI(1,1,:))-angle(CDI(1,2,:)))/52;
      
end




%force it to be SU-MIMO/MU-MIMO
if modeID == 1
    numRxNow = 1;
elseif modeID == 2
    numRxNow = bfRxNode;
elseif modeID == 3
    numRxNow = 1;
end

%if SU-MIMO, round-robin user selection
    if (modeID == 1)
        if SU_focus == 0               
            %ceil(numPrecodeRounds/2)
            if (mod(ceil(numPrecodeRounds/2),numRxNode)==0)
                rxID = numRxNode;
            else
                rxID = mod(ceil(numPrecodeRounds/2),numRxNode);
            end  
       
        else
            rxID = SU_focus;        
        end 
    end


%% ====== precoding

clear DATA_up2;

% not polling packet
if pktType ~= 1
    
    %random beamforming, the core is randmly generated
    if useRBF
        for subc = 7:59
            if subc==33 continue; end
            if length(selectedUser)==1
                CSI(selectedUser,:,subc)=rand(1,numTxAntenna)+1i*rand(1,numTxAntenna);
            end
        end
    end


% compute weights for each subcarrier

%predefine the metrix to make SU-MIMO and MU-MIMO compatiable

if modeID == 1
    iPrecodeM = zeros(numTxAntenna,64);    
elseif modeID == 2
    %TODO: initialize for OPUS
    if pktType == 3
        iPrecodeM = zeros(numTxAntenna, bfRxNode,64);   
    elseif pktType == 2
        if probAllBeam
            targetSize = numTxAntenna;
        else
            targetSize = numTxAntenna - length(selectedUser);
        end
        iPrecodeM = zeros(numTxAntenna, targetSize,64);   
    end
end




if (modeID == 1)||(modeID == 2)
    % determine the power scale for subcarriers including pilot
     for subc = 7:59 
             if (subc == 33) continue; end
             for rxn = 1:bfRxNode        
                chGain(rxn,subc)= CSI(rxn,:,subc)*CSI(rxn,:,subc)';
             end
             
             for rxn = 1:bfRxNode
                 
                 if usePowerAllocation
                    PA_scale(rxn,subc) = chGain(rxn,subc)/sum(chGain(:,subc));
                 else
                    PA_scale(rxn,subc) = 1/bfRxNode;
                 end
                 
                    
                    %chGain_record(rxn,subc,numValidDataPkt)= chGain(rxn,subc);             
                    %PA_record(rxn,subc,numValidDataPkt)=PA_scale(rxn,subc);
               
                
             end
             
            
     end
 end



 %precoding weight matrix
for subc = 7:59
    
    if (subc == 33) continue; end
 

   
    %if (isPrecodingOn == 1)
        
        if (modeID == 1 && numTxAntenna > 1)%SU-MIMO beamforming

                iPrecodeM(:,subc) = exp(-1i*angle(CSI(rxID,:,subc)));
        end

        if (modeID == 2 && numTxAntenna > 1) % MU-MIMO
            
            %Moore-Penrose pseudoinverse            
            %iPrecodeM(:,:,subc) = CDI(selectedUser,:,subc)' * inv(CDI(selectedUser,:,subc) * CDI(selectedUser,:,subc)');

            %TODO: beamform to optimal directions (these directions are orthogonal with each other) 
            %since these optimal directions are semi-orthogonal matrix, its pinv is its transpose
            if pktType == 3
                iPrecodeM(:,:,subc) = pinv(CSI(selectedUser,:,subc));
            elseif pktType == 2
                if probAllBeam
                    probDirection = null(CSI(selectedUser,:,subc)).';
                    CSI_fake(:,:,subc) = [CSI(selectedUser,:,subc); probDirection ];
                    iPrecodeM(:,:,subc) = pinv(CSI_fake(:,:,subc));
                else
                    iPrecodeM(:,:,subc) = null(CSI(selectedUser,:,subc));
                end
            end
            
    
        end  

        if (modeID == 3)
            iPrecodeM(:,:,subc) = ones(numRxAntenna,numTxAntenna);
        end

        
    %else	
        %{
        %disable precoding        
        iPrecodeM(:,:,subc) = ones(bfRxNode,numTxAntenna);
        %}
    %end

end


%{
Power Allocation
    for subc = 7:59        
         for rxn = 1:bfRxNode 
             
             if (subc == 33) continue; end
             
             if modeID == 1 %SU-MIMO
                
            precodeM(:,rxn,subc) = iPrecodeM(:,rxn,subc)*PA_scale(rxn,subc);  
          
         end        
         end
%}

    precodeM = iPrecodeM;

if modeID ==2 && pktType == 3 
    CNum_w(:,numDataPkt)=zeros(1,64); 
    for subc = 7:59
            if (subc==0+33)||(subc==-21+33)||(subc==-7+33)||(subc==7+33)||(subc==21+33)            
                continue;
            end        
                CNum_w(subc,numDataPkt)=20*log10(cond(precodeM(:,:,subc)));

    end   
end

% ------ construct VLTF
VLTF_coded = zeros(64,numTxAntenna);
for subc = 1:64
       
    for txa = 1:numTxAntenna   
        
        if modeID == 2 %MU-MIMO		
            if pktType == 3
                VLTF_coded(subc, txa) = precodeM(txa,:,subc) * VLTF_freq(1:bfRxNode,subc);	
            elseif pktType == 2
                VLTF_coded(subc, txa) = precodeM(txa,:,subc) * VLTF_freq(1:targetSize,subc);	
            end
        elseif modeID == 1 %SU-MIMO
            VLTF_coded(subc, txa) = precodeM(txa,subc) * VLTF_freq(rxID,subc);	
        end
		
        if (isnan(VLTF_coded(subc, txa))) 
            VLTF_coded(subc, txa) = 0;
        end
    end 
end 

% generate time domain VLTF
maxamp=zeros(numTxAntenna,1);
for txa = 1:numTxAntenna
    vt = ifft(fftshift(VLTF_coded(:,txa)), 64);
    VLTF_t(txa, :) = [vt(64-16+1:64); vt].';
    VLTF_t_up(txa, :) = interp(VLTF_t(txa, :), osamp);
	maxamp(txa) = max([ max(abs(real(VLTF_t_up(txa,:)))), max(abs(imag(VLTF_t_up(txa,:))))]);
end

% VLTF Normalization
if modeID == 2
    %MU-MIMO must use total normalization to 
    VLTF_t_up = VLTF_t_up/max(maxamp);  
elseif modeID == 1
    %SU-MIMO can have antenna level normalization
    for txa = 1:numTxAntenna   
        VLTF_t_up(txa,:) = VLTF_t_up(txa,:) / maxamp(txa); 
    end
end


% *** Zero-forcing precoding for data in each subcarrier ***

for os = 1:numOFDMsymb
    
    for subc = 1:64 
            
        %line mode ID
        sind = mod(subc+31,64)+1;
        
        %folded mode ID
        bindex = (os-1)*64+subc; %data bit index
        for txa = 1:numTxAntenna   
            if numTxAntenna>1
                %FIXME: map s to subc
                if modeID == 2 %MU-MIMO (different user's data are combined here)                   
                    %freqDataBlock(bindex, txa) = precodeM(txa,:,sind) * symbIn(bindex,:).';
                    
                     % TODO: control the real precoding part
                     
                     if pktType == 3
                         usrSymb = precodeM(txa,:,sind) .* symbIn(bindex,selectedUser);
                     elseif pktType == 2
                         usrSymb = precodeM(txa,:,sind) .* symbIn(bindex,1:targetSize);
                     end

                     freqDataBlock(bindex, txa) = sum(usrSymb);  
                     %mRatioMatrix = abs(precodeM(txa,1,sind))/abs(precodeM(txa,2,sind));
                     %mRatio = abs(usrSymb(1))/abs(usrSymb(2));
                elseif modeID == 1            %SU-MIMO                     
                    freqDataBlock(bindex, txa) = precodeM(txa,sind) * symbIn(bindex,rxID).';
                elseif modeID == 3            %open loop MIMO                     
                    freqDataBlock(bindex, txa) = symbIn(bindex,txa).';
                end    

            else %SISO
                
                freqDataBlock(bindex, txa) = symbIn(bindex,rxID).';
            end
                
        end 
        
       
    end
end
    
% *** IFFT: generate time domain data signal 
for txa = 1:numTxAntenna
    outCount = 0;
    
    for ds = 1:numOFDMsymb
        ifftOut = ifft(freqDataBlock((ds-1)*64+1:ds*64,txa), 64);	
        % *** Add guard interval 
		output_t(outCount+1:outCount+80, txa) = [ifftOut(64-16+1:64); ifftOut];
		outCount = outCount + 80;
    end
     
	DATA_up2(:,txa) = interp(output_t(:,txa), osamp); % Upsample
    
    maxamp1(txa) = max([max(abs(real(DATA_up2(:, txa)))), max(abs(imag(DATA_up2(:, txa))))]);
    
    txa_ratio_real(txa) = max(abs(real(DATA_up2(:, txa))))/min(abs(real(DATA_up2(:, txa))));
    txa_ratio_imag(txa) = max(abs(imag(DATA_up2(:, txa))))/min(abs(imag(DATA_up2(:, txa))));
    
    txa_ratio_power(txa) = max(abs(DATA_up2(:, txa)))/min(abs(DATA_up2(:, txa)));
    
end

if pktType == 3


    normW(numDataPkt) = 0;

    if modeID == 2
    %calculate the norm of W
        for rxn = 1:bfRxNode
            for subc = 7:59
                if subc == 33 continue; end

                normW_subc(subc) = norm(precodeM(:,rxn,subc));

            end
            normW_aver = sum(normW_subc)/52;
            normW(numDataPkt) = normW(numDataPkt) + normW_aver;
        end
    end
end


if pktType ~= 1
    den(numPrecodeRounds) = max(maxamp1);
end
    
  
%Data Normalization
if modeID == 2
    %MU-MIMO must use total normalization
    DATA_up2 = pScale*DATA_up2 / max(maxamp1);  
elseif modeID == 1
    %SU-MIMO can have antenna level normalization
    for txa = 1:numTxAntenna   
        DATA_up2(:,txa) = pScale*DATA_up2(:,txa) / maxamp1(txa); 
    end
elseif modeID == 3
    %open loop MIMO also use total normalization?
    DATA_up2 = pScale*DATA_up2 / max(maxamp1);  
end


else %Polling packet does not have data
    for txa = 1:numTxAntenna
        DATA_up2(:,txa) = zeros(length(LTF_up2), 1); 
    end
end 


% ------- Concatenate all tx samples ---------
zero_vector = zeros(1, length(LTF_up2));

clear TxData;
txn = 1;
if numTxAntenna == 1
td = [STF_up2, LTF_up2, VLTF_t_up(1,:), DATA_up2(:,1).'];
elseif numTxAntenna == 2
td = [STF_up2, LTF_up2, zero_vector, VLTF_t_up(1,:), DATA_up2(:,1).'];
elseif numTxAntenna == 3
td = [STF_up2, LTF_up2, zero_vector, zero_vector, VLTF_t_up(1,:), DATA_up2(:,1).'];
elseif numTxAntenna == 4
td = [STF_up2, LTF_up2, zero_vector,zero_vector, zero_vector, VLTF_t_up(1,:), DATA_up2(:,1).'];
end

if (TxLength-length(td) < 0)
	fprintf(1, 'ERROR!!! data overflow!\n');
	return;
end

if (USESIM==0)
    TxData(:,1,txn) = [td zeros(1,max(1,TxLength-length(td)))];
else
    TxData(:,1,txn) = [1e-5*zeros(1,80) td zeros(1,max(1,TxLength-length(td)))];    
end


td = [];
if (numTxAntenna >= 2)
if numTxAntenna == 2
td = [STF_up2, zero_vector, LTF_up2, VLTF_t_up(2,:), DATA_up2(:,2).'];
elseif numTxAntenna == 3
td = [STF_up2, zero_vector, LTF_up2, zero_vector, VLTF_t_up(2,:), DATA_up2(:,2).'];
elseif numTxAntenna == 4
td = [STF_up2, zero_vector, LTF_up2, zero_vector, zero_vector, VLTF_t_up(2,:), DATA_up2(:,2).'];
end
end % end txantenna 2
if (USESIM==0)
    TxData(:,2,txn) = [td zeros(1,max(1,TxLength-length(td)))];
else
    TxData(:,2,txn) = [1e-5*zeros(1,80) td zeros(1,max(1,TxLength-length(td)))];
    %TxData(:,2,txn) = [1e-5*zeros(1,80) td];
end

% tx antenna 3
td = [];
if (numTxAntenna >= 3)
if (numTxAntenna==3)
td = [STF_up2, zero_vector, zero_vector, LTF_up2, VLTF_t_up(3,:), DATA_up2(:,3).'];
elseif (numTxAntenna==4)
td = [STF_up2, zero_vector, zero_vector, LTF_up2, zero_vector, VLTF_t_up(3,:), DATA_up2(:,3).'];
end
if (USESIM==0)
    TxData(:,3,txn) = [td zeros(1,max(1,TxLength-length(td)))];
else
    TxData(:,3,txn) = [1e-5*zeros(1,80) td zeros(1,max(1,TxLength-length(td)))];
    %TxData(:,3,txn) = [1e-5*zeros(1,80) td];
end
end

% tx antenna 4
td = [];
if (numTxAntenna >= 4)
td = [STF_up2, zero_vector, zero_vector, zero_vector, LTF_up2, VLTF_t_up(4,:), DATA_up2(:,4).'];
if (USESIM==0)
    TxData(:,4,txn) = [td zeros(1,max(1,TxLength-length(td)))];
else
    TxData(:,4,txn) = [1e-5*zeros(1,80) td zeros(1,max(1,TxLength-length(td)))];
    %TxData(:,4,txn) = [1e-5*zeros(1,80) td];
end
end



if USESIM == 1       
	
	%% start: add this part for trace-driven simulation
%per subcarrier channel distortion
clear rxData;

%polling packet, polling_prob packet and the following data packet should share the same CSI
%information (be distorted by the same channel)


if pktType ~= 3 %data packet will naturall + 1
   CH_Sequence_Index = numDataPkt + skipTrace + 1;
else
   CH_Sequence_Index = numDataPkt + skipTrace;
end

for rxn = 1:streamNum    

	for txa = 1:numTxAntenna
		CH_cond = zeros(1,64); 
        
          %load the files     
          %t1 stands for transmit antenna 1 and 2
          %t2 stands for transmit antenna 3 and 4

	  %map Rx index
	  mapRxIndex = mod(rxn,20);
	  if mapRxIndex == 0
		  mapRxIndex = 20;
	  end

      if useFake
          for subc = 1:64
              CH_cond(subc) = fakeCH(rxn,txa,subc);
          end
      else         
            %channel_file = load(sprintf('OFDM_CHGaintrace_r%d_t%d.txt',rxn,ceil(txa/2)));
            channel_file = load(sprintf('OFDM_CHGaintrace_r%d_t%d.txt',mapRxIndex,ceil(txa/2)));
	    tmpSize = size(channel_file);
	   
		%get collected channel data
		for subc = 1:64
			%which row in the file of the corresponding subcarrier
		    rowIndex = (CH_Sequence_Index - 1) * numTxAntenna*64 + (1-mod(txa,2))*64 + subc;
		    %wrap around when the execution exceeds the amount of traces
		    rowIndex = mod(rowIndex,tmpSize(1));
		    if rowIndex == 0
			    rowIndex = tmpSize(1);
		    end

            CH_cond(subc) = channel_file(rowIndex,3) + 1i*channel_file(rowIndex,4);

                 if distortForMany
                             %add user distortion to differentiate users
                             CH_cond(subc) =   CH_cond(subc)*userDistortion(rxn,txa,subc); 
                 end
   
        end 
      end
        
		%Antenna data before distortion
		data_before = squeeze(TxData(:,txa,txID));
		%channel distortion by convolution
        
        %copy this part to all other folders
		% Reverse engineering STF, distort, remodulate
		%data_after = conv(data_before, ifft(fftshift(CH_cond)));
		data_after = zeros(length(data_before)+length(CH_cond)-1,1);
		%data_after(160-63:160) = (fft(data_before(160-63:160))).'.*fftshift(CH_cond);
        data_after(160-63:160) = (fft(data_before(160-63:160))).'.*ones(size(CH_cond));% test 0311
		data_after(160-63:160) = ifft(data_after(160-63:160));
		data_after(160-79:160-64) = data_after(160-15:160);
		%data_after(240-63:240) = (fft(data_before(240-63:240))).'.*fftshift(CH_cond);
        data_after(240-63:240) = (fft(data_before(240-63:240))).'.*ones(size(CH_cond));% test 0311
		data_after(240-63:240) = ifft(data_after(240-63:240));
		data_after(240-79:240-64) = data_after(240-15:240); 
		%scale = 1/max([max(real(data_after(81:240))), max(imag(81:240))]);
		%data_after(81:240) = scale * data_after(81:240); 
        
		% Reverse engineering LTF, distort, remodulate        
		for txaa=1:numTxAntenna
            data_after(240+(txaa-1)*160+96-63:240+(txaa-1)*160+96) = (fft(data_before(240+(txaa-1)*160+96-63:240+(txaa-1)*160+96))).'.*fftshift(CH_cond);
            data_after(240+(txaa-1)*160+96-63:240+(txaa-1)*160+96) = ifft(data_after(240+(txaa-1)*160+96-63:240+(txaa-1)*160+96)); 
            data_after(240+(txaa-1)*160+96-64-31:240+96-64) = data_after(240+(txaa-1)*160+96-31:240+96); 
            data_after(240+(txaa-1)*160+96+1:240+(txaa-1)*160+96+64) = data_after(240+(txaa-1)*160+96-63:240+(txaa-1)*160+96); 
		%scale = 1/max([max(real(data_after(240+(txaa-1)*160+1:240+(txaa-1)*160+160))), max(imag(data_after(240+(txaa-1)*160+1:240+(txaa-1)*160+160)))]);
		%data_after(240+(txaa-1)*160+1:240+(txaa-1)*160+160) = scale * data_after(240+(txaa-1)*160+1:240+(txaa-1)*160+160); 
		%fprintf(1, 'scale=%g\n', scale);
		end
		tmp = 80+160+160*numTxAntenna;
        
		% Reverse engineering data, distort, remodulate
		for d = tmp+16:80:length(data_before)-79
			db = ifft(fft(data_before(d+1:d+64)).'.*fftshift(CH_cond));
			data_after(d+1:d+64) = db;
			data_after(d-15:d) = db(49:64); 
			%data_after(d-15:d+64) = scale * data_after(d-15:d+64);
		end
		
		%initialize rx data
		if (txa==1) rxData(rxn,:) = zeros(1,length(data_before)); end
        
		%antenna data after distortion
		%the second index of TxData is antenna ID		
		rxData(rxn,:) = rxData(rxn,:) + data_after(1:length(data_before)).';
	end %end for txa = 1:numTxAntenna
	
    
        sigma = 1e-6;
        
    

	% Noise floor depends on receiver
        rxData(rxn,:) = rxData(rxn,:) + ...
            normrnd(0, sqrt(sigma), [1 length(rxData(rxn,:))])/1.414...
            + 1i*normrnd(0, sqrt(sigma), [1 length(rxData(rxn,:))])/1.414;% Add Gaussian noise 
    
    
	
	%store to file
	% generate "standard" format for complex numbers [real imag, real imag,...]
	fname = sprintf('tmp_sim_modulated_Rx%d.txt', rxn);
    fo = fopen(fname, 'wb'); 
	o1 = [];
	a = real(rxData(rxn,:));
	b = imag(rxData(rxn,:));
	for k = 1:length(rxData(rxn,:))
		o1 = [o1 a(k) b(k)];
	end
	fwrite(fo, o1, 'float32');
	fclose(fo);

end %end for rxn = 1:numRxNow
%% end: add this part for trace-driven simulation


else % USESIM=0, run on WARP

%% ============= Start tx and rx =======================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit and receive signal using WARPLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign the data to transmit
for txn = 1:numTxNode
        wl_basebandCmd(node_tx(txn),RF_vector, 'write_IQ', TxData(:,:,txn));
end

% start transmission and reception
wl_interfaceCmd(node_tx,sum(RF_vector),'tx_en');
wl_interfaceCmd(node_rx,sum(RF_vector),'rx_en');

wl_basebandCmd(node_tx,sum(RF_vector),'tx_buff_en');
wl_basebandCmd(node_rx,sum(RF_vector),'rx_buff_en');

eth_trig.send();

% get the received data
for rxn = 1:numRxNode
        rxData_board(:,:,rxn) = wl_basebandCmd(node_rx(rxn),RF_vector(1:numRxAntenna),'read_IQ', 0, TxLength+TxDelay);
end

end % end if USESIM==1

% ============ Performance Record ============
if (modeID == 1 || modeID == 2)  
Subc_SINR=zeros(bfRxNode,64);
else
Subc_SINR=zeros(numRxAntenna,64);
end

if pktType == 2
    %note that the size of the matrix is always changing
    overhearSINR=zeros(numRxNode,length(selectedUser));
end

%All the receivers should feedback their CSI

for rxn = 1:numRxNode 
    
        if pktType == 1      %polling packet  
            
            [~, ~] = MUMIMOrx(rxn,rxn);
            
        elseif pktType == 2  %polling_prob packet, overhear
            
                if any(selectedUser==rxn)==0 % for those unselected users


                    fprintf(1, '***** RXN %d not selected \n', rxn);         

                    %TODO: overhear part for PUS
                    if probAllBeam
                        startID = length(selectedUser)+1;
                    else
                        startID = 1;
                    end

                    for dirID = startID:targetSize
                    %for those unselected users,check correlation with opt direction
                    [decodingSINR, avgThru] = MUMIMOrx(rxn,dirID);
                                %store the overheared SINR
                                overhearSINR(rxn,dirID)=decodingSINR;
                                %overhearSINR(rxn,dirID)=decodingSINR-10*log10(CQI(rxn));                                
                    end
                    % also record the optimal direction index
                    [avgOverhearSINR(rxn),optDirID(rxn)]=max(overhearSINR(rxn,:));

                else

                    fprintf(1, '***** RXN %d selected \n', rxn);         
                    %selected users will not be added to the group
                    avgOverhearSINR(rxn)=0;

                end
            
        else            %data packet
            [decodingSINR, avgThru] = MUMIMOrx(rxn,rxn);
        end
        
        % only record performance for data packet
        if pktType == 3 
            
            if modeID == 2 %MUMIMO
                if find(selectedUser == rxn)>0
                    pktSINR(numDataPkt,rxn) = decodingSINR;
                    thruRecord(numDataPkt,rxn) = avgThru; 
                else
                    pktSINR(numDataPkt,rxn) = 0;
                    thruRecord(numDataPkt,rxn) = 0; 
                end
            elseif modeID == 1 %SU-MIMO
                if rxn == rxID
                    pktSINR(numDataPkt,rxn) = decodingSINR;
                    thruRecord(numDataPkt,rxn) = avgThru; 
                else
                    pktSINR(numDataPkt,rxn) = 0;
                    thruRecord(numDataPkt,rxn) = 0; 
                end
            elseif modeID == 3 %open loop MIMO
                    pktSINR(numDataPkt,:) = decodingSINR;
                    thruRecord(numDataPkt,:) = avgThru; 
            end

                
                for subc = 7:59    
                    if (subc == 33) ||(subc==-21+33)||(subc==-7+33)||(subc==7+33)||(subc==21+33)
                        continue;
                    end
                    if (modeID == 1)||(modeID == 2)
                        Subc_SINR(rxn,subc)=sinr_subc(otherTo48(subc));  
                    elseif modeID == 3
                        Subc_SINR(:,subc)=sinr_subc(:,otherTo48(subc));  
                    end
                end            


        end 
end




% for the probing probe packet, select one more user
if pktType == 2 && useOSel
    
    
    % in the one users case, the user's channel magnitude can mask the
    % correlation, use a scaling factor to control
        preferMetric = avgOverhearSINR * targetSize/numTxAntenna
    % the one with the largest overhearing SINR is the most parallel to one of the optimal directions
    if randBKOF == 0

        % the time domain enery burst mechanism
        preferBits = 8;
        bitOverhead = 3;
        totalLevels = 2^preferBits;
        %preferMax = 45;
        %pResolution = preferMax/totalLevels
        pResolution = max(preferMetric)/totalLevels
        preferMetric_q = round(preferMetric/pResolution)
        bkof_amount(numDataPkt+1)=preferBits*bitOverhead;
        
        %collision
                if length(find(preferMetric_q==max(preferMetric_q)))>1

                    colLoc = [colLoc,numPrecodeRounds];
                    %numUserLeft = 0;
                    colFlag = 1;

                end
                
                %TODO: beamform to part
               
                if(max(preferMetric)<FeedbackThreshold)&&usePartial
                    %beamform to part of the users
                    fprintf(1, '!!! all users left have bad alignment\n');
                    %jump to beamforming
                    %overhearSINR(numPrecodeRounds)=max(preferMetric);
                    numUserLeft = 0;
                else
                    %add one more user


				 if useFair && numDataPkt>=1
				     %calculate average throughput for proportional fairness
				     alpha = 0.1;
				     for rxn = 1:numRxNode
					 if numDataPkt == 1
					    %PF_Thru(rxn) = 10e6;
					    PF_Thru(rxn) = 35e6;
					 else
					    PF_Thru(rxn) = (1 - alpha) * PF_Thru(rxn) + alpha*thruRecord(numDataPkt-1,rxn);
					 end
				     end
				     preferMetric = preferMetric./(PF_Thru/10e6);
				 end



                    [~,userToAdd] = max(preferMetric);    
                    bfRxNode = bfRxNode + 1;
                    selectedUser = [selectedUser, userToAdd]
                end




                %random beamforming, no feedback               
                if useRBF
                    noCSINum = length(selectedUser);

                    %TODO: replace CSI with its fake direction ID
                    for subc = 7:59
                        if subc==33 continue; end;
                        CSI(userToAdd,:,subc)=iPrecodeM(:,optDirID(userToAdd),subc).';     
                    end

                    fprintf(1, '------ Selected RX %d has good alignment with direction %d, no CSI feedback\n', userToAdd, optDirID(userToAdd)); 
                end
    else

        % TODO: implement the random backoff mechanism
        
        alignThreshold = 4;
        feedbackGroup = find(preferMetric>alignThreshold);


        if ~isempty(feedbackGroup)


            % implement the priority random backoff
            bkofVal = [];
            CWmax = 60;
            for ita = 1:length(feedbackGroup)
                CW = 2*max(CWmax-round(preferMetric(feedbackGroup(ita))),1);
                %CW = 32;
                bkofVal(ita) = randi([1 CW]);
            end




            %check for collision
            if useFakeCol
                % use the fake collision evaluation for a large number of users
                totalFakeNum = 100;
                goodRatio = 0.1;
                P_col = 1-(1-1/(CWmax/2*2))^(totalFakeNum*goodRatio-1)
                if rand<P_col

                    colLoc = [colLoc,numPrecodeRounds];
                    numUserLeft = 0;
                    colFlag = 1;

                end
            else
                if length(find(bkofVal==min(bkofVal)))>1

                    colLoc = [colLoc,numPrecodeRounds];
                    % if collision happened, the BS beamforms to already selected users
                    numUserLeft = 0;
                    colFlag = 1;

                end
            end

            if colFlag == 0
                %userToAdd = feedbackGroup(randi([1 length(feedbackGroup)]));
                [bkof_incre,userInGroup] = min(bkofVal);    
                userToAdd = feedbackGroup(userInGroup);
                % store the bkof performance
                bkof_amount(numDataPkt+1)=bkof_amount(numDataPkt+1)+bkof_incre;

                bfRxNode = bfRxNode + 1;
                selectedUser = [selectedUser, userToAdd]
                

                %random beamforming, no feedback               
                if useRBF
                    noCSINum = length(selectedUser);

                    %TODO: replace CSI with its fake direction ID
                    for subc = 7:59
                        if subc==33 continue; end;
                        CSI(userToAdd,:,subc)=iPrecodeM(:,optDirID(userToAdd),subc).';     
                    end

                    fprintf(1, '------ Selected RX %d has good alignment with direction %d, no CSI feedback\n', userToAdd, optDirID(userToAdd)); 
                end
                
            end
        else
            % if we cannot collect M users, then go beamforming to already selected users
            numUserLeft = 0;
        end
    end
end


%% plot for every data packet & Reset
if(pktType == 3)





                %TODO: evaluate the throughput based on capacity and overhead
                ts_duration = 9e-6;
                %maybe we can assume packet aggregation here?
                %when the users is not selected, it can aggregate packets
                %pktDuration = 5e-3;
                pktLength = 1500*8;
                pktDuration = pktLength/20e6;
                numQbits = 8;
                Time_preamble = 20e-6; 
                % the duration of each probing sequence
                probeDuration = 20e-6;

                % overhead of each CSI feedback
                pollingTransRate = 6.5e6; %the transmission rate for polling and CSI is 6Mbps
                %pktCSILength= numTxAntenna * 2 *(numQbits+3)*52; 
                pktCSILength= numTxAntenna *(numQbits+3)*52; 
                pktCSITime= pktCSILength/pollingTransRate;

                dataAmount = sum(thruRecord(numDataPkt,:))*pktDuration;
 
                % the overhead is divided for each of the packets between two user selections
                if useOSel
                    timeAmount = (Time_preamble + (length(selectedUser)-noCSINum)*pktCSITime  ...
                                               + (length(selectedUser)-1)*probeDuration ...
                                               + bkof_amount(numDataPkt)*ts_duration)/pktInterval ...
                                               + pktDuration;
                    if colFlag 
                        %detect collision only after transmitting this CSI
                        timeAmount = timeAmount + pktCSITime + probeDuration;
                        colFlag = 0;
                    end
                else
                    timeAmount = (Time_preamble + length(selectedUser)*pktCSITime)/pktInterval  ...
                                               + pktDuration;
                end

                    sumRate(numDataPkt)=dataAmount/timeAmount;

                    bkof_fraction(numDataPkt) = bkof_amount(numDataPkt)*ts_duration/timeAmount;
                    
                    mDataAmount(selectedUser)= mDataAmount(selectedUser)+thruRecord(numDataPkt,selectedUser)*pktDuration;
                    mTimeAmount = mTimeAmount + timeAmount;



    clear legendInfo
    
       
    cc=hsv(2*streamNum);
    if modeID ==1   %SU-MIMO
        
        figure(4)
        for rxIta = 1:numRxNode
            if rxIta~=rxID
                continue;
            else
                rxn = rxIta;
                legendInfo=['SINR Rx' num2str(rxn)];
            end
            plot(1:64,Subc_SINR(rxn,:),'color',cc(2*rxn,:)) 
            hold on 
        end   
            hold off
            legend(legendInfo)
            xlabel('subc ID')
            ylabel('dB')
            
    elseif modeID == 2 %MU-MIMO
        
        figure(4)
        %plot SINR
        %subplot(1,2,1)
        for rxIta = 1:bfRxNode
            rxn = selectedUser(rxIta);
            legendInfo{rxIta}=['SINR Rx' num2str(rxn)];
            plot(1:64,Subc_SINR(rxn,:),'color',cc(2*rxn,:))        
            hold on 
        end  
            hold off
            legend(legendInfo)
            xlabel('subc ID')
            ylabel('dB')
            
            %{
        %plot condition number
        %subplot(1,2,2)
        figure(40)
            plot(1:64,CNum(:,numDataPkt),'r')     
            %hold on
            %plot(1:64,CNum_w(:,numDataPkt),'b')    
            hold off
            %legend('CN_H','CN_w')
            legend('CN')
            xlabel('subc ID')
            ylabel('dB')
            %}
            
    elseif modeID == 3 %open loop MIMO
        
        figure(4)
        %plot SINR
        %subplot(1,2,1)
        for rxa = 1:numRxAntenna
            legendInfo{rxa}=['SINR Rx antenna' num2str(rxa)];
            plot(1:64,Subc_SINR(rxa,:),'color',cc(2*rxa,:))        
            hold on 
        end  
            hold off
            legend(legendInfo)
            xlabel('subc ID')
            ylabel('dB')
    end


end



%% end: determine whether to request polling 



end % end for numPrecodeRounds = 1:totalDataPkts

%t = cputime-startTime;
%fprintf(1, 'WARP delay: %g\n', t/totalDataPkts/2/1.5);


% ============ Reset and disable the boards ==========
if USESIM == 0
    wl_basebandCmd(nodes,sum(RF_vector),'tx_rx_buff_dis');
    wl_interfaceCmd(nodes,sum(RF_vector),'tx_rx_dis');
end

%% start: performance evaluateion
clear legendInfo;
figure(1)

for rxn = 1:numRxNode
        legendInfo{rxn}=['SINR Rx' num2str(rxn)];
        plot(1:numDataPkt, pktSINR(:,rxn)','color',cc(2*rxn,:))
        hold on 
end

legend(legendInfo)
xlabel('packet ID')
ylabel('SINR')

if showAlignment
figure(20)
for rxn = 1:numRxNode
        legendInfo{rxn}=['hk-wk correlation Rx' num2str(rxn)];
        plot(1:numPrecodeRounds, align_store(rxn,:),'color',cc(2*rxn,:))
        hold on 
end

legend(legendInfo)
xlabel('precode round')
ylabel('alignment with corresponding w')
end

%{
figure(2)
for rxn = 1:numRxNode
        legendInfo{rxn}=['CQI Rx' num2str(rxn)];
        plot(1:numDataPkt, CQI_store(:,rxn)','color',cc(2*rxn,:))
        hold on 
end

legend(legendInfo)
xlabel('packet ID')
ylabel('CQI')

figure(3)
for rxn = 1:numRxNode
        legendInfo{rxn}=['realCQI Rx' num2str(rxn)];
        plot(1:numDataPkt, realCQI_store(:,rxn)','color',cc(2*rxn,:))
        hold on 
end

legend(legendInfo)
xlabel('packet ID')
ylabel('realCQI')




figure(21)
        plot(1:numPrecodeRounds, cond_store,'r')

xlabel('precode round')
ylabel('CN')

%}


%{
figure(23)
        plot(1:numPrecodeRounds, den,'r')

xlabel('precode round')
ylabel('maxamp')
%}
%{
figure(2)
for rxn = 1:bfRxNode
    if rxn == 1
        plot(1:numDataPkt, thruRecord(:,1)','r');
        hold on
    else
        plot(1:numDataPkt, thruRecord(:,2)','b');    
    end
end
legend('Rx1','Rx2')
xlabel('packet ID')
ylabel('Throughput')
%}

%{
figure(3)
for rxn = 1:bfRxNode
    if rxn == 1
        plot(1:numValidDataPkt, validThruRecord(:,1)','r');
        hold on
    else
        if length(pktSINR(rxn,:))>1
            plot(1:numValidDataPkt, validThruRecord(:,2)','b');   
        end
    end
end
legend('Rx1','Rx2')
xlabel('packet ID')
ylabel('Valid Throughput')
%}


    %{
    figure(5)
    for rxn = 1:bfRxNode
        PAplot(1,:) = PA_record(rxn,:,totalValidDataPkts);
        if rxn == 1        
            plot(1:59, PAplot,'r');
            hold on
        else
            plot(1:59, PAplot,'b');    
        end
    end
    legend('Rx1','Rx2')
    xlabel('packet ID')
    ylabel('Power fraction')
    %}
   %{
    figure(5)
    for rxn = 1:bfRxNode
        CHplot_t(1,:) = chGain_record(rxn,20,:);
        if rxn == 1        
            plot(1:length(CHplot_t), CHplot_t,'r');
            hold on
        else
            plot(1:length(CHplot_t), CHplot_t,'b');    
        end
    end
    legend('Rx1','Rx2')
    xlabel('packet ID')
    ylabel('ch Gain')
   
    figure(6)
    for rxn = 1:bfRxNode
        CHplot(1,:) = chGain_record(rxn,:,totalValidDataPkts);
        if rxn == 1        
            plot(1:59, CHplot,'r');
            hold on
        else
            plot(1:59, CHplot,'b');    
        end
    end
    legend('Rx1','Rx2')
    xlabel('subc ID')
    ylabel('ch Power Gain')
   
    
    figure(8)
    for rxn = 1:bfRxNode
        CHplot8(1,:) = abs(avgCHout(rxn,1,:));
        if rxn == 1        
            plot(1:64, CHplot8,'r');
            hold on
        else
            plot(1:64, CHplot8,'b');    
        end
    end
    title('ch Gain txa1')
    legend('Rx1','Rx2')
    xlabel('subc ID')
    ylabel('ch Gain')
     
    
    figure(9)
    for rxn = 1:bfRxNode
        CHplot9(1,:) = abs(avgCHout(rxn,2,:));
        if rxn == 1        
            plot(1:64, CHplot9,'r');
            hold on
        else
            plot(1:64, CHplot9,'b');    
        end
    end
    title('ch Gain txa2')
    legend('Rx1','Rx2')
    xlabel('subc ID')
    ylabel('ch Gain')
    %} 
    %{
    figure(7)
    
              
            plot(1:numDataPkt, CNum(20,:),'r');
        
    
    
    xlabel('num data pkt')
    ylabel('CN (not dB)')
    %}
%{
if userSelMode==3 && modeID>1
    
    figure(10)
    plot(1:numDataPkt, corrVal,'r')
    hold on;
    plot(1:numDataPkt, corrVal2,'g')
    hold on;
    plot(1:numDataPkt, corrVal3,'b')
    legend('abs','real','imag')
    
    figure(12)
        plot(1:numDataPkt, logMagDiff,'r')   
    legend('logMagDiff')
end

    figure(11)
    plot(1:numDataPkt, den,'r')   
    legend('maxamp')

    figure(12)
    plot(1:numDataPkt, normW,'r')   
    legend('normW')

    if modeID>1
        figure(13)
        plot(1:length(CNum_aver), CNum_aver,'r')   
        legend('CN aver')
    end
    
   figure(14)
    plot(1:length(phaseDiff), phaseDiff,'r')   
    legend('phaseDiff')    
    xlabel('t')
    ylabel('phase difference')

    
figure(15)
clear legendInfo

mTmp = size(CNum_cross);
numLines = mTmp(2);
cc=hsv(numLines);
for pairID = 1:numLines
        legendInfo{pairID}=['Group: ' num2str(pairStore(pairID,1)) ' and ' num2str(pairStore(pairID,2))];
        plot(1:mTmp(1), CNum_cross(:, pairID)','color',cc(pairID,:))
        hold on 
end

legend(legendInfo)
xlabel('packet ID')
ylabel('pair CN')

pairStore(:,3) = (sum(CNum_cross)/mTmp(1))'
%}

selectRecord

mCapacity = mean(thruRecord)

sumCapacity = sum(mCapacity)

avgSumRate = mean(sumRate)

mThroughput = mDataAmount/mTimeAmount

JainIndex = sum(mThroughput)^2/(numRxNode*sum(mThroughput.^2))
