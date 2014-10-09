% Implementation of IEEE 802.15.4 PHY layer
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments: This file contains the 802.15.4 PHY layer implementation. This
% file contains the transmitter side encoding and modulation function. For
% receiver side decoding and demodulation function please check ZigbeeRx.m
% file
%
% Version: 0.0.2
% Last modified: 01/14/2014
% 
% Comments: Implementing a WARP board Tx-Rx communication. Generalizing
% from 1 packet to multiple packets. Working for multiple packets now.
% Plotting bit error positions and BER.
% 
% Version: 0.0.3
% Last modified: 01/17/2014
% 
% Comments: Debugging WARP board transmission problems. Adding carrier
% sensing to remove extraneous signals at the beginning which might cause
% trouble for FM demodulator, DC correction and/or M&M circuit.
% Implementing a frequency offset compensation mechanism at the receiver
% side. Using LTF style preamle at the transmitter side which receiver used
% to compensate for the frequency offset. Debugged the WARP transmission
% problem. Removed the LTF frequency offset compensation. Also, removed the
% parameters for clock recovery. Added the condition to discard packets
% having much higher BER. Discarded packets are accounted for packet loss
% rate.
% WARPLab transmission is now working. Tested with multiple packets.
% Adding the code to reset and disable boards after Tx/Rx is complete.
% 
% 



clc;
clear all;

global CHIP_ENCODE_MAPPING CHIP_DECODE_MAPPING DEFAULT_ERROR_THRESHOLD...
    USESIM osamp PACKET_SYNC_VECTOR numBitsTotal numSymbolsTotal...
    MM_FILTER_NSTEPS MM_FILTER_NTAPS FILTER_TAPS SIMSNR numTxAntenna... 
    numRxAntenna numRxNode numTxNode useFakeChannel DEBUG_ON...
    VERBOSE1 VERBOSE2 numZeroPreamble AllowBERperPacket...
    CARRIER_SENSE_WIN_SIZE RSSthresh padding_size STF_up2 LTF_up2...
    SM_CHAN useSISO;


%% Basic mode control
USESIM = 0; % whether to use simulation
useSISO = 1; % use SISO transmission or Spatial Modulation
% Show constellation or not
DEBUG_ON = 1;
% First level of print outs
VERBOSE1 = 1;
% Second level and detailed print outs
VERBOSE2 = 0;
totalDataPkts = 50; % total data packets in experiment

useFakeChannel = 0; % Use a fake channel to Debug


%% WARP Board PHY parameter configuration
TxDelay = 1000;
TxLength = 2^14-TxDelay; % Length of transmission. In [0:2^14-1-TxDelay]
CarrierChannel = 13; % Channel in the 2.4 GHz band. In [1:14] (avoid...
                     % 1 to 11); 5GHz in [15:37]
TxGain_RF = 35; % Tx RF Gain. In [0:63] 
TxGain_BB = 1; % Tx Baseband Gain. In [0:3]
RxGain_BB = 10; % Rx Baseband Gain. In [0:31]
RxGain_RF = 1; %2; % Rx RF Gain. In [1:3]


%% ====== implicit controls

numTxNode = 1;
numRxNode = 1;
if useSISO
    numTxAntenna = 1; % number of tx antennas      
    numRxAntenna = 1; % number of rx antennas
else % useSISO = 0, Use spatial modulation mode
    % At least 2 transmit antenna
    numTxAntenna = 2; % number of tx antennas
    numRxAntenna = 1; % number of rx antennas
end

if USESIM 
    osamp = 1;
else
    % Upsample by 2 for WARP board transmission
    osamp = 1;
end


%% Assigning values to global variables

numZeroPreamble = 16;

if USESIM
    SIMSNR = 10;
end

% Here each symbol is of 4 bits
numSymbolsTotal = 120;

% Always take a chunk of 4bit at a time
numBitsTotal = numSymbolsTotal*4;

% Packet synchronization vector, used in SHR
PACKET_SYNC_VECTOR = hex2dec('A7');

% Allowable number of error bits in a 32-bit chip sequence
% For simulation with Fake channel there should not be any error, however 
% for Gaussian channel or WARP board communication, allow upto 5 chips 
% and 10 chips error respectively
if USESIM
    if useFakeChannel
        DEFAULT_ERROR_THRESHOLD = 0;
    else
        DEFAULT_ERROR_THRESHOLD = 5;
    end
else
    DEFAULT_ERROR_THRESHOLD = 10;
end

% Allowable BER per packet, default 1e-3
AllowBERperPacket = 1e-3;

% Carrier sensing window size
CARRIER_SENSE_WIN_SIZE = 8;

% Choose proper RSSThresh to help carrier sense the input signal
RSSthresh = 0.05;

% Chips encoding sequence, every 4 bit is mapped to 32 bit
% For Chip mapping, please check,
% http://nesl.ee.ucla.edu/fw/thomas/thomas_project_report.pdf

CHIP_ENCODE_MAPPING = ...
                [3653456430
                3986437410
                786023250
                585997365
                1378802115
                891481500
                3276943065
                2620728045
                2358642555
                3100205175
                2072811015
                2008598880
                125537430
                1618458825
                2517072780
                3378542520];
    
% ====== Read global variables 
ZigbeeglobalVars;   

            
%% Implicit controls
numDataPkt = 0;
            

%% Construct packets & samples
fg = fopen('databits.dat', 'w');
for k = 1:numBitsTotal             
	if (rand < 0.5)
        fprintf(fg, '1 ');
    else
        fprintf(fg, '0 ');
    end
    fprintf(fg, '\n');
end
fclose(fg);

% Input data bits
inb = load('databits.dat');
% Take 4 bit chunk at a time
ins = reshape(inb, 4, numSymbolsTotal).';
% Convert to symbol
ins = bi2de(ins, 'left-msb');


%% ========== Initialization of WARP radio parameters ========
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

if useSISO
    RF_vector = [RFA];
else
    RF_vector = [RFA,RFB];
end

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


% Add 8 contiguous '0' symbols followed by Start Frame Delimiter i.e.
% packet sync vector
sync_vector = [zeros(1,numZeroPreamble) hex2dec('A') 7];

% Concatenate sync vector with data
ins = [sync_vector'; ins];
% To test with only 2 '0' symbols at the beginning
% ins = [0; 0; 0; 0; 0; 0; 0; 0];

% Find the value in the chip encoding table indexed by the symbol
inc = encode_symbols(ins);

input_symbols = [];
for ns = 1:length(inc)
%     fprintf('%x\n', inc(ns));
    % Convert 32-bit integer into bits
    temp = de2bi(inc(ns), 'left-msb');
    % Pad 0 in front if length less than 32 bits
    temp = [zeros(1 , 32-length(temp)) temp];
    % Get every 2 bits
    temp = reshape(temp, 2, 16).';
    % Convert every 2 bits to decimal value, this will be input to qpsk 
    % modulator
    input_symbols = [input_symbols bi2de(temp,'left-msb').'];
end

input_symbols;

% O-QPSK Modulation
% Convert input samples with QPSK half-sine pulse shape symbols
psk_mod_sym = [];
for k = 1:length(input_symbols)
    psk_mod_sym = [psk_mod_sym qpsk_half_sine_mod(input_symbols(k))];
end

psk_mod_sym = psk_mod_sym.';

% Add 2 samples delay in the Q-phase
delay = 2;
real_psk_mod_sym = real(psk_mod_sym);
imag_psk_mod_sym = imag(psk_mod_sym);

imag_psk_mod_sym = [zeros(1, delay).'; imag_psk_mod_sym];
real_psk_mod_sym = [real_psk_mod_sym; zeros(1, delay).'];
psk_mod_sym_delay = real_psk_mod_sym + sqrt(-1)*imag_psk_mod_sym;

if USESIM
    TxData_board = psk_mod_sym_delay;
else
    % Interpolate samples by osamp time
    TxData_board = interp(psk_mod_sym_delay, osamp);
end


if useSISO == 0 % useSISO = 0, Spatial Modulation
    % Add Channel estimator for Spatial Modulation scheme
    TxData_board = [SM_CHAN.'; TxData_board];
end


% Add some padding in front of the samles, FM demodulator and 
% Mueller and Müller (M&M) clock recovery component needs extra 8 taps at 
% the beginning and at the end
if USESIM
    padding_size = 32;
    padding = zeros(1, padding_size);
    TxData_board = [padding.'; TxData_board; padding.'];
    
else % USESIM = 0, pad for WARP board
    padding_size = 128;
    padding = zeros(1, padding_size);
    TxData_board = [padding.'; TxData_board; padding.'];
end

if DEBUG_ON
    figure(201);
    plot(1:length(TxData_board), real(TxData_board), 'b-', ...
         1:length(TxData_board), imag(TxData_board), 'r-');
    xlabel('Sample count');
    ylabel('Tx samples');
    title('Transmitted samples');
end


%% ================= Start transmission rounds =====================

accBER = 0.0;
BER = 0.0;
accCER = 0.0;
CER = 0.0;
correctDataPacket = 0;
header_lost = 0;

while numDataPkt < totalDataPkts
    
fprintf(1, '==== Transmission start for %d-th DATA packet ====\n',...
                                                    numDataPkt+1);
if USESIM % Simulate
        
    if useFakeChannel
        % Consider no channel, just a line transmission :)
        % This is very helpful for debugging purpose
        rxData_board = TxData_board;
        
	else % No fake channel, add additive white Gaussian noise
        
        % AWGN Channel to produce variable SNR condition at the 
        % receiver side, SIMSNR variable controls the parameter. 
        % Signal power is assumed to be of 0 W
        rxData_board = awgn(TxData_board, SIMSNR);
    end

else % USESIM=0, Run on WARP
    
%% ============= Start tx and rx =======================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit and receive signal using WARPLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign the data to transmit
wl_basebandCmd(node_tx(1), RF_vector, 'write_IQ', TxData_board);

% start transmission and reception
wl_interfaceCmd(node_tx, sum(RF_vector),'tx_en');
wl_interfaceCmd(node_rx, sum(RF_vector),'rx_en');

wl_basebandCmd(node_tx, sum(RF_vector),'tx_buff_en');
wl_basebandCmd(node_rx, sum(RF_vector),'rx_buff_en');

eth_trig.send();

% get the received data
rxData_board = wl_basebandCmd(node_rx(1), RF_vector, 'read_IQ', 0,...
                                                        TxLength+TxDelay);
                                                    
end % end if USESIM==1
    
if DEBUG_ON
    figure(202);
    plot(1:length(rxData_board), real(rxData_board), 'b-', ...
         1:length(rxData_board), imag(rxData_board), 'r-');
    xlabel('Sample count');
    ylabel('Rx samples');
    title('Received samples');
end

[BER, CER, syncFound] = ZigbeeRx(rxData_board);

% Check whether the packet has correct header, otherwise discard the BER
% calculation
if syncFound
    if BER < 0.05 % Consider packets which has BER less than 5%, otherwise,
        % treat as a packet loss scenario
        accBER = accBER + BER;
        accCER = accCER + CER;
        correctDataPacket = correctDataPacket + 1;
    end
else
    header_lost = header_lost + 1;    
end

numDataPkt = numDataPkt + 1;

end

%% ============ Reset and disable the boards ==========
if USESIM == 0
    wl_basebandCmd(nodes, sum(RF_vector), 'tx_rx_buff_dis');
    wl_interfaceCmd(nodes, sum(RF_vector), 'tx_rx_dis');
end


%% Performance metrics

fprintf ('Average BER = %0.6g%%, CER = %0.6g%%\n', ...
           (accBER/correctDataPacket)*100, (accCER/correctDataPacket)*100);
% fprintf('Packet loss rate = %0.6g%%\n', ...
%     (numDataPkt-correctDataPacket)/numDataPkt);
% fprintf('Packet header loss rate = %0.6g%%\n', header_lost/numDataPkt);
% fprintf('Packet loss due to BER = %0.6g%%\n', ...
%     (numDataPkt-correctDataPacket-header_lost)/numDataPkt);

fprintf('\n\n');                                           