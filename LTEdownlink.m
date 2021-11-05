% Graduation Project_LTE Technology
% Downlink Physical Layer
% edit by LTE Team
% Setup
% Define parameters.
clear all; clc;
nFFT = 256; % fft size
nDSC = nFFT; % number of data subcarriers
CP = 1/4;
CP_start = (nFFT-(nFFT*CP))+1;
CP_end = (nFFT*CP)+1;
Total_Lingth = (nFFT*CP)+nFFT;
M = 64; % Size of signal constellation
k = log2(M); % Number of bits per symbol
n_FRM = 300; 
FRM = 6144;
FRM = FRM-24;
K1 = FRM+24;
codeRate = 1/3;
f1=263; f2=480; inx=0:(K1-1);
Indices = rem((f1*inx)+(f2*(inx.^2)),K1)+1;
Trellis = poly2trellis(4, [13 15], 13);
TurboEnc = comm.TurboEncoder('TrellisStructure',Trellis,'InterleaverIndices',Indices);
TurboDec = comm.TurboDecoder('TrellisStructure',Trellis,'InterleaverIndices',Indices,'NumIterations',6);
zm = comm.RectangularQAMModulator(M, 'SymbolMapping', 'Gray', 'BitInput', true,'NormalizationMethod','Average power','OutputDataType','double');
zd = comm.RectangularQAMDemodulator(M, 'SymbolMapping', 'Gray', 'BitOutput', true,'NormalizationMethod','Average power','DecisionMethod','Log-likelihood ratio','VarianceSource','Input port');
CRC_Tx = comm.CRCGenerator('Polynomial',[1 1 zeros(1, 16) 1 1 0 0 0 1 1]);
CRC_Rx = comm.CRCDetector('Polynomial', [1 1 zeros(1, 16) 1 1 0 0 0 1 1]);
if M==2
EbNo = [0:8];
elseif M==4
EbNo = [0:8];
elseif M==16
EbNo = [0:8];
elseif M==64
EbNo = [0:8];
elseif M==256
EbNo = [0:8];
end
bit_error_rate=[];
resource_block=12;
No_Code_Blocks=21;
for ii = 1:length(EbNo)
    ii
 Total_number_of_errors=0;
 for jj=1:n_FRM
%% Signal Source
t_data = randi([0 1], FRM,1); % Random binary data stream
%% CRC 
t_data_crc = step(CRC_Tx,t_data);
%% Turbo Coding
cod_data=step(TurboEnc,t_data_crc);
%% Interleaving
D = K1 +4;
cod_data_Interleaver = reshape(cod_data,3,D);
colTcSb = 32;
rowTcSb = ceil(D/colTcSb);
Kpi = colTcSb * rowTcSb;
Nd = Kpi - D;
d =(1:D)';
colPermPat =[0,16,8,24,4,20,12,28,2,18,10,26,6,22,14,30,1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31];
y = [NaN*ones(Nd, 1); d]; 
inpMat = reshape(y, colTcSb, rowTcSb).';
permInpMat = inpMat(:, colPermPat+1);
Index1 = permInpMat(:);
% For 3rd stream only
pai = zeros(colTcSb*rowTcSb, 1);
for i = 1 : length(pai)
pai(i) = colPermPat(floor((i-1)/rowTcSb)+1) + colTcSb*(mod(i-1, rowTcSb)) + 1;
end
ytemp = inpMat.';
y2 = ytemp(:);
Index2 = y2(pai);
 for inv=1:3
 if inv==3
 Index=Index2;
 else
 Index=Index1;
 end
 IndexG=find(~isnan(Index)==1);
 IndexB=find(isnan(Index)==1);
 cod_data_Interleaver_temp=cod_data_Interleaver(inv,:);
 cod_data_Interleaver_temp1=zeros(size(Index));
 
cod_data_Interleaver_temp1(IndexG)=cod_data_Interleaver_temp(Index(IndexG));
 Nd=numel(IndexB);
 cod_data_Interleaver_temp1(IndexB)=nan*ones(Nd,1);
 
cod_data_Interleaver_temp2=cod_data_Interleaver_temp1(isfinite(cod_data_Interleaver_temp1));
cod_data_Interleaver_temp3(inv,:)= cod_data_Interleaver_temp2.';
 end
Interleaved_data =reshape(cod_data_Interleaver_temp3,D*3,1); % Buffer 
%% Rate Matching
if codeRate==(1/3); ek_Length=length(Interleaved_data); else
ek_Length=k*ceil((D/k)/codeRate);end
Interleaved_data_matched=Interleaved_data(1:ek_Length);
%% Modulation
mod_data = step(zm, Interleaved_data_matched); 
%% Data Mapping
allocated_resource_block=No_Code_Blocks*resource_block;
no_OFDM_Symbol=ceil((ek_Length/k)/allocated_resource_block);
Dummy_zeros= zeros((allocated_resource_block*no_OFDM_Symbol)-(ek_Length/k),1);
All_Data=[mod_data;Dummy_zeros];
%% serial to parallel conversion
par_data = reshape(All_Data,allocated_resource_block,no_OFDM_Symbol).';
clear mod_data;
%% DFT-precoding.
par_data_freq = fft(par_data,[],2);
%% Subcarrier mapping.
par_data_freq_LFDMA = zeros(no_OFDM_Symbol,nDSC);
subband=0;
par_data_freq_LFDMA(:,[1:allocated_resource_block]+(allocated_resource_block*subband)) = par_data_freq; % Localized mapping.
%% fourier transform time doamain data and normalizing the data
IFFT_data = (nFFT/sqrt(nDSC))*ifft(par_data_freq_LFDMA,nFFT,2);
clear par_data;
%% addition cyclic prefix
cylic_add_data = [IFFT_data(:,[CP_start:nFFT]) IFFT_data];
clear IFFT_data;
%% parallel to serial coversion
ser_dataIN = reshape(cylic_add_data.',(Total_Lingth)*no_OFDM_Symbol,1);
clear cylic_add_data; 
%% AWGN Channel
snr(ii) = EbNo(ii)+(10*log10(k))+ 10*log10(codeRate)-10*log10(nDSC/(allocated_resource_block));
chan_awgn = awgn(ser_dataIN,snr(ii),'measured');
clear ser_dataIN;
%% serial to parallel coversion
para_dataOUT= reshape(chan_awgn,Total_Lingth,no_OFDM_Symbol).'; 
%% removing cyclic prefix
cyclic_pre_rem = para_dataOUT(:,[CP_end:(nFFT+(nFFT*CP))]);
clear para_dataOUT;
%% converting to frequency domain
FFT_recdata= (sqrt(nDSC)/nFFT)*fft(cyclic_pre_rem,nFFT,2); 
clear cyclic_pre_rem;
%% Subcarrier de-mapping.
FFT_recdata_LFDMA = FFT_recdata(:,[1:allocated_resource_block]+(allocated_resource_block*subband));
%% DFT-Despreading
FFT_recdata_time = ifft(FFT_recdata_LFDMA,[],2);
%% parallel to serial coversion
ser_data= reshape(FFT_recdata_time.',allocated_resource_block*no_OFDM_Symbol,1);
clear FFT_recdata;
%% Data de-mapping.
ser_data_demapping = ser_data(1:(ek_Length/k),:);%% Subcarrier demapping.
%% Demodulation
demod_Data = step(zd, ser_data_demapping,10^(-snr(ii)/10)); 
%% Rate Dematching
demod_Data_Dematched=zeros(3*D,1);
demod_Data_Dematched(1:numel(demod_Data))=demod_Data;
%% De-Interleaving
R_data_Denterleaver = reshape(demod_Data_Dematched,3,D); % Buffer
 for dnv=1:3
if dnv==3
 Index=Index2;
 else
 Index=Index1;
 end
 IndexG=find(~isnan(Index)==1); 
 IndexOut=Index(IndexG);
 R_data_Denterleaver_temp=(R_data_Denterleaver(dnv,:));
 
R_data_Denterleaver_temp2=zeros(size(R_data_Denterleaver_temp));
 R_data_Denterleaver_temp2(IndexOut)=R_data_Denterleaver_temp;
 R_data_Denterleaver_temp3(dnv,:)= R_data_Denterleaver_temp2.';
 end
DeInterleaved_data =reshape(R_data_Denterleaver_temp3,D*3,1);
%% Turbo - Decoding
R_data=step(TurboDec,-DeInterleaved_data);
%% CRC
R_data_crc = step(CRC_Rx,R_data);
%% BER Computation
[number_of_errors,bit_error_rate_temp] = biterr(t_data,R_data_crc);
clear t_data ; clear demod_Data ;
Total_number_of_errors=Total_number_of_errors+number_of_errors; 
 end
 bit_error_rate(ii)=Total_number_of_errors/(FRM*n_FRM);
 clear Total_number_of_errors
 bit_error_rate
end
%% ploating
semilogy(EbNo,bit_error_rate,'--sr','linewidth',2);
hold on;
if M==2 % Theoretical BER for 
BPSK
BER_theory = berawgn(EbNo,'psk',M,'nondiff')
else % Theoretical BER for QAM
BER_theory = berawgn(EbNo,'qam',M)
end
semilogy(EbNo,BER_theory,'-b','linewidth',2);
legend('simulated','Theoritical without Coding'); grid on;
xlabel('E_b/N_o'); ylabel('BER')
title('BER Vs. E_b/N_o for M-QAM OFDM in AWGN channel');
axis([0 EbNo(end) 10^-5 1]);
