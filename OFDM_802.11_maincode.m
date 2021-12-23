close all;
packet=randi([0 1],1,4160);    %Random sequence generation

% BPSK Modulation
bpsk_packet=2*packet-1;        
bpsk_packet=[bpsk_packet -1*ones(1,16)]; %padding zeros to make it divisible by 48
group=zeros(87,48);

%grouping 
j=0;

%making a 2D matrix where each row has 48 bits 
for i=1:48:size(bpsk_packet)-47
    j=j+1;
    group(j,:)=bpsk_packet(i:(i+47));   
end

X=zeros(87,64); %ifft input OFDM

%creating the input for IFFT block, inserting pilots and nulls accordingly 
for i=1:1:87
     
     %nulls
     X(i,1)=0; 
     X(i,28:38)=0; 
     
     %pilots
     X(i,8)=1;
     X(i,22)=1;
     X(i,58)=1;
     X(i,44)=1;

     %bpsk bits
     X(i,2:7)=group(i,25:30);
     X(i,9:21)=group(i,31:43);
     X(i,23:27)=group(i,44:48);
     X(i,39:43)=group(i,1:5);
     X(i,45:57)=group(i,6:18);
     X(i,59:64)=group(i,19:24);
end 

x=zeros(87,64);

%calculating the ifft and storing it in 2D array x
for i=1:1:87
    temp=X(i,:);
    x(i,:)=ifft(temp,64);
end 

%add cyclic prefix to each OFDM symbol
x_signal_with_cp=zeros(87,80);

for i=1:1:87
    x_signal_with_cp(i,1:16)=x(i,49:64);
    x_signal_with_cp(i,17:80)=x(1,1:64);
end 

%Preamble: Addition of STF and LTF

stf={ 0+0*1i, 0+0*1i, 1+ 1i, 0+0*1i, 0+0*1i, 0+0*1i, -1-1i, 0+0*1i, 0+0*1i, 0+0*1i, 1+1i, 0+0*1i, 0+0*1i, 0+0*1i, -1-1i, 0+0*1i, 0+0*1i, 0+0*1i, -1-1i, 0+0*1i, 0+0*1i, 0+0*1i, 1+1i, 0+0*1i, 0+0*1i,0+0*1i, 0+0*1i, 0+0*1i, 0+0*1i, 0+0*1i, -1-1i, 0+0*1i, 0+0*1i, 0+0*1i, -1-1i, 0+0*1i, 0+0*1i,0+0*1i, 1+1i, 0+0*1i, 0+0*1i, 0+0*1i, 1+1i, 0+0*1i, 0+0*1i, 0+0*1i, 1+1i, 0+0*1i, 0+0*1i, 0+0*1i, 1+1i, 0+0*1i,0+0*1i};
stf=sqrt(13/6)*cell2mat(stf);

X=zeros(1,64);

%mapping the stf bits to the IFFT input 
X(1)=stf(27);
X(2:27)=stf(28:53);
X(28:38)=0;%null
X(39:64)=stf(1:26);

time_stf=ifft(X,64); %finding time domain representation of STF

stf_periodical=[time_stf time_stf time_stf(1:32)]; %repeated sequence of stf = 8 us

%windowing (according to the standard)
stf_periodical(1)=stf_periodical(1)/2;
stf_periodical(160)=stf_periodical(160)/2;

figure (1)
stem(abs(stf_periodical))
title('1 (c) Magnitude of STF (time domain)')
xlabel('bit index')
ylabel('Magnitude (STF)')

%long-training sequence
ltf={1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1};
ltf=cell2mat(ltf);
 
X_ltf=zeros(1,64);

%mapping LTF frequency domain values for giving input to IFFT block
X_ltf(1)=ltf(27);
X_ltf(2:27)=ltf(28:53);
X_ltf(28:38)=0;%null
X_ltf(39:64)=ltf(1:26);

time_ltf=ifft(X_ltf,64); %taking IFFT of the LTF

%repeating LTF
ltf_periodical=[time_ltf(33:64) time_ltf time_ltf];

%windowing LTF (according to the standard)
ltf_periodical(1)=ltf_periodical(1)/2;
ltf_periodical(160)=ltf_periodical(160)/2;

%preamble
preamble=[stf_periodical ltf_periodical];

%opening an OFDM 2D array => 1D
i=1;
OFDM_signal=zeros(1,6960);
for k=1:80:6881
        OFDM_signal(k:k+79)=x_signal_with_cp(i,:);
        i=i+1;
end 

%ADDING PREAMBLE TO THE ENTIRE OFDM PACKET
OFDM_signal_with_preamble=[preamble OFDM_signal];

figure (2)
W=fft(OFDM_signal_with_preamble,64);
stem(abs(W.^2))
title('1(c) PSD of OFDM')
xlabel('frequency')
ylabel('Power')

%Adding zeros for idle time
 OFDM_signal_with_idle=[zeros(1,100) OFDM_signal_with_preamble]; %for idle time
  
 %magnitude distortion
 OFDM_magnitude_distorted=OFDM_signal_with_idle*10^(-5);

%phase shift 
 OFDM_phase_shift=zeros(1,size(OFDM_magnitude_distorted,2));
 for i=1:1:size(OFDM_magnitude_distorted,2)
     OFDM_phase_shift(i)=OFDM_magnitude_distorted(i)*exp(-3*1i*pi/4);
 end

 %frequency shift  
OFDM_freq_shift=zeros(1,size(OFDM_magnitude_distorted,2));

 for k=1:1:size(OFDM_magnitude_distorted,2)
     OFDM_freq_shift(k)=OFDM_phase_shift(k)*exp(-1i*2*pi*k*0.00017); %******************DOUBT***************
end 


OFDM_noisy=zeros(1,size(OFDM_magnitude_distorted,2));
%add channel noise
 for i=1:1:size(OFDM_noisy,2)
     OFDM_noisy(i)=OFDM_freq_shift(i)+normrnd(0,sqrt(10^(-11)));
 end 
% 

figure (3)
% (2) PLOTTING THE noisy STF
stf_plot=OFDM_noisy(101:261);
stem(abs(stf_plot))
title('noisy STF plot (time-domain)')
ylabel('Magnitude(STF)')
xlabel('sample index')


% %self-correlating algorithm to detect the STF
R=zeros(1,size(OFDM_noisy,2));
E=zeros(1,size(OFDM_noisy,2));
auto_corr_signal=0;
energy = 0;

%self-correlation function
for m=17:1:size(OFDM_noisy,2)-16
      for i=m:1:(m+16-1)
           auto_corr_signal=auto_corr_signal+OFDM_noisy(i)*conj(OFDM_noisy(i-16));
           energy = energy+abs(OFDM_noisy(i)*(OFDM_noisy(i)));
      end 
      R(m)=auto_corr_signal;
      E(m)=energy;
      energy = 0;
      auto_corr_signal=0;
end 
%plotting autocorrelation results
figure(4)
plot(abs(R))
title('(3) Autocorrelation function for Packet detection')
xlabel('sample index m')
ylabel('R(m)')



R_diff=diff(abs(R));
E_diff=diff(E);

for m=17:1:size(OFDM_noisy,2)-1
    if((R_diff(m)-E_diff(m))<1*10^(-18) && abs(R(m))>2*10^(-11))
        break;
    end 
end 

%packet sync ---------> Cross Correlation
cross_corr=zeros(1,size(OFDM_noisy,2));

stf_correlation_segment = stf_periodical(17:32);
cross_corr_signal = 0;
for i=1:1:size(OFDM_noisy,2)-15
    ofdm_correl_segment = OFDM_noisy(i:i+15);
    for j=1:1:16
        cross_corr_signal = cross_corr_signal + ofdm_correl_segment(j)*conj(stf_correlation_segment(j));
    end
    cross_corr(i) = cross_corr_signal;
    cross_corr_signal = 0;
end 

%plotting crosscorrelation results
figure(5)
plot(abs(cross_corr))
title('(4) Crosscorrelation function for Packet sync')
xlabel('sample index m')
ylabel('Cross_Corr')
% 
[peak_value,peak_index]=maxk(abs(cross_corr),10);
peak_index=sort(peak_index);
STF_start = peak_index(1);
STF_end_minus16 = peak_index(10);
LTF_start = STF_end_minus16 + 48;   %without LTF CP
LTF_end = LTF_start + 127;

% %****************Channel estimation and packet decoding ****************
% 
% %frequency-drift correction
count = 0;
del_freq=zeros(1,64);
gain_sum = 0;

for i=1:1:64
    count = count+1;
    h_channel=OFDM_noisy(LTF_start+(i-1)+64)/OFDM_noisy(LTF_start+(i-1));
    %imag_part = imag(h_channel);
    del_freq(i)=imag(h_channel)/(2*pi*64);
    gain_sum = gain_sum+del_freq(i);
end

delta_f_avg = gain_sum/count;
OFDM_freq_corrected=zeros(1,7380);

for i=1:1:7380
    OFDM_freq_corrected(i) = OFDM_noisy(i)*(exp(-1*1j*2*pi*i*delta_f_avg));
end

% Phase and Magnitude Correction : divide freq corrected signal with known
% LTF for calculating channel gain
OFDM_noisy_LTF_FFT = fft(OFDM_freq_corrected(LTF_start:LTF_start+63),64);

channel_gain=zeros(1,64);

for i=1:1:64
    channel_gain(i) = OFDM_noisy_LTF_FFT(i)/X_ltf(i);
end
figure(6)
stem(abs(channel_gain))
title('(6) Channel distortion')

% %************DECODING********************
OFDM_per_symbol=zeros(87,64);
OFDM_per_symbol_angle=zeros(87,64);

p=1;
for i=421:80:5540
        OFDM_per_symbol(p,:)=fft(OFDM_freq_corrected(i+16:i+16+63),64)/channel_gain(3);
        OFDM_per_symbol_angle(p,:) = angle(OFDM_per_symbol(p,:));
        p=p+1;
end 

OFDM_final_bitmap_with_pilot_null = zeros(87,64);
for i=1:1:87
    for j = 1:1:64
        if(OFDM_per_symbol_angle(i,j)>0)
            OFDM_final_bitmap_with_pilot_null(i,j)=-1;
        else
            OFDM_final_bitmap_with_pilot_null(i,j)=1;
        end
    end
end
% 
OFDM_final_bitmap=zeros(87,48);

for i = 1:1:87
    OFDM_final_bitmap(i,:) = [OFDM_final_bitmap_with_pilot_null(i,2:7) OFDM_final_bitmap_with_pilot_null(i,9:21) OFDM_final_bitmap_with_pilot_null(i,23:27) OFDM_final_bitmap_with_pilot_null(i,39:43) OFDM_final_bitmap_with_pilot_null(i,45:57) OFDM_final_bitmap_with_pilot_null(i,59:64)];
end

%de-map
OFDM_baseband=zeros(87,48);
for i=1:1:87
    OFDM_baseband(i,:)= [group(i,25:30) group(i,31:43) group(i,44:48) group(i,1:5) group(i,6:18) group(i,19:24)];
end

%computing the number of bit-errors after decoding 
error=0;
for i=1:1:87
    for j=1:1:48
        if (OFDM_baseband(i,j)~=OFDM_final_bitmap(i,j))
            error=error+1;
        end 
    end 
end 
