# SFND Radar Target Generation and Detection
Report  - Writeup

This project uses Matlab to introduce frequency modulated continuous-wave (FMCW) radar and related post-processing techniques. The topics covered include:
1) Fast Fourier transforms (FFT) and 2D FFT
2) Clutter v. target discrimination
3) Sizing chirp bandwith to meet system requirements for range resolution
4) Phased array beam steering to determine angle of arrival (AoA)
5) Constant false alarm rate (CFAR) noise suppression
6) Signal-to-noise ratio (SNR) and dynamic thresholding

---
#### CRITERIA 1. FMCW Waveform Design 
```Matlab
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc= 77e9;             %carrier freq
rangeResolution = 1;
c = 3e8;
maxRange = 200;
rttScale = 5.5;
Bsweep = c/(2*rangeResolution);  %Bandwidth
Tchirp = rttScale*2*maxRange/c;  %chirp time
slope = Bsweep/Tchirp;           %slope of chirps
```

#### CRITERIA 2. Simulation Loop

```Matlab
d0 = 80;         %initial position 
v0 = -50;        %initial velocity

Nd=128;                %# of doppler cells OR # of sent periods % number of chirps
Nr=1024;               % for length of time OR # of range cells      

t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
L = length(t);
%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,L);    %transmitted signal
Rx=zeros(1,L);    %received signal
Mix = zeros(1,L); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,L);
td=zeros(1,L);
for i=1:L         
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = d0 + v0*t(i);
    td(i) = 2*r_t(i)/c;
    
    %For each time sample we need update the transmitted and received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+slope*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + (slope*(t(i)-td(i))^2)/2));
    
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i).*Rx(i);  %Beat Signal
    
end
```

#### CRITERIA 3. Range FFT (1st FFT)

```Matlab
sig_fft = fft(Mix,Nr)./Nr;

sig_fft = abs(sig_fft);       % Take the absolute value of FFT output
sig_fft = sig_fft(1:(Nr/2));  % only half of the samples

figure ('Name','Range from First FFT')  %plotting the range
plot(sig_fft);                          % plot FFT output 
axis ([0 200 0 0.5]);
xlabel('measured range');

```
[[/results/Figure1-Range from First FFT.png|Range from First FFT]]


* doppler FFT (2st FFT)

```Matlab
Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','2D FFT Range Doppler Map');
surf(doppler_axis,range_axis,RDM);
```
[[/results/Figure2-2D FFT Range Doppler Map.png|2D FFT Range Doppler Map]]

#### CRITERIA 4. 2D CFAR Implementation

Determine the number of Training cells for each dimension. Similarly, pick the number of guard cells.

```Matlab
Tr = 10;   % Training (range dimension)
Td = 4     % Training cells (doppler dimension)
Gr = 5;    % Guard cells (range dimension)
Gd = 2;    % Guard cells (doppler dimension)
      
offset = 1.4;  % offset the threshold by SNR value in dB
```
Calculate the total number of training and guard cells
```Matlab
gridSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
trainingCellsNum = gridSize-(2*Gr+1)*(2*Gd+1);
```
Create a vector to store noise_level for each iteration on training cells.
```Matlab
noise_level = zeros(Nr/2-2*(Td+Gd),Nd-2*(Tr+Gr));
%
```

Slide the cell under test across the complete matrix. 
Make sure the CUT has margin for Training and Guard cells from the edges.

```Matlab
CFAR = zeros(size(RDM));
for range_index=1:Nd-2*(Tr+Gr)
    for doppler_index=1:Nr/2-2*(Td+Gd)
        ...
    end
end
```

For every iteration,  convert the value from logarithmic to linear using db2pow function, then sum the signal level within all the training cells, and get the mean value. then convert it back to db using pow2db.

```Matlab
        trainingCellsPatch = db2pow(RDM(doppler_index:doppler_index+2*(Td+Gd),range_index:range_index+2*(Gr+Tr)));
        trainingCellsPatch(Td+1:end-Td,Tr+1:end-Tr) = 0;
        noise_level(doppler_index,range_index) = pow2db(sum(sum(trainingCellsPatch))/trainingCellsNum);
```

Add the offset to it to determine the SNR threshold.

```Matlab
        threshold = noise_level(doppler_index,range_index)*offset;
```

Apply the threshold to the CUT

```Matlab
        if RDM(doppler_index+(Td+Gd),range_index+(Td+Gr))>threshold
            CFAR(doppler_index+(Td+Gd),range_index+(Td+Gr)) = 1;
        else
            CFAR(doppler_index+(Td+Gd),range_index+(Td+Gr)) = 0;
        end
```


[[/results/Figure3-2D CA-CFAR Filtered RDM.png|2D CA-CFAR Filtered RDM]]
