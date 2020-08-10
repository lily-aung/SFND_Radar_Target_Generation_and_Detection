clear
clc
close all

%% Radar Specs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
% speed of light = 3e8
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity remains contant
d0 = 80;
v0 = -50;

%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
rangeResolution = 1;
c = 3*10^8;
maxRange = 200;
rttScale = 5.5;

Bsweep = c/(2*rangeResolution);  %Bandwidth
Tchirp = rttScale*2*maxRange/c;  %chirp time
slope = Bsweep/Tchirp;           %slope of chirps
fc= 77e9;                        %carrier frequency of Radar 
                                                       
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

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

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
%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.

%run the FFT on the beat signal along the range bins dimension (Nr) andnormalize.
sig_fft = fft(Mix,Nr)./Nr;

 % *%TODO* :
sig_fft = abs(sig_fft);       % Take the absolute value of FFT output

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sig_fft = sig_fft(1:(Nr/2));            % half of the samples

figure ('Name','Range from First FFT')  %plotting the range
plot(sig_fft);                          % plot FFT output 
axis ([0 200 0 0.5]);
xlabel('measured range');

%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

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

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;   % Training (range dimension)
Td = 4     % Training cells (doppler dimension)

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 5;    % Guard cells (range dimension)
Gd = 2;    % Guard cells (doppler dimension)

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 1.4;

% *%TODO* :
% Calculate the total number of training and guard cells
gridSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
trainingCellsNum = gridSize-(2*Gr+1)*(2*Gd+1);
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(Nr/2-2*(Td+Gd),Nd-2*(Tr+Gr));

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

CFAR = zeros(size(RDM));

   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR
for range_index=1:Nd-2*(Tr+Gr)
    for doppler_index=1:Nr/2-2*(Td+Gd)
        %I want to extract only the Training cells. To do that, I at first get
        %the sliding patch and, after converting it to power from decibel,
        %I set to zero whatever it is not in the position of
        %training cells; the zero submatrix will not contribute to the sum
        %over the whole patch and, by dividing for the number of training 
        %cells I will get the noise mean level
        trainingCellsPatch = db2pow(RDM(doppler_index:doppler_index+2*(Td+Gd),range_index:range_index+2*(Gr+Tr)));
        trainingCellsPatch(Td+1:end-Td,Tr+1:end-Tr) = 0;
        noise_level(doppler_index,range_index) = pow2db(sum(sum(trainingCellsPatch))/trainingCellsNum);
        % Use the offset to determine the SNR threshold
        threshold = noise_level(doppler_index,range_index)*offset;
        % Apply the threshold to the CUT
        if RDM(doppler_index+(Td+Gd),range_index+(Td+Gr))>threshold
            CFAR(doppler_index+(Td+Gd),range_index+(Td+Gr)) = 1;
        else
            CFAR(doppler_index+(Td+Gd),range_index+(Td+Gr)) = 0;
        end
           
    end
end


% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 

%not necessary

% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','2D CA-CFAR Filtered RDM');
surf(doppler_axis,range_axis,CFAR);
colorbar;
