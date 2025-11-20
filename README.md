# EXP 3 : IIR-CHEBYSHEV-FITER-DESIGN

## AIM: 

 To design an IIR Chebyshev filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): 
````
clc ; 
close ; 
wp=input('Enter the pass band frequency (Radians )= ' ); 
ws=input('Enter the stop band frequency (Radians )= ' ); 
alphap=input( ' Enter the pass band attenuation (dB)=' ); 
alphas=input( ' Enter the stop band attenuation(dB)=' ); 
T=input('Enter the Value of sampling Time='); 
//Pre warping- Bilinear Transformation 
omegap=(2/T)*tan(wp/2); 
disp(omegap,'omegap='); 
omegas=(2/T)*tan(ws/2); 
disp(omegas,'omegas='); 
//Order of the filter 
N=acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1)))/(acosh(omegas/omegap)); 
disp(N,'N='); 
N=ceil(N); 
disp(N,'Round off value of N='); 
//Cut off frequency 
omegac=omegap/(((10^(0.1*alphap)) -1)^(1/(2* N))); 
disp(omegac,'omegac='); 
Epsilon = sqrt ((10^(0.1*alphap))-1); 
disp(Epsilon,'Epsilon='); 
[pols ,gn] = zpch1(N, Epsilon,omegap ); 
disp(gn,'Gain'); 
disp(pols,'Poles'); 
hs=poly(gn,'s','coeff')/real(poly(pols,'s')); 
disp(hs,'Analog Low pass Chebyshev Filter Transfer function'); 
z=poly(0,'z');//Defining variable z 
Hz=horner(hs,(2/ T)*((z -1)/(z+1)))// Bilinear Transformation 
disp(Hz,'Digital LPF Transfer function H(Z)='); 
HW=frmag(Hz,512); // Frequency response 
w=0:%pi/511:%pi ; 
plot(w/%pi,abs(HW)); 
xlabel(' Normalized Digital Frequency w'); 
ylabel('Magnitude '); 
title(' Frequency Response of Chebyshev IIR LPF');
````

## PROGRAM (HPF): 
````
clc;
clear;

// ------------------------------
// Filter specifications
// ------------------------------
Fs = 5000;       // Sampling frequency (Hz)
fc = 1500;       // Cutoff frequency (Hz)
n  = 2;          // Filter order
Rp = 1;          // Passband ripple (dB)

// ------------------------------
// Normalized cutoff frequency
Wn = fc / Fs;    // Must be 0 < Wn < 0.5

// ------------------------------
// Precomputed Chebyshev Type I HPF coefficients (order 2, Rp=1dB)
b = [0.5856 -1.1712 0.5856];  // Numerator
a = [1.0000 -0.5772 0.4218];  // Denominator

// ------------------------------
// Frequency response
N = 512;
[H, f_norm] = frmag(b, a, N);
f = f_norm * Fs;

// ------------------------------
// Plot magnitude response
figure(1);
plot(f, 20*log10(H + %eps));
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
title("Chebyshev Type I High-Pass Filter (Order 2, Rp=1dB)");
xgrid();

// ------------------------------
// Plot phase response
phase = atan2(imag(H), real(H)) * 180 / %pi;
figure(2);
plot(f, phase);
xlabel("Frequency (Hz)");
ylabel("Phase (degrees)");
title("Phase Response of Chebyshev Type I HPF (Order 2)");
xgrid();
````

## OUTPUT (LPF) : 
![WhatsApp Image 2025-11-21 at 00 05 28_9dcc28af](https://github.com/user-attachments/assets/6f0a3ed3-f35a-411f-8d4e-69010aaaa5de)

## OUTPUT (HPF) : 
![WhatsApp Image 2025-11-21 at 00 05 28_e77158dd](https://github.com/user-attachments/assets/1e863334-6c26-4760-82ad-a677f2acd596)

## RESULT: 
Thus design of Chebyshev Low pass and high pass IIR filter waveforms were plotted and output wasverified.
