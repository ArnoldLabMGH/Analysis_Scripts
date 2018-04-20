Freq=31.25;
N_freq=31.25/2;
Time=1/Freq;
Length=size(subj1_accel,1);
t=linspace(0,Length,Length)*Time;

figure(1)
plot(time,subj1_accel)
grid

FFT_subj1_accel=fft(subj1_accel)/Length;
Fv=linspace(0,1,fix(Length/2)+1)*N_freq;
Iv=1:length(Fv);

figure(2)
semilogx(Fv,abs(FFT_subj1_accel(Iv))*2);
grid

Wp=1/N_freq;
Ws=1.1/N_freq;
Rp=10;
Rs=50;
[n,Ws]=cheb2ord(Wp,Ws,Rp,Rs);
[z,p,k]=cheby2(n,Rs,Ws);
[sosbp,gbp]=zp2sos(z,p,k);

figure(3)
freqz(sosbp,2^16,Freq)

s_filt=filtfilt(sosbp,gbp,subj1_accel);

figure(4)
plot(t,subj1_accel,'-b')
hold on
plot(t,s_filt,'-r','Linewidth',1.5);
hold off
xlabel('Time')
ylabel('Amplitude')
legend('Original','Lowpass Filtered')
