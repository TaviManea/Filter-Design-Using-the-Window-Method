close all;
alfa=39.2973;
beta=4.5730;
L=7;
M=23;
r=86.4324;

omega_c = 1.2346;
freq_c=omega_c/pi;
K=5000;

w1=boxcar(M);

w2=triang(M);

w3=blackman(M);

w4=chebwin(M,r);
w4_1=chebwin(M,r-5);
w4_2=chebwin(M,r+5);

w5=hamming(M);

w6=hanning(M);

w7=kaiser(M,beta);
w7_1=kaiser(M,beta-1);
w7_2=kaiser(M,beta+1);

w8=tukeywin(M,alfa/100);
w8_1=tukeywin(M,24.2973/100);
w8_2=tukeywin(M,54.2973/100);

w9=zeros(1,M);
w9_1=zeros(1,M);
w9_2=zeros(1,M);

for n=0:M-1
w9(n+1)= (sin(2*pi*((2*(n)-M+1) / (2*(M-1))))/(2*pi *((2*(n)-M+1) / (2*(M-1))))).^L;
w9_1(n+1)=(sin(2*pi*((2*(n)-M+1) / (2*(M-1))))/(2*pi *((2*(n)-M+1) / (2*(M-1))))).^(L-1);
w9_2(n+1)=(sin(2*pi*((2*(n)-M+1) / (2*(M-1))))/(2*pi *((2*(n)-M+1) / (2*(M-1))))).^(L+1);
end

 w9=w9';
 w9(12)=1;


w9_1=(w9_1)';
w9_1(12)=1;


w9_2=w9_2';
w9_2(12)=1;


figure;
sgtitle('Ferestre fara parametri');
subplot(2,2,1)
stem(w2);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra triunghiulara');
subplot(2,2,2)
stem(w3);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra blackman');
subplot(2,2,3)
stem(w5);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra hamming');
subplot(2,2,4)
stem(w6);
title('Fereastra hanning');   
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
   

figure;
sgtitle('Fereastra Chebyshev si Kaiser cu valori diferite ale parametrilor')
subplot(2,3,1)
stem(w4_1);
text(10, 1.2, sprintf(' r = %.2f dB', r-5));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Chebyshev');

subplot(2,3,2)
stem(w4);
text(10, 1.2, sprintf(' r = %.2f dB', r));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Chebyshev');


subplot(2,3,3)
stem(w4_2);
text(10, 1.2, sprintf(' r = %.2f dB', r+5));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Chebyshev');

subplot(2,3,4)
stem(w7_1);
text(10, 1.2, sprintf('\\beta = %.2f dB',beta-1 ));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Kaiser');

subplot(2,3,5)
stem(w7);
text(10, 1.2, sprintf('\\beta = %.2f dB', beta));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Kaiser');


subplot(2,3,6)
stem(w7_2);
text(10, 1.2, sprintf('\\beta = %.2f dB', beta+1));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Kaiser');

figure;
sgtitle('Fereastra Lanczos si Tukey cu valori diferite ale parametrilor')
subplot(2,3,1)
stem(w9_1');
text(10, 1.2, sprintf(' L = %.2f ', L-1));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Lanczos');


subplot(2,3,2)
stem(w9);
text(10, 1.2, sprintf(' L = %.2f ', L));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Lanczos');


subplot(2,3,3)
stem(w9_2);
text(10, 1.2, sprintf(' L = %.2f ', L+1));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Lanczos');

subplot(2,3,4)
stem(w8_1);
text(10, 1.2, sprintf('\\alpha = %.2f %%',alfa-15));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Tukey');

subplot(2,3,5)
stem(w8);
text(10, 1.2, sprintf('\\alpha = %.2f %%', alfa));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Tukey');

subplot(2,3,6)
stem(w8_2);
text(10, 1.2, sprintf('\\alpha = %.2f %%', alfa+15));
ylim([-0.2 1.5]);
xlabel('Indicele coeficientului');
ylabel('Valoarea coeficientului');
title('Fereastra Tukey');


w1=w1/sum(w1);
[W1,om_1]=freqz(w1,1,5000);

figure;
plot(om_1/pi,20*log10(abs(W1)));
title('Spectrul ferestrei dreptunghiulare');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');

w2=w2/sum(w2);
[W2,om_2]=freqz(w2,1,5000);

w3=w3/sum(w3);
[W3,om_3]=freqz(w3,1,5000);

w5=w5/sum(w5);
[W5,om_5]=freqz(w5,1,5000);

w6=w6/sum(w6);
[W6,om_6]=freqz(w6,1,5000);

w4=w4/sum(w4);
w4_1=w4_1/sum(w4_1);
w4_2=w4_2/sum(w4_2);

w7=w7/sum(w7);
w7_1=w7_1/sum(w7_1);
w7_2=w7_2/sum(w7_2);

w8=w8/sum(w8);
w8_1=w8_1/sum(w8_1);
w8_2=w8_2/sum(w8_2);

w9=w9/sum(w9);
w9_1=w9_1/sum(w9_1);
w9_2=w9_2/sum(w9_2);

[W4,om_4]=freqz(w4,1,5000);
[W4_1,om_4_1]=freqz(w4_1,1,5000);
[W4_2,om_4_2]=freqz(w4_2,1,5000);

[W7,om_7]=freqz(w7,1,5000);
[W7_1,om_7_1]=freqz(w7_1,1,5000);
[W7_2,om_7_2]=freqz(w7_2,1,5000);

[W8,om_8]=freqz(w8,1,5000);
[W8_1,om_8_1]=freqz(w8_1,1,5000);
[W8_2,om_8_2]=freqz(w8_2,1,5000);

[W9,om_9]=freqz(w9,1,5000);
[W9_1,om_9_1]=freqz(w9_1,1,5000);
[W9_2,om_9_2]=freqz(w9_2,1,5000);


figure;
sgtitle("Spectrele ferestrelor fara parametri");
subplot(2,2,1)
plot(om_2/pi , 20*log10(abs(W2)));
title('Spectrul ferestrei triunghiulare');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
ylim([-350 0])

subplot(2,2,2)
plot(om_3/pi,20*log10(abs(W3)));
title('Spectrul ferestrei Blackman');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
ylim([-350 0])

subplot(2,2,3)
plot(om_5/pi,20*log10(abs(W5)));
title('Spectrul ferestrei hamming');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
ylim([-350 0])

subplot(2,2,4)
plot(om_6/pi,20*log10(abs(W6)));
title('Spectrul ferestrei hanning');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
ylim([-350 0])

figure;
sgtitle('Spectrele ferestrelor Chebyshev si Kaiser cu diferite valori ale parametrilor')
subplot(2,3,1);
plot(om_4_1/pi,20*log10(abs(W4_1)));
title('Spectrul ferestrei Chebyshev');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -50, sprintf(' r = %.2f dB', r-5));
ylim([-180 0])

subplot(2,3,2);
plot(om_4/pi,20*log10(abs(W4)));
title('Spectrul ferestrei Chebyshev');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -50, sprintf(' r = %.2f dB', r));
ylim([-180 0])

subplot(2,3,3);
plot(om_4_2/pi,20*log10(abs(W4_2)));
title('Spectrul ferestrei Chebyshev');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -50, sprintf(' r = %.2f dB', r+5));
ylim([-180 0])

subplot(2,3,4)
plot(om_7_1/pi,20*log10(abs(W7_1)));
title('Spectrul ferestrei Kaiser');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -20, sprintf('\\beta = %.2f dB',beta-1 ));

subplot(2,3,5)
plot(om_7/pi,20*log10(abs(W7)));
title('Spectrul ferestrei Kaiser');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -20, sprintf('\\beta = %.2f dB',beta ));

subplot(2,3,6)
plot(om_7_2/pi,20*log10(abs(W7_2)));
title('Spectrul ferestrei Kaiser');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -20, sprintf('\\beta = %.2f dB',beta+1 ));

figure;
sgtitle('Spectrele ferestrelor Lanczos si Tukey cu diferite valori ale parametrilor')
subplot(2,3,1);
plot(om_9_1/pi,20*log10(abs(W9_1)));
title('Spectrul ferestrei Lanczos');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -50, sprintf(' L = %.2f ', L-1));
ylim([-250 0])

subplot(2,3,2);
plot(om_9/pi,20*log10(abs(W9)));
title('Spectrul ferestrei Lanczos');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -50,sprintf(' L = %.2f ', L)); 

subplot(2,3,3);
plot(om_9_2/pi,20*log10(abs(W9_2)));
title('Spectrul ferestrei Lanczos');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -50, sprintf(' L = %.2f ', L+1));

subplot(2,3,4)
plot(om_8_1/pi,20*log10(abs(W8_1)));
title('Spectrul ferestrei Tukey');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -15, sprintf('\\alpha = %.2f %%', alfa-15));
ylim([-140 0])

subplot(2,3,5)
plot(om_8/pi,20*log10(abs(W8)));
title('Spectrul ferestrei Tukey');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -15, sprintf('\\alpha = %.2f %%', alfa));
ylim([-140 0])

subplot(2,3,6)
plot(om_8_2/pi,20*log10(abs(W8_2)));
title('Spectrul ferestrei Tukey');
xlabel('Frecvenţă normalizată ');
ylabel('Amplitudine [dB]');
text(0.4, -15, sprintf('\\alpha = %.2f %%', alfa+15));

omega=linspace(0,pi,K);

h1=fir1(M-1,freq_c,boxcar(M));
h2=fir1(M-1,freq_c,triang(M));
h3=fir1(M-1,freq_c,blackman(M));
h5=fir1(M-1,freq_c);
h6=fir1(M-1,freq_c,hanning(M));

figure;
sgtitle('Secventele pondere ale ferestrelor neparametrizate')
subplot(1,5,1)
stem(h1);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Dreptunghiulara');

subplot(1,5,2)
stem(h2);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Triunghiulara');

subplot(1,5,3)
stem(h3);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Blackman');

subplot(1,5,4)
stem(h5);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Hamming');

subplot(1,5,5)
stem(h6);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Hanning');

h4=fir1(M-1,freq_c,chebwin(M,r));
h4_1=fir1(M-1,freq_c,chebwin(M,r-5));
h4_2=fir1(M-1,freq_c,chebwin(M,r+5));

h7=fir1(M-1,freq_c,kaiser(M,beta));
h7_1=fir1(M-1,freq_c,kaiser(M,beta-1));
h7_2=fir1(M-1,freq_c,kaiser(M,beta+1));

h8=fir1(M-1,freq_c,tukeywin(M,alfa/100));
h8_1=fir1(M-1,freq_c,tukeywin(M,(alfa-15)/100));
h8_2=fir1(M-1,freq_c,tukeywin(M,(alfa+15)/100));

for n=0:M-1
w9(n+1)= (sin(2*pi*((2*(n)-M+1) / (2*(M-1))))/(2*pi *((2*(n)-M+1) / (2*(M-1)))))^L;
w9_1(n+1)=(sin(2*pi*((2*(n)-M+1) / (2*(M-1))))/(2*pi *((2*(n)-M+1) / (2*(M-1)))))^(L-1);
w9_2(n+1)=(sin(2*pi*((2*(n)-M+1) / (2*(M-1))))/(2*pi *((2*(n)-M+1) / (2*(M-1)))))^(L+1);
end

w9=w9';
w9(12)=1;
w9_1=(w9_1)';
w9_1(12)=1;
w9_2=w9_2';
w9_2(12)=1;


h9=fir1(M-1,freq_c,w9);
h9_1=fir1(M-1,freq_c,w9_1);
h9_2=fir1(M-1,freq_c,w9_2);


figure;
sgtitle('Secventele pondere ale filtrelor pentru ferestrele Chebyshev si Kaiser')
subplot(2,3,1)
stem(h4_1);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Chebyshev r = 81.4324 dB')

subplot(2,3,2)
stem(h4);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Chebyshev r = 86.4324 dB')

subplot(2,3,3)
stem(h4_2);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Chebyshev r = 91.4324 dB')

subplot(2,3,4)
stem(h7_1);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Kaiser \beta = 3.5730 dB')

subplot(2,3,5)
stem(h7);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Kaiser \beta = 4.5730 dB')

subplot(2,3,6)
stem(h7_2);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Kaiser \beta = 5.5730 dB')

figure;
sgtitle('Secventele pondere ale filtrelor pentru ferestrele Lanczos si Tukey')
subplot(2,3,1)
stem(h9_1);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Lanczos L = 6')

subplot(2,3,2)
stem(h9);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Lanczos L = 7')

subplot(2,3,3)
stem(h9_2);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Lanczos L = 8')

subplot(2,3,4)
stem(h8_1);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Tukey \alpha = 24.2973%')

subplot(2,3,5)
stem(h8);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Tukey \alpha = 39.2973%')

subplot(2,3,6)
stem(h8_2);
xlabel('Indicele Coeficientului');
ylabel('Valoarea Coeficientului')
title('Tukey \alpha = 54.2973%')

[H1,omeg_1]=freqz(h1,1,K);
[H2,omeg_2]=freqz(h2,1,K);
[H3,omeg_3]=freqz(h3,1,K);
[H5,omeg_5]=freqz(h5,1,K);
[H6,omeg_6]=freqz(h6,1,K);

[H4,omeg4]=freqz(h4,1,K);
[H4_1,omeg4_1]=freqz(h4_1,1,K);
[H4_2,omeg4_2]=freqz(h4_2,1,K);

[H7,omeg7]=freqz(h7,1,K);
[H7_1,omeg7_1]=freqz(h7_1,1,K);
[H7_2,omeg7_2]=freqz(h7_2,1,K);

[H8,omeg8]=freqz(h8,1,K);
[H8_1,omeg8_1]=freqz(h8_1,1,K);
[H8_2,omeg8_2]=freqz(h8_2,1,K);

[H9,omeg9]=freqz(h9,1,K);
[H9_1,omeg9_1]=freqz(h9_1,1,K);
[H9_2,omeg9_2]=freqz(h9_2,1,K);

figure;
sgtitle('Spectrele si fazele filtrelor obtinute pentru ferestrele neparametrice')
subplot(2,5,1)
plot(omeg_1/pi,20*log10(abs(H1)));
title(' fer.dreptunghiulara');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot(2,5,6)
plot(omeg_1/pi,(angle(H1)));
title(' fer.dreptunghiulara');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,5,2)
plot(omeg_2/pi,20*log10(abs(H2)));
title('fer.triunghiulara');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot(2,5,7)
plot(omeg_2/pi,(angle(H2)));
title(' fer.triunghiulara');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,5,3)
plot(omeg_3/pi,20*log10(abs(H3)));
title('fer.blackman');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot(2,5,8)
plot(omeg_3/pi,(angle(H3)));
title('fer.blackman');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,5,4)
plot(omeg_5/pi,20*log10(abs(H5)));
title('fer.hamming');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot(2,5,9)
plot(omeg_5/pi,(angle(H5)));
title('fer.hamming');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,5,5)
plot(omeg_6/pi,20*log10(abs(H6)));
title('fer.hanning');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot(2,5,10)
plot(omeg_6/pi,(angle(H6)));
title('fer.hanning');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

figure;
sgtitle('Spectre si faze pentru filtrele realizate cu ferestrele Chebyshev si Kaiser')

subplot(3,4,1)
plot(omeg4_1/pi,20*log10(abs(H4_1)));
title('Chebyshev r = 81.4324');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-250 0]);

subplot(3,4,2)
plot(omeg4_1/pi,(angle(H4_1)));
title('Chebyshev r = 81.4324');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,3)
plot(omeg7_1/pi,20*log10(abs(H7_1)));
title('Kaiser \beta = 3.5730 ');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-250 0]);

subplot(3,4,4)
plot(omeg7_1/pi,(angle(H7_1)));
title('Kaiser \beta = 3.5730');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,5)
plot(omeg4/pi,20*log10(abs(H4)));
title('Chebyshev r=86.4324');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-250 0]);

subplot(3,4,6)
plot(omeg4/pi,(angle(H4)));
title('Chebyshev r = 86.4324');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,7)
plot(omeg7/pi,20*log10(abs(H7)));
title('Kaiser \beta = 4.5730 ');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-250 0]);

subplot(3,4,8)
plot(omeg7/pi,(angle(H7)));
title('Kaiser \beta = 4.5730');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,9)
plot(omeg4_2/pi,20*log10(abs(H4_2)));
title('Chebyshev r=91.4324');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-250 0]);

subplot(3,4,10)
plot(omeg4_2/pi,(angle(H4_2)));
title(' Chebyshev r = 91.4324');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,11)
plot(omeg7_2/pi,20*log10(abs(H7_2)));
title('Kaiser \beta = 5.5730 ');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-250 0]);

subplot(3,4,12)
plot(omeg7_2/pi,(angle(H7_2)));
title('Kaiser \beta = 5.5730');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

figure;
sgtitle('Spectre si faze pentru filtrele realizate cu ferestrele Lanczos si Tukey')

subplot(3,4,1)
plot(omeg9_1/pi,20*log10(abs(H9_1)));
title('Lanczos L=6');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-200 0]);

subplot(3,4,2)
plot(omeg9_1/pi,(angle(H9_1)));
title(' Lanczos L=6');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,3)
plot(omeg8_1/pi,20*log10(abs(H8_1)));
title(' Tukey \alpha = 24.2973 ');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-200 0]);

subplot(3,4,4)
plot(omeg8_1/pi,(angle(H8_1)));
title(' Tukey \alpha = 24.2973');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,5)
plot(omeg9/pi,20*log10(abs(H9)));
title('Lanczos L=7');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')

subplot(3,4,6)
plot(omeg9/pi,(angle(H9)));
title(' Lanczos L=7');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,7)
plot(omeg8/pi,20*log10(abs(H8)));
title(' Tukey \alpha = 39.2973');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-200 0]);

subplot(3,4,8)
plot(omeg8/pi,(angle(H8)));
title('Tukey \alpha = 39.2973');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,9)
plot(omeg9_2/pi,20*log10(abs(H9_2)));
title('Lanczos L=8');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-200 0]);

subplot(3,4,10)
plot(omeg9_2/pi,(angle(H9_2)));
title('Lanczos L=8');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(3,4,11)
plot(omeg8_2/pi,20*log10(abs(H8_2)));
title(' Tukey \alpha = 54.2973 ');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-200 0]);

subplot(3,4,12)
plot(omeg8_2/pi,(angle(H8_2)));
title(' Tukey \alpha = 54.2973');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })




hkaiser1=fir1(34,freq_c,kaiser(35,beta-1));
hkaiser2=fir1((2*M)-1,freq_c,kaiser((2*M),beta-1));
hcheb1=fir1(34,freq_c,chebwin(35,r-5));
hcheb2=fir1(2*M-1,freq_c,chebwin(2*M,r-5));
[HK1,omega_k1]=freqz(hkaiser1,1,K);
[HK2,omega_k2]=freqz(hkaiser2,1,K);
[HC1,omega_c1]=freqz(hcheb1,1,K);
[HC2,omega_c2]=freqz(hcheb2,1,K);


figure;
sgtitle('Faza 2 subpunctul b pentru #1 din clasament')
subplot(2,3,1)
plot(omeg7_1/pi,20*log10(abs(H7_1)));
title('Kaiser \beta = 3.5730 ordin M');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot (2,3,2)
plot(omega_k1/pi,20*log10(abs(HK1)));
title('Kaiser \beta = 3.5730 ordin M+M/2');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot(2,3,3)
plot(omega_k2/pi,20*log10(abs(HK2)));
title('Kaiser \beta = 3.5730 ordin 2M');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-160 0]);

subplot(2,3,4)
plot(omeg7_2/pi,(angle(H7_2)));
title('Kaiser \beta = 3.5730 ordin M');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,3,5)
plot(omega_k1/pi,(angle(HK1)));
title('Kaiser \beta = 3.5730 ordin M+M/2');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,3,6)
plot(omega_k2/pi,(angle(HK2)));
title('Kaiser \beta = 3.5730 ordin 2M');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })




figure;
sgtitle('Faza 2 subpunctul b pentru #9 din clasament')
subplot(2,3,1)
plot(omeg4_1/pi,20*log10(abs(H4_1)));
title('Chebyshev r = 81.4324 ordin M');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-180 0]);

subplot(2,3,2)
plot(omega_c1/pi,20*log10(abs(HC1)));
title('Chebyshev r = 81.4324 ordin M+M/2');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-180 0]);

subplot(2,3,3)
plot(omega_c2/pi,20*log10(abs(HC2)));
title('Chebyshev r = 81.4324 ordin 2M');
xlabel('Frecventa normalizata')
ylabel('Amplitudine(dB)')
ylim([-180 0]);

subplot(2,3,4)
plot(omeg4_1/pi,(angle(H4_1)));
title('Chebyshev r = 81.4324 ordin M');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,3,5)
plot(omega_c1/pi,(angle(HC1)));
title('Chebyshev r = 81.4324 ordin M+M/2');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })

subplot(2,3,6)
plot(omega_c2/pi,(angle(HC2)));
title('Chebyshev r = 81.4324 ordin 2M');
xlabel('Frecventa normalizata')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })



