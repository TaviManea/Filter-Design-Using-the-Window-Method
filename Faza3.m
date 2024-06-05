
clc;
omega_p = 1.2346;

omega_s = 1.5521;

Delta_p = 3.9297;

Delta_s = 3.9297;

beta=4.5730;

Delta_pn=Delta_p/100;
Delta_pn1=20*log10(1+Delta_pn);
Delta_pn2=20*log10(1-Delta_pn);
Delta_sn=20*log10(Delta_s/100);

M=23;
M1=35;
M2=45;
K=5000;

puls_c1=0.75*omega_p+0.25*omega_s;
puls_c2=0.5*omega_p+0.5*omega_s;
puls_c3=0.25*omega_p+0.75*omega_s;

freq_c1=puls_c1/pi;
freq_c2=puls_c2/pi;
freq_c3=puls_c3/pi;

h1=fir1(M-1,freq_c1,kaiser(M,beta-1));
h2=fir1(M-1,freq_c2,kaiser(M,beta-1));
h3=fir1(M-1,freq_c3,kaiser(M,beta-1));
h4=fir1(M1-1,freq_c1,kaiser(M1,beta-1));
h5=fir1(M1-1,freq_c2,kaiser(M1,beta-1));
h6=fir1(M1-1,freq_c3,kaiser(M1,beta-1));
h7=fir1(M2-1,freq_c1,kaiser(M2,beta-1));
h8=fir1(M2-1,freq_c2,kaiser(M2,beta-1));
h9=fir1(M2-1,freq_c3,kaiser(M2,beta-1));

 [H1,omega1]=freqz(h1,1,K);
 [H2,omega2]=freqz(h2,1,K);
 [H3,omega3]=freqz(h3,1,K);
 [H4,omega4]=freqz(h4,1,K);
 [H5,omega5]=freqz(h5,1,K);
 [H6,omega6]=freqz(h6,1,K);
 [H7,omega7]=freqz(h7,1,K);
 [H8,omega8]=freqz(h8,1,K);
 [H9,omega9]=freqz(h9,1,K);


 figure;
 sgtitle('Secventele pondere ale filtrelor in ordine descrescatoare ca performanta')
 subplot(3,3,8)
 stem(h1);
 title('Ordinul M = 23 ');
 subplot(3,3,6)
 stem(h2);
 title('Ordinul M = 23');
 subplot(3,3,9)
 stem(h3);
 title('Ordinul M = 23');
 subplot(3,3,7)
 stem(h4);
 title('Ordinul M = 35');
 subplot(3,3,5)
 stem(h5);
 title('Ordinul M = 35');
 subplot(3,3,3)
 stem(h6);
 title('Ordinul M = 35');
 subplot(3,3,2)
 stem(h7);
 title('Ordinul M = 45');
 subplot(3,3,4)
 stem(h8);
 title('Ordinul M = 45');
 subplot(3,3,1)
 stem(h9);
 title('Ordinul M = 45');

figure;
sgtitle('Spectrele filtrelor')
subplot(3,3,1)
plot(omega9,20*log10(abs(H9)))
title('M=45 \Delta_{pr} = 17.02% , \Delta_{sr} = 0.58%')
xline([omega_p puls_c3 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,2)
plot(omega7,20*log10(abs(H7)))
title('M=45 \Delta_{pr} = 17.19% , \Delta_{sr} = 0.55%')
xline([omega_p puls_c1 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,3)
plot(omega6,20*log10(abs(H6)))
title('M=35 \Delta_{pr} = 1.07% , \Delta_{sr} = 23.53%')
xline([omega_p puls_c3 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,4)
plot(omega8,20*log10(abs(H8)))
title('M=45 \Delta_{pr} = 1.03% , \Delta_{sr} = 1.4%')
xline([omega_p puls_c2 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,5)
plot(omega5,20*log10(abs(H5)))
title('M=35 \Delta_{pr} = 6.15% , \Delta_{sr} = 6.32%')
xline([omega_p puls_c2 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,6)
plot(omega2,20*log10(abs(H2)))
title('M=23 \Delta_{pr} = 17.78% , \Delta_{sr} = 17.48%')
xline([omega_p puls_c2 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,7)
plot(omega4,20*log10(abs(H4)))
title('M=35 \Delta_{pr} = 23.62% , \Delta_{sr} = 0.8%')
xline([omega_p puls_c1 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,8)
plot(omega1,20*log10(abs(H1)))
title('M=23 \Delta_{pr} = 32.07% , \Delta_{sr} = 7.16%')
xline([omega_p puls_c1 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');

subplot(3,3,9)
plot(omega3,20*log10(abs(H3)))
title('M=23 \Delta_{pr} =7.52 , \Delta_{sr} = 32.09%%')
xlabel('Frecventa');
ylabel('Amplitudine(dB)');
xline([omega_p puls_c3 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');


figure;
sgtitle('Fazele filtrelor')
subplot(3,3,1)
plot(omega9,angle(H9))
xline([omega_p puls_c3 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c3 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=45 \Delta_{pr} = 17.02% , \Delta_{sr} = 0.58%')

subplot(3,3,2)
plot(omega7,angle(H7))
xline([omega_p puls_c1 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c1 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=45 \Delta_{pr} = 17.19% , \Delta_{sr} = 0.55%')

subplot(3,3,3)
plot(omega6,angle(H6))
xline([omega_p puls_c3 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c3 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=35 \Delta_{pr} = 1.07% , \Delta_{sr} = 23.53%')

subplot(3,3,4)
plot(omega8,angle(H8))
xline([omega_p puls_c2 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c2 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=45 \Delta_{pr} = 1.03% , \Delta_{sr} = 1.4%')

subplot(3,3,5)
plot(omega5,angle(H5))
xline([omega_p puls_c2 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c2 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=35 \Delta_{pr} = 6.15% , \Delta_{sr} = 6.32%')

subplot(3,3,6)
plot(omega2,angle(H2))
xline([omega_p puls_c2 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c2 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=23 \Delta_{pr} = 17.78% , \Delta_{sr} = 17.48%')

subplot(3,3,7)
plot(omega4,angle(H4))
xline([omega_p puls_c1 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c1 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=35 \Delta_{pr} = 23.62% , \Delta_{sr} = 0.8%')

subplot(3,3,8)
plot(omega1,angle(H1))
xline([omega_p puls_c1 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c1 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=23 \Delta_{pr} = 32.07% , \Delta_{sr} = 7.16%')

subplot(3,3,9)
plot(omega3,angle(H3))
xline([omega_p puls_c3 omega_s],'red');
xlabel('Frecventa')
ylabel('Faza(rad)')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi' ,'-\pi/2',0,'\pi/2','\pi' })
xticks([omega_p puls_c3 omega_s ]);
xticklabels({'\omega_{p}','\omega_{c}','\omega_{s}'})
title('M=23 \Delta_{pr} =7.52 , \Delta_{sr} = 32.09%%')




 [Delta_pr1,Delta_sr1]=Functie(h1,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr1);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr1);

 [Delta_pr2,Delta_sr2]=Functie(h2,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr2);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr2);

 [Delta_pr3,Delta_sr3]=Functie(h3,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr3);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr3);

 [Delta_pr4,Delta_sr4]=Functie(h4,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr4);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr4);

 [Delta_pr5,Delta_sr5]=Functie(h5,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr5);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr5);

 [Delta_pr6,Delta_sr6]=Functie(h6,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr6);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr6);

 [Delta_pr7,Delta_sr7]=Functie(h7,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr7);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr7);

 [Delta_pr8,Delta_sr8]=Functie(h8,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr8);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr8);

 [Delta_pr9,Delta_sr9]=Functie(h9,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr9);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr9);




function [Delta_pr,Delta_sr] = Functie(h,omega_p,omega_s)


[H,omega]=freqz(h,1,1000);

 om_p= find(omega <= omega_p);
 om_s=find(omega >= omega_s);

 H_p=abs(H(om_p));
 H_s=abs(H(om_s));

 Delta_pr=max(abs(1-abs(H_p)));
 Delta_sr=max(abs(H_s));


Delta_pr=Delta_pr*100;
Delta_sr=Delta_sr*100;

end





