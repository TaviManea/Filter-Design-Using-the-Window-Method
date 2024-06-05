omega_p = 1.2346;
omega_s = 1.5521;
Delta_p = 3.9297;
Delta_s = 3.9297;

beta=2;

M=29;

omega_c=1.39369;

M1=31;

omega_t1=1.4;

M2=33;
beta2=2;
omega_t2=1.395;

w1=kaiser(M,beta);
w2=kaiser(M1,beta);
w3=kaiser(M2,beta);

Delta_pn=Delta_p/100;
Delta_pn1=20*log10(1+Delta_pn);
Delta_pn2=20*log10(1-Delta_pn);
Delta_sn=20*log10(Delta_s/100);


h1 = fir1(M-1,omega_c/pi,kaiser(M,beta));
h2=fir1(M1-1,omega_t1/pi,kaiser(M1,beta));
h3=fir1(M2-1,omega_t2/pi,kaiser(M2,beta));

[H1,omega1]=freqz(h1,1,5000);
[H2,omega2]=freqz(h2,1,5000);
[H3,omega3]=freqz(h3,1,5000);


[Delta_pr1,Delta_sr1]=Functie(h1,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr1);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr1);
[Delta_pr1,Delta_sr1]=Functie(h2,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr1);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr1);
[Delta_pr1,Delta_sr1]=Functie(h3,omega_p,omega_s);
 fprintf('Toleranta realizata in banda de trecere: %.2f%%\n', Delta_pr1);
 fprintf('Toleranta realizata in banda de stopare: %.2f%%\n', Delta_sr1);




save MANEA_Octavian_F#4 h1 omega_p omega_c omega_s Delta_p Delta_s

figure;
subplot(3,3,1)
plot(omega1,20*log10(abs(H1)));
title("M=29,\omega_{c}=1.39369,kaiser,\beta=2");
xline([omega_p omega_c omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');
subplot(3,3,2)
plot(omega2,20*log10(abs(H2)));
title("M=31,\omega_{c}=1.4,kaiser,\beta=2");
xline([omega_p omega_t1 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');
subplot(3,3,3)
plot(omega1,20*log10(abs(H1)));
title("M=33,\omega_{c}=1.395,kaiser,\beta=2");
xline([omega_p omega_t2 omega_s],'red');
line([0,omega_p],[Delta_pn1,Delta_pn1],'Color','green');
line([0,omega_p],[Delta_pn2,Delta_pn2],'Color','green');
line([omega_s,pi],[Delta_sn,Delta_sn],'Color','green');
xlabel('Frecventa');
ylabel('Amplitudine(dB)');
subplot(3,3,4)
plot(omega1,angle(H1))
title('Faza Gold')
xlabel('Frecventa')
ylabel('Faza(rad)')
subplot(3,3,5)
plot(omega2,angle(H2))
title('Faza Silver')
xlabel('Frecventa')
ylabel('Faza(rad)')

subplot(3,3,6)
plot(omega3,angle(H3))
title('Faza Bronze')
xlabel('Frecventa')
ylabel('Faza(rad)')

subplot(3,3,7)
stem(h1)
title('Secventa pondere Gold')
subplot(3,3,8)
stem(h2)
title('Secventa pondere Silver')

subplot(3,3,9)
stem(h3)
title('Secventa pondere Bronze')


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