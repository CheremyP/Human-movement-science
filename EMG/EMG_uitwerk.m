clear;
% close all


%% OPDRACHT 1 %%

load kalibreren.mat

x = [0 48.3672 100.2916 148.3625 195.9704 245.3716];
y = [mean(nulmeting) mean(vijfkg) mean(tienkg) mean(vijftienkg) mean(twintigkg) mean(vijfentwintigkg)];
coef = polyfit(x,y,1);
kalibratie_y = coef(1)*x + coef(2);
figure(1)
plot(x,kalibratie_y,'-o')
title('Kalibratielijn')
xlabel('Krachten [N]')
ylabel('Spanning [mV]')

% Check formule spanning_kracht

%% OPDRACHT 3 %% 

% close all
load MVC_flexoren.mat
load MVC_extensoren.mat
MVC_flexoren = double(MVC_flexoren);
MVC_extensoren = double(MVC_extensoren);

N = length(MVC_flexoren);
k = [0:N-1];
fs = 2000;
dt = 1/fs;
t_flexoren = k*dt;

N1 = length(MVC_extensoren);
k1 = [0:N1-1];
fs1 = 2000;
dt1 = 1/fs1;
t_extensoren = k1*dt1;

figure(2)
subplot(2,1,1)
plot(t_flexoren,MVC_flexoren)
title('MVC flexoren')
xlabel('Tijd [s]')
ylabel('Spanning [mV]')
legend('krachtopnemer', 'biceps' ,'brachioradialis' ,'triceps lange kop','triceps laterale kop')
hold on
subplot(2,1,2)
plot(t_extensoren,MVC_extensoren)
title('MVC extensoren')
xlabel('Tijd [s]')
ylabel('Spanning [mV]')
legend('krachtopnemer', 'biceps' ,'brachioradialis' ,'triceps lange kop','triceps laterale kop')


[PSD,f] = pwelch(MVC_flexoren-mean(MVC_flexoren),[],[],[],fs);
figure(3)
plot(f,PSD)
ylim([0 125000])

Wn1 = [10 , 500] /(fs/2);
[B,A] = butter(2,Wn1);
MVC_flexoren_filt = filtfilt(B,A,MVC_flexoren);

figure(4)
plot(t_flexoren,MVC_flexoren_filt(:,2:3),t_flexoren,MVC_flexoren(:,2:3))
legend('biceps gefilterd' ,'brachioradialis gefilterd', 'biceps' ,'brachioradialis')

Wn2 = [10 , 500] /(fs1/2);
[B,A] = butter(2,Wn2);
MVC_extensoren_filt = filtfilt(B,A,MVC_extensoren);

figure(5)
plot(t_extensoren,MVC_extensoren_filt(:,4:5),t_extensoren,MVC_extensoren(:,4:5))

legend('tricep lang gefilterd' ,'triceps lateraal gefilterd', 'triceps lang' ,'triceps lateraal')


MVC_flexoren_rv = abs(MVC_flexoren_filt(:,2:3));
MVC_flexoren_arv = mean(abs(MVC_flexoren_filt(:,2:3)));

MVC_extensoren_rv = abs(MVC_extensoren_filt(:,4:5));
MVC_extensoren_arv = mean(abs(MVC_extensoren_filt(:,4:5)));

% close all

a = ones([1,length(t_flexoren)]);
a = a';
MVC_flexoren_arv = MVC_flexoren_arv .* a;

b = ones([1,length(t_extensoren)]);
b = b';
MVC_extensoren_arv = MVC_extensoren_arv .* b;


for i = 1:N-999
    
    MVC_flexoren_arv_biceps(i) = mean(abs(MVC_flexoren_filt(i:i+999,2)));
    MVC_flexoren_arv_brachioradialis(i) = mean(abs(MVC_flexoren_filt(i:i+999,3)));
       
end

for i = 1:N1-999
     MVC_extensoren_arv_triceplang(i) = mean(abs(MVC_extensoren_filt(i:i+999,4)));
    MVC_extensoren_arv_triceplateraal(i) = mean(abs(MVC_extensoren_filt(i:i+999,5)));
end


MVC_flexoren_arv_biceps_max = max(MVC_flexoren_arv_biceps);
MVC_flexoren_arv_biceps_max_lijn = MVC_flexoren_arv_biceps_max .*a;
MVC_flexoren_arv_brachioradialis_max = max(MVC_flexoren_arv_brachioradialis);
MVC_flexoren_arv_brachioradialis_max_lijn = MVC_flexoren_arv_brachioradialis_max .*a;
MVC_extensoren_arv_triceplang_max = max(MVC_extensoren_arv_triceplang);
MVC_extensoren_arv_triceplang_max_lijn = MVC_extensoren_arv_triceplang_max .* b;
MVC_extensoren_arv_triceplateraal_max = max(MVC_extensoren_arv_triceplateraal);
MVC_extensoren_arv_triceplateraal_max_lijn = MVC_extensoren_arv_triceplateraal_max.* b;


figure (6)
plot(t_flexoren, MVC_flexoren_rv, t_flexoren,MVC_flexoren_arv,t_flexoren,MVC_flexoren_arv_biceps_max_lijn,t_flexoren,MVC_flexoren_arv_brachioradialis_max_lijn)
legend('biceps gefilterd' ,'brachioradialis gefilterd', 'biceps gemiddeld alles' ,'brachioradialis gemiddeld alles' , 'biceps gemiddeld max' , 'brachioradialis gemiddeld max')

figure (7)
plot(t_extensoren, MVC_extensoren_rv, t_extensoren,MVC_extensoren_arv,t_extensoren,MVC_extensoren_arv_triceplang_max_lijn,t_extensoren,MVC_extensoren_arv_triceplateraal_max_lijn)
legend('tricep lang gefilterd' ,'tricep lateraal gefilterd', 'tricep lang gemiddeld alles' ,'tricep lateraal gemiddeld alles' , 'tricep lang gemiddeld max' , 'tricep lateraal gemiddeld max')

%% OPDRACHT 4 %%

% close all
load opdr4_05Hz.mat
load opdr4_2Hz.mat

N2 = length(opdr4_05Hz);
k2 = [0:N2-1];
dt2 = 1/fs;
t2 = k2*dt2;
figure(8)
subplot(2,1,1)
plot(t2,opdr4_05Hz(:,2:3),t2,opdr4_05Hz(:,5))
title('0.5 Hz')
xlabel('Tijd [s]')
ylabel('Spanning [mV]')
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal')
hold on
N3 = length(opdr4_2Hz);
k3 = [0:N3-1];
dt3 = 1/fs;
t3 = k3*dt3;
subplot(2,1,2)
plot(t3,opdr4_2Hz(:,2:3),t3,opdr4_2Hz(:,5))
title('2 Hz')
xlabel('Tijd [s]')
ylabel('Spanning [mV]')
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal')

[PSD1,f1] = pwelch(opdr4_05Hz-mean(opdr4_05Hz),[],[],[],fs);
figure(9)
plot(f1,PSD1(:,2:3),f1,PSD1(:,5))
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal') 
ylim([0 20000])

[PSD2,f2] = pwelch(opdr4_2Hz-mean(opdr4_2Hz),[],[],[],fs);
figure(10)
plot(f2,PSD2(:,2:3),f2,PSD2(:,5))
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal') 

opdr4_05Hz = double(opdr4_05Hz);
opdr4_2Hz = double(opdr4_2Hz);
[B1,A1] = butter(2,[10/(fs/2),500/(fs/2)]);
opdr4_05Hz_filt = filtfilt(B1,A1,opdr4_05Hz);
opdr4_2Hz_filt = filtfilt(B1,A1,opdr4_2Hz);

[PSD3,f3] = pwelch(opdr4_05Hz_filt,[],[],[],fs);
figure(11)
plot(f3,PSD3(:,2:3),f3,PSD3(:,5))
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal') 

[PSD4,f4] = pwelch(opdr4_2Hz_filt,[],[],[],fs);
figure(12)
plot(f4,PSD4(:,2:3),f4,PSD4(:,5))
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal') 

% close all

opdr4_05Hz_gelijk_norm = [abs(opdr4_05Hz_filt(:,2))/MVC_flexoren_arv_biceps_max abs(opdr4_05Hz_filt(:,3))/MVC_flexoren_arv_brachioradialis_max abs(opdr4_05Hz_filt(:,5))/MVC_extensoren_arv_triceplateraal_max]*100;
figure(13)
plot(t2,opdr4_05Hz_gelijk_norm)
title('Genormaliseerde gelijkgerichte 0.5 Hz')
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal') 

[B2,A2] = butter(2,4/(fs/2));
opdr4_05Hz_gelijk_norm_omhul = filtfilt(B2,A2,opdr4_05Hz_gelijk_norm);
figure(14)
plot(t2,opdr4_05Hz_gelijk_norm_omhul)

% opdr4_05Hz_gelijk_omhul2 = filtfilt(B2,A2,abs(opdr4_05Hz_filt));
% opdr4_05Hz_gelijk_norm_omhul2 = [opdr4_05Hz_gelijk_omhul2(:,2)/MVC_flexoren_arv_biceps_max opdr4_05Hz_gelijk_omhul2(:,3)/MVC_flexoren_arv_brachioradialis_max opdr4_05Hz_gelijk_omhul2(:,5)/MVC_extensoren_arv_triceplateraal_max]*100;
% figure(15)
% plot(t2,opdr4_05Hz_gelijk_norm_omhul2)

opdr4_2Hz_gelijk_norm = [abs(opdr4_2Hz_filt(:,2))/MVC_flexoren_arv_biceps_max abs(opdr4_2Hz_filt(:,3))/MVC_flexoren_arv_brachioradialis_max abs(opdr4_2Hz_filt(:,5))/MVC_extensoren_arv_triceplateraal_max]*100;
figure(16)
plot(t3,opdr4_2Hz_gelijk_norm)
title('Genormaliseerde gelijkgerichte 2 Hz')
legend('Biceps' , 'Brachioradialis' , 'Triceps lateraal') 

opdr4_2Hz_gelijk_norm_omhul = filtfilt(B2,A2,opdr4_2Hz_gelijk_norm);
figure(17)
plot(t3,opdr4_2Hz_gelijk_norm_omhul)

% [B3,A3] = butter(2,2/(fs/2));
% opdr4_05Hz_gelijk_norm_omhul3 = filtfilt(B3,A3,opdr4_05Hz_gelijk_norm);
% figure(18)
% plot(t2,opdr4_05Hz_gelijk_norm_omhul3)

opdr4_05Hz_kracht = spanning_kracht(opdr4_05Hz(:,1));
opdr4_2Hz_kracht = spanning_kracht(opdr4_2Hz(:,1));
figure(19)
subplot(2,1,1)
plot(t2,opdr4_05Hz_kracht)
title('Kracht 0.5 Hz')
xlabel('Tijd [s]') 
ylabel('Kracht [N]')
hold on
subplot(2,1,2)
plot(t3,opdr4_2Hz_kracht)
title('Kracht 2 Hz')
xlabel('Tijd [s]') 
ylabel('Kracht [N]')

% close all

figure (20)
yyaxis left
plot(t2,opdr4_05Hz_kracht)
yyaxis right
plot(t2,opdr4_05Hz_gelijk_norm_omhul)
title('kracht en EMG 0.5 Hz')

figure (21)
yyaxis left
plot(t3,opdr4_2Hz_kracht)
yyaxis right
plot(t3,opdr4_2Hz_gelijk_norm_omhul)
title('kracht en EMG 2 Hz')

% je ziet dat het verloop veel op elkaar lijkt, maar je hebt een kleine
% vertraging in het krachtsignaal door EMD
% close all
[c,lags] = xcov(opdr4_05Hz_kracht,opdr4_05Hz_gelijk_norm_omhul(:,3));
figure (22)
tau = lags*dt2;
ishow = find(tau>=-1 & tau<=1);
plot(tau(ishow),c(ishow))

[kruiscorrmax , imax] = max(c);
EMD_brachioradialis_05 = tau(imax);

[c1,lags1] = xcov(opdr4_2Hz_kracht,opdr4_2Hz_gelijk_norm_omhul(:,3));
figure (23)
tau1 = lags1*dt2;
ishow1 = find(tau1>=-1 & tau1<=1);
plot(tau1(ishow1),c1(ishow1))

[kruiscorrmax1 , imax1] = max(c1);
EMD_brachioradialis_2 = tau1(imax1);



%% OPDRACHT 5 %%

load opdr5.mat



opdr5 = double(opdr5);

N5 = length(opdr5);
k5 = [0:N5-1];
dt5 = 1/fs;
t5 = k5*dt5;

kracht_opdr5 = spanning_kracht(opdr5(:,1));


[PSD5,f5] = pwelch(opdr5(:,2)-mean(opdr5(:,2)),[],[],[],fs);
figure(50)
plot(f5,PSD5)
ylim([0 8000])

Wn5 = [10 , 500] /(fs/2);
[B5,A5] = butter(2,Wn5);
opdr5_biceps_filt = filtfilt(B5,A5,opdr5(:,2));

figure (51)
subplot(2,1,1)
plot(t5,opdr5(:,2))
hold on
subplot(2,1,2)
plot(t5,opdr5_biceps_filt)

%  power = opdr5(:,2).^2;

mediaanfreq = [];
gemfreq = [];
opdr5_eff = [];
opdr5_arv = [];

for j = 1:20000:240000
   % j = 1:12
   %[PSD5,f] = pwelch(opdr5(1+((j-1)*20000):j*20000,2)-mean(opdr5(1+((j-1)*20000):j*20000,2)),[],[],[],fs);
    
    [PSD5,f5] = pwelch(opdr5_biceps_filt(j:j+19999)-mean(opdr5_biceps_filt(j:j+19999)),[],[],[],fs);

    mediaanfreq = [mediaanfreq f5(sum(cumsum(PSD5)<=sum(PSD5)/2))];
   
    gemfreq = [gemfreq (sum(f5.*PSD5)/sum(PSD5))];
    
    
    opdr5_eff = [opdr5_eff rms(opdr5_biceps_filt(j:j+19999))];
    opdr5_arv = [opdr5_arv mean(abs(opdr5_biceps_filt(j:j+19999)))];

end

figure (52)
subplot(1,2,1)
bar(mediaanfreq);
subplot(1,2,2)
bar(gemfreq);

figure (53)
subplot(1,2,1)
bar(opdr5_eff);
subplot(1,2,2)
bar(opdr5_arv);

close all


%% OPDRACHT 6 %% 

load opdr6.mat
opdr6_lang = double(opdr6_lang);
opdr6_beide = double(opdr6_beide);

N6 = length(opdr6_lang);
k6 = [0:N6-1];
dt6 = 1/fs;
t6 = k6*dt6;

N7 = length(opdr6_beide);
k7 = [0:N7-1];
dt7 = 1/fs;
t7 = k7*dt7;


figure(60)
subplot(2,1,1)
plot(t6,opdr6_lang)
subplot(2,1,2)
plot(t7,opdr6_beide)


Wn6 = [10 , 500] /(fs1/2);
[B6,A6] = butter(2,Wn6);

opdr6_lang_filt = filtfilt(B6,A6,opdr6_lang);
opdr6_beide_filt = filtfilt(B6,A6,opdr6_beide);

figure(61)
subplot(2,1,1)
plot(t6,opdr6_lang_filt)
subplot(2,1,2)
plot(t7,opdr6_beide_filt)

opdr6_lang_lang_norm = abs(opdr6_lang_filt(:,1)/MVC_extensoren_arv_triceplang_max);
opdr6_lang_lateraal_norm = abs(opdr6_lang_filt(:,2)/MVC_extensoren_arv_triceplateraal_max);

opdr6_beide_lang_norm = abs(opdr6_beide_filt(:,1)/MVC_extensoren_arv_triceplang_max);
opdr6_beide_lateraal_norm = abs(opdr6_beide_filt(:,2)/MVC_extensoren_arv_triceplateraal_max);

figure(62)
subplot(2,1,1)
plot(t6,opdr6_lang_lang_norm, t6 ,opdr6_lang_lateraal_norm)
subplot(2,1,2)
plot(t7,opdr6_beide_lang_norm, t7,opdr6_beide_lateraal_norm)

opdr6_lang_omhul = [filtfilt(B2,A2,opdr6_lang_lang_norm) filtfilt(B2,A2,opdr6_lang_lateraal_norm)] ;
opdr6_beide_omhul = [filtfilt(B2,A2,opdr6_beide_lang_norm) filtfilt(B2,A2,opdr6_beide_lateraal_norm)] ;

figure(63)
subplot(2,1,1)
plot(t6,opdr6_lang_omhul)
subplot(2,1,2)
plot(t7,opdr6_beide_omhul)

opdr6_lang_gem = mean(opdr6_lang_omhul);
opdr6_beide_gem = mean(opdr6_beide_omhul);


opdr6_verschil_1 = opdr6_lang_gem(1) - opdr6_lang_gem(2);
opdr6_verschil_2 = opdr6_beide_gem(1) - opdr6_beide_gem(2);

