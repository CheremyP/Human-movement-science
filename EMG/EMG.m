%%%Uitwerkingen EMG
clear
close all
%% bestanden laden
load('Kalibratie1_0kg-01-26Feb2020at09h58m18s.mat')
kalibratie0kg = DataMinute001;

load('Kalibratie2_5kg-01-26Feb2020at10h02m42s.mat')
kalibratie5kg = DataMinute001;

load('Kalibratie3_10kg-02-26Feb2020at10h04m13s.mat')
kalibratie10kg = DataMinute001;

load('Kalibratie4_15kg-04-26Feb2020at10h08m09s.mat')
kalibratie15kg = DataMinute001;

load('Kalibratie5_20kg-05-26Feb2020at10h10m26s.mat')
kalibratie20kg = DataMinute001;
 
load('Kalibratie6_25kg-06-26Feb2020at10h12m27s.mat')
kalibratie25kg = DataMinute001;
 
%% dubbele matix maken
% De datamatrixen staan in single precisions getallen, wat inhoudt dat ze
% maar 4 bytes aan opslag vragen. Echter gaan we bij de uitwerkingen
% functies gebruiken die double precisions getallen nodig hebben, die 8
% bytes aan opslag vragen. We moeten dus de data omzetten naar double
% precisions getallen dmv de functie double. 
kalibratie0kg = double(kalibratie0kg);
kalibratie5kg = double(kalibratie5kg);
kalibratie10kg = double(kalibratie10kg);
kalibratie15kg = double(kalibratie15kg);
kalibratie20kg = double(kalibratie20kg);
kalibratie25kg = double(kalibratie25kg);
 
%% Opracht 1: kalibratielijn van de krachtopmeter

%vector massa maken
%touwtje om de gewichten vast te maken was 30.4 gram
m1 = 0.0;
m2 = 0.0304 + 4.982;
m3 = 0.0304 + 9.650;
m4 = 0.0304 + 4.982 + 4.941 + 4.858;
m5 = 0.0304 + 4.982 + 4.941 + 4.858 + 4.828;
m6 = 0.0304 + 24.982;

%Kracht vector maken
m = [m1 m2 m3 m4 m5 m6];
g = 9.81;
F = m*g;

%spanning vector maken
spanning0kg = mean(kalibratie0kg(:,1));
spanning5kg = mean(kalibratie5kg(:,1));
spanning10kg = mean(kalibratie10kg(:,1));
spanning15kg = mean(kalibratie15kg(:,1));
spanning20kg = mean(kalibratie20kg(:,1));
spanning25kg = mean(kalibratie25kg(:,1));
spanning = [spanning0kg spanning5kg spanning10kg spanning15kg spanning20kg spanning25kg];

%kalibratielijn maken
coef = polyfit(F,spanning,1);
Y = coef(1)*F + coef(2);
% Y = -11.0153*F -15.6063

% plot maken kalibratielijn met gebruikte data punten
figure(1)
plot(F,Y)
hold on
plot(F,spanning, 'o')
title('Kalibratielijn')
xlabel('Kracht (N)')
ylabel('Spanning(mV)')
legend('lijn', 'gebruikte data')
% ja er is spraken van een ofset zo groot als de waarde van coef(2) dus
% -15.6063

% kalibratie lijn omzetten
F_spier = (Y-coef(2))/coef(1);

%% Opdracht 3: Maximale vrijwillige contractie

%Data laden
%arm flexie
load('opdracht3_MVC_Bibr_Brac2-10-26Feb2020at10h58m58s.mat')
BibrBrac = DataMinute001;
%arm extentie
load('opdracht3_MVC_TRLO_TRLA2-15-26Feb2020at11h04m58s.mat')
TrloTrla = DataMinute001;

%Dubbele matrix maken
BibrBrac = double(BibrBrac);
TrloTrla = double(TrloTrla);

%Nu heb je de data van alle poorten, poort 1 is de kracht poort 2 de
%bicepbrachii poort 3 de brachioradialis poort 4 de tricep brachii lange
%kop en poort 5 de tricep brachii laterale kop

Bicepbrachii = BibrBrac(:,2);
Brachioradialis = BibrBrac (:,3);
Tricepsbrachii_lange_kop = TrloTrla(:,4);
Tricepsbrachii_laterale_kop = TrloTrla(:,5);

%tijdas maken
fs = 2000;
dt = 1/fs;
N1 = length(Bicepbrachii);
k1 = [0:N1-1];
t1 = dt*k1;

%plot de gemeten EMG-signalen tegen de tijd
figure(2)
plot(t1, Bicepbrachii)
hold on 
plot(t1, Brachioradialis)
title('MVC onbewerkt')
xlabel('tijd (s)')
ylabel('Spanning (mV)')
legend('Bicepbrachii', 'Brachioradialis')

%tijdas maken omdat de signalen niet dezelfde lengte hebben
fs = 2000;
dt = 1/fs;
N2 = length(Tricepsbrachii_lange_kop);
k2 = [0:N2-1];
t2 = dt*k2;

figure(3)
plot(t2, Tricepsbrachii_lange_kop)
hold on 
plot(t2, Tricepsbrachii_laterale_kop)
title('MVC onbewerkt')
xlabel('tijd (s)')
ylabel('Spanning (mV)')
legend('Tricepsbrachii lange kop', 'Tricepsbrachii laterale kop')


%gegevens filteren
% hoogdoorlaat filter met lage fc (20) tegen bewegingsartefacten
% laagdoorlaat filter met fs/5 afsnijfrequentie als anti-aliasing filter
N = 2;
fc = [20 250];
Wn = fc/(fs/2);
[B, A] = butter(N,Wn,'bandpass');

Bicepbrachii_filt = filtfilt(B,A,Bicepbrachii);
Brachioradialis_filt = filtfilt(B,A,Brachioradialis);
Tricepsbrachii_lange_kop_filt = filtfilt(B,A,Tricepsbrachii_lange_kop);
Tricepsbrachii_laterale_kop_filt = filtfilt(B,A,Tricepsbrachii_laterale_kop);

%plot gefilterde signaal in zelfde plot als onbewerkt signaal
figure(4)
plot(t1, Bicepbrachii)
hold on 
plot(t1, Brachioradialis)
hold on 
plot(t1,Bicepbrachii_filt)
hold on 
plot(t1,Brachioradialis_filt)
title('MVC')
xlabel('tijd (s)')
ylabel('Spanning (mV)')
legend('Bicepbrachii onbewerkt', 'Brachioradialis onbewerkt','Bicepbrachii bewerkt', 'Brachioradialis bewerkt')

figure(5)
plot(t2, Tricepsbrachii_lange_kop)
hold on 
plot(t2, Tricepsbrachii_laterale_kop)
hold on 
plot(t2,Tricepsbrachii_lange_kop_filt)
hold on 
plot(t2,Tricepsbrachii_laterale_kop_filt)
title('MVC')
xlabel('tijd (s)')
ylabel('Spanning (mV)')
legend('Tricepsbrachii lange kop onbewerkt', 'Tricepsbrachii laterale kop onbewerkt', 'Tricepsbrachii lange kop bewerkt', 'Tricepsbrachii laterale kop bewerkt')
%de offset van de signalen is door het filteren verwijdert. 

%Bereken de gelijkgerichte waarde voor elke spier
% formule = mean(abs(Y))

jSt = [0 0];
hGem = [0 0];

for i = 1:2
    if i == 1
        MVC_filt = Bicepbrachii_filt;
    elseif i == 2
        MVC_filt = Brachioradialis_filt; 
    end
    
    for j = 1:N1-999
        gem = mean(abs(MVC_filt(j:j+999)));
        
        if gem>hGem(i)
            jSt(i)=j;
            hGem(i)=gem;
        end
    end
    
    Begintijd(i) = t1(jSt(i));
    Eindtijd(i) = t1(jSt(i)+999);
end

jSt2 = [0 0];
hGem2 = [0 0];

for i2 = 1:2
    if i2 == 1 
        MVC_filt2 = Tricepsbrachii_lange_kop_filt;
    elseif i2 == 2
        MVC_filt2 = Tricepsbrachii_laterale_kop_filt;
    end
    
    for j2 = 1:N2-999
        gem2 = mean(abs(MVC_filt2(j2:j2+999)));
        
        if gem2>hGem2(i2)
            jSt2(i2)=j2;
            hGem2(i2)=gem2;
        end
    end
    
    Begintijd2(i2) = t2(jSt2(i2));
    Eindtijd2(i2) = t2(jSt2(i2)+999);
end

XarvBicepbrachii_filt = hGem(1);
XarvBrachioradialis_filt = hGem(2);
XarvTricepsbrachii_lange_kop_filt = hGem2(1);
XarvTricepsbrachii_laterale_kop_filt = hGem2(2);

%% Opdracht 4: Elektromechanische vertraging
load('opdracht4_0.5HZ_1-19-26Feb2020at11h16m40s.mat')
EMG_05Hz = DataMinute001; 
load('opdracht4_2HZ_1-20-26Feb2020at11h17m54s.mat')
EMG_2Hz = DataMinute001; 

%door een kapotte kabel aangesloten op de tricepsbrachii laterale kop is
%dit signaal niet bruikbaar deze is dan ook elke keer met % ervoor 

EMG_05Hz = double(EMG_05Hz);
EMG_2Hz = double(EMG_2Hz);

F_05Hz_V = EMG_05Hz(:,1);
Bicepbrachii_05 = EMG_05Hz(:,2);
Brachioradialis_05 = EMG_05Hz(:,3);
Tricepsbrachii_laterale_kop_05 = EMG_05Hz(:,5);

F_2Hz_V = EMG_2Hz(:,1);
Bicepbrachii_2 = EMG_2Hz(:,2);
Brachioradialis_2 = EMG_2Hz(:,3);
Tricepsbrachii_laterale_kop_2 = EMG_2Hz(:,5);

%plot onbewerkte signalen tegen de tijd 0.5Hz
%tijdas maken
N3 = length(Bicepbrachii_05);
k3 = [0:N3-1];
t3=dt*k3;

figure(6)
plot(t3, Bicepbrachii_05)
hold on 
plot(t3, Brachioradialis_05)
%hold on 
%plot(t3, Tricepsbrachii_laterale_kop_05)
title('EMG-signaal 0.5Hz')
xlabel('tijd (s)')
ylabel('spanning (mV)')
legend('Bicepbrachii', 'Brachioradialis')

%plot onbewerkte signalen tegen de tijd 2Hz
%tijdas maken
N4 = length(Bicepbrachii_2);
k4 = [0:N4-1];
t4=dt*k4;

figure(7)
plot(t4, Bicepbrachii_2)
hold on 
plot(t4, Brachioradialis_2)
%hold on 
%plot(t4, Tricepsbrachii_laterale_kop_2)
title('EMG-signaal 2Hz')
xlabel('tijd (s)')
ylabel('spanning (mV)')
legend('Bicepbrachii', 'Brachioradialis')

%filter het signaal
Bicepbrachii_05_filt = filtfilt(B,A,Bicepbrachii_05);
Brachioradialis_05_filt = filtfilt(B,A,Brachioradialis_05);
%Tricepsbrachii_laterale_kop_05_filt = filtfilt(B,A,Tricepsbrachii_laterale_kop_05);
Bicepbrachii_2_filt = filtfilt(B,A,Bicepbrachii_2);
Brachioradialis_2_filt = filtfilt(B,A,Brachioradialis_2);
%Tricepsbrachii_laterale_kop_2_filt = filtfilt(B,A,Tricepsbrachii_laterale_kop_2);

%bepaal de genormaliseerde gelijkgreichte EMG-signalen
%formule = abs(x)/Xarv *100

Bicepbrachii_05_filt_norm = (abs(Bicepbrachii_05_filt)/XarvBicepbrachii_filt)*100;
Brachioradialis_05_filt_norm = (abs(Brachioradialis_05_filt)/XarvBrachioradialis_filt)*100;
%Tricepsbrachii_laterale_kop_05_filt_norm =(abs(Tricepsbrachii_laterale_kop_05_filt)/XarvTricepsbrachii_lange_kop_filt)*100;
Bicepbrachii_2_filt_norm = (abs(Bicepbrachii_2_filt)/XarvBicepbrachii_filt)*100;
Brachioradialis_2_filt_norm = (abs(Brachioradialis_2_filt)/XarvBrachioradialis_filt)*100;
%Tricepsbrachii_laterale_kop_2_filt_norm =(abs(Tricepsbrachii_laterale_kop_2_filt)/XarvTricepsbrachii_lange_kop_filt)*100;

%Omhullende genormaliseerde gelijkgerichte EMG-signalen (dit doe je door genormaliseerde signaal te filteren met laag doorlaat filter, 
%we gebruiken een afsnij frequentie van 4Hz die is gegeven)
fc2 = 4;
Wn2 = fc2/(fs/2);
[B2,A2] = butter(N,Wn2);

Bicepbrachii_05_filt_norm_filt = filtfilt(B2,A2,Bicepbrachii_05_filt_norm);
Brachioradialis_05_filt_norm_filt = filtfilt(B2,A2,Brachioradialis_05_filt_norm);
%Tricepsbrachii_laterale_kop_05_filt_norm_filt = filtfilt(B2,A2,Tricepsbrachii_laterale_kop_05_filt_norm);
Bicepbrachii_2_filt_norm_filt = filtfilt(B2,A2,Bicepbrachii_2_filt_norm);
Brachioradialis_2_filt_norm_filt = filtfilt(B2,A2,Brachioradialis_2_filt_norm);
%Tricepsbrachii_laterale_kop_2_filt_norm_filt = filtfilt(B2,A2,Tricepsbrachii_laterale_kop_2_filt_norm);

%Plot het resultaat
figure (8)
plot(t3,Bicepbrachii_05_filt_norm_filt)
hold on
plot(t3,Brachioradialis_05_filt_norm_filt)
%hold on
%plot(t3,Tricepsbrachii_laterale_kop_05_filt_norm_filt)
title('genormaliseerd EMG-signaal 0.5Hz')
xlabel('tijd (s)')
ylabel('Activiteit (%MVC)')
legend('Bicepbrachii', 'Brachioradialis')

figure (9)
plot(t4,Bicepbrachii_2_filt_norm_filt)
hold on
plot(t4,Brachioradialis_2_filt_norm_filt)
%hold on
%plot(t3,Tricepsbrachii_laterale_kop_2_filt_norm_filt)
title('genormaliseerd EMG-signaal 2Hz')
xlabel('tijd (s)')
ylabel('Activiteit (%MVC)')
legend('Bicepbrachii', 'Brachioradialis')

%zet het signaal om in Newton met de kalibratielijn, de kracht is tijdens
%het practicum gemeten in volt en we willen het omzetten naar Newton

F_05_N = (F_05Hz_V - coef(2))/coef(1);
F_2_N = (F_2Hz_V - coef(2))/coef(1);

%plot de kracht in Newton tegen de tijd
figure(10)
plot(t3,F_05_N)
title('Kracht tijdens contractie van 0.5 Hz')
xlabel('tijd (s)')
ylabel('kracht (N)')

figure(11)
plot(t4,F_2_N)
title('Kracht tijdens contractie van 2 Hz')
xlabel('tijd (s)')
ylabel('kracht (N)')

%vergelijken van de kracht met de genormaliseerde waarden

figure (12)
yyaxis right
plot(t3,F_05_N,'g')
ylabel('Kracht (N)')
yyaxis left
plot(t3,Bicepbrachii_05_filt_norm_filt,'r')
hold on 
plot(t3,Brachioradialis_05_filt_norm_filt,'b')
ylabel('Activiteit (%MVC)')
xlabel('tijd (s)')
title('vergelijking kracht met genormaliseerde EMG 0.5Hz')
legend ('Bicepbrachii','Brachioradialis')

figure (13)
yyaxis right
plot(t4,F_2_N)
ylabel('Kracht (N)')
yyaxis left
plot(t4,Bicepbrachii_2_filt_norm_filt,'b')
hold on 
plot(t4,Brachioradialis_2_filt_norm_filt,'r')
ylabel('Activiteit (%MVC)')
xlabel('tijd (s)')
title('vergelijking kracht met genormaliseerde EMG 2Hz')
legend ('Bicepbrachii', 'Brachioradialis')

%De pieken van het krachtsignaal komen net iets later dan de pieken in het
%EMG signaal, dit is het EMD (electromechanical delay). Als we kijken naar
%de pieken rond de 2 seconde van de brachioradialis op 0.5HZ zien we de piek van het EMG signaal op 1.742s
%en de piek van de kracht zit op 1.907 de EMD is dus ongeveer 0.165
%seconden. 

%Bereken de EMD van de 0.5 HZ van de brachioradialis met behulp van de kruiscorrelatie en plot deze
[kruiscorr, lags]= xcov(F_05_N,Brachioradialis_05_filt_norm_filt,'coeff');
tau = lags*dt;
%emg loop voor op kracht

figure(14)
plot(tau,kruiscorr)
title('kruiscorrelatie tegen tau 0.5Hz')
xlabel('tijdsverschuiving (s)')
ylabel('kruiscorrelatie')
xlim([-1 1])

%bepalen op welk punt de kruiscorrelatie het hoogst is, dat punt is imax
%De EMD is dus op het punt imax 
[kruiscorrmax, imax] = max(kruiscorr);
EMD = tau(imax);

%Bereken de EMD van de 2 HZ van de brachioradialis met behulp van de kruiscorrelatie en plot deze
[kruiscorr2, lags2]= xcov(F_2_N,Brachioradialis_2_filt_norm_filt,'coeff');
tau2 = lags2*dt;
figure(15)
plot(tau2,kruiscorr2)
title('kruiscorrelatie tegen tau 2Hz')
xlabel('tijdsverschuiving (s)')
ylabel('kruiscorrelatie')

[kruiscorrmax2, imax2] = max(kruiscorr2);
EMD2 = tau2(imax2);

%ja de contractiefrequentie heeft invloed op de EMD, deze is aanzienlijk
%lager bij een hogere contractiefrequentie. 

%% Opdracht 5: Bepaal het effect van vermoeidheid op het EMG-signaal
load('opdracht5_60-22-26Feb2020at11h25m24s.mat')
EMG_vermoeid = [DataMinute001; DataMinute002];

EMG_vermoeid = double(EMG_vermoeid);
Bicepbrachii_vermoeid = EMG_vermoeid(:,2); 

%plot het krachtsignaal

%tijdas maken
fs = 2000;
dt = 1/fs;
N5 = length(Bicepbrachii_vermoeid);
k5 = 0:N5-1;
t5 = k5*dt; 

figure(16)
plot(t5,Bicepbrachii_vermoeid)
title('EMG bicepbrachii')
xlabel('tijd (s)')
ylabel('spanning (mV)')

%filter het EMG-signaal (ik maak eerst powerspectrum om te kijken waar de ruis zit)
[Pxx,F]= pwelch(Bicepbrachii_vermoeid - mean(Bicepbrachii_vermoeid),[],[],[],fs);
figure (17)
plot(F,Pxx)
title('Powerspectrum')
xlabel('Frequentie (Hz)')
ylabel('PSD (V^2/Hz)')
xlim([0 500])

N=2;
fc3 = [20 500];
Wn3 = fc3/(fs/2);
[B3,A3] = butter(N,Wn3,'bandpass');
Bicepbrachii_vermoeid_filt = filtfilt(B3,A3,Bicepbrachii_vermoeid);

%plot gefilterde signaal
figure(18)
plot(t5,Bicepbrachii_vermoeid_filt)
title('EMG bicepbrachii gefilterd')
xlabel('tijd (s)')
ylabel('spanning (mV)')

%Verdeel het signaal in tijdsintervallen  van 10 seconden
samples_interval = 10/dt; %in 10 sec zitten 20000 samples 
intervallen = length(Bicepbrachii_vermoeid_filt)/samples_interval;

ti1 = Bicepbrachii_vermoeid_filt(1:20000);
ti2 = Bicepbrachii_vermoeid_filt(20001:40000);
ti3 = Bicepbrachii_vermoeid_filt(40001:60000);
ti4 = Bicepbrachii_vermoeid_filt(60001:80000);
ti5 = Bicepbrachii_vermoeid_filt(80001:100000);
ti6 = Bicepbrachii_vermoeid_filt(100001:120000);
ti7 = Bicepbrachii_vermoeid_filt(120001:140000);
ti8 = Bicepbrachii_vermoeid_filt(140001:160000);

%bereken voor elke tijdsinterval de mediaanfrequentie en de gemiddelde
%frequentie. (dit moet kunnen in een loop maar ik weet niet hoe je die moet opstellen, 
%dus ik doe het voor elk tijdsinterval appart)

%formule gemiddelde frequentie = sum(f.*p)/sum(p) (f = frequentie p=power)

%eerst frequentie en power berekenen met pwelch voor elk tijdsinterval
[p1,f1]= pwelch(ti1- mean(ti1),[],[],[],fs);
[p2,f2]= pwelch(ti2- mean(ti2),[],[],[],fs);
[p3,f3]= pwelch(ti3- mean(ti3),[],[],[],fs);
[p4,f4]= pwelch(ti4- mean(ti4),[],[],[],fs);
[p5,f5]= pwelch(ti5- mean(ti5),[],[],[],fs);
[p6,f6]= pwelch(ti6- mean(ti6),[],[],[],fs);
[p7,f7]= pwelch(ti7- mean(ti7),[],[],[],fs);
[p8,f8]= pwelch(ti8- mean(ti8),[],[],[],fs);

%gemiddelde frequentie berekenen voor elk tijdsinterval
gemf1 = sum(f1.*p1)/sum(p1);
gemf2 = sum(f2.*p2)/sum(p2);
gemf3 = sum(f3.*p3)/sum(p3);
gemf4 = sum(f4.*p4)/sum(p4);
gemf5 = sum(f5.*p5)/sum(p5);
gemf6 = sum(f6.*p6)/sum(p6);
gemf7 = sum(f7.*p7)/sum(p7);
gemf8 = sum(f8.*p8)/sum(p8);

%formule mediaanfrequentie = f(sum(cumsum(p)<=sum(p)/2))
%mediaanfrequentie berekenen 
mf1 = f1(sum(cumsum(p1)<=sum(p1)/2));
mf2 = f2(sum(cumsum(p2)<=sum(p2)/2));
mf3 = f3(sum(cumsum(p3)<=sum(p3)/2));
mf4 = f4(sum(cumsum(p4)<=sum(p4)/2));
mf5 = f5(sum(cumsum(p5)<=sum(p5)/2));
mf6 = f6(sum(cumsum(p6)<=sum(p6)/2));
mf7 = f7(sum(cumsum(p7)<=sum(p7)/2));
mf8 = f8(sum(cumsum(p8)<=sum(p8)/2));

%plot de mediaanfrequenties en de gemiddelde frequentie in een staafdiagram

%eerst vector maken van de mediaanfrequenties en de gemiddelde frequentie
gemf = [gemf1 gemf2 gemf3 gemf4 gemf5 gemf6 gemf7 gemf8];
mf = [mf1 mf2 mf3 mf4 mf5 mf6 mf7 mf8];

%maak vector voor de volgnummers van het tijdsinterval
vn = 1:1:8;

%maak staafdiagram van mediaanfrequenties tegen het volgnummer
figure(19)
bar(vn,mf)
title('mediaanfrequenties van de bicepbrachii')
xlabel('tijdsintervallen')
ylabel('mediaanfrequentie (Hz)')

%maak staafdiagram van de gemiddelde frequentie tegen het volgnummer
figure(20)
bar(vn,gemf)
title('gemiddelde frequenties van de bicepbrachii')
xlabel('tijdsintervallen')
ylabel('gemiddelde frequenties (Hz)')

%bereken voor elk tijdsinterval de effectieve waarde en de gemiddelde
%gelijkgerichte waarde
%formule effectieve waarde = wortel van de gemiddelde power in matlab rms(x)
% formule gelijkgerichte waarde = gemiddelde van de absolute waarde in
% matlab mean(abs(x))

% effectieve waarde berekennen
ef1 = rms(ti1);
ef2 = rms(ti2);
ef3 = rms(ti3);
ef4 = rms(ti4);
ef5 = rms(ti5);
ef6 = rms(ti6);
ef7 = rms(ti7);
ef8 = rms(ti8);

ef = [ef1 ef2 ef3 ef4 ef5 ef6 ef7 ef8];

%gemiddelde gelijkgerichte waarde berekennen
ggw1 =  mean(abs(ti1));
ggw2 =  mean(abs(ti2));
ggw3 =  mean(abs(ti3));
ggw4 =  mean(abs(ti4));
ggw5 =  mean(abs(ti5));
ggw6 =  mean(abs(ti6));
ggw7 =  mean(abs(ti7));
ggw8 =  mean(abs(ti8));

ggw = [ggw1 ggw2 ggw3 ggw4 ggw5 ggw6 ggw7 ggw8];

%plot de effectieve waarde en de gemiddelde gelijkgreichte waarde tegen het
%volgnummer van het tijdsinterval
figure(21)
bar(vn,ef)
title('effectieve waarde van het EMG-signaal van de bicepbrachii')
xlabel('tijdsintervallen')
ylabel('effectieve waarde (Hz)')

figure(22)
bar(vn,ggw)
title('gemiddelde gelijkgerichte waarde van het EMG-signaal van de bicepbrachii')
xlabel('tijdsintervallen')
ylabel('gemiddelde gelijkgerichte waarde (Hz)')

%% Opdracht 6: Verschillende koppen van dezelfde spier
load('opdracht6_beide_-24-26Feb2020at11h45m37s.mat')
EMG_extensie_elleboog = [DataMinute001; DataMinute002];
load('opdracht6_langeretro_1-23-26Feb2020at11h41m54s.mat')
EMG_retroflexie = [DataMinute001; DataMinute002];

EMG_extensie_elleboog = double (EMG_extensie_elleboog);
EMG_retroflexie = double(EMG_retroflexie);

EMG_extensie_elleboog_lange_kop = EMG_extensie_elleboog(:,4);
EMG_extensie_elleboog_laterale_kop = EMG_extensie_elleboog(:,5);
EMG_retroflexie_lange_kop = EMG_retroflexie(:,4);
EMG_retroflexie_laterale_kop = EMG_retroflexie(:,5);

%Filter het signaal 
N=2;
fc4 = [20 500];
Wn4 = fc4/(fs/2);
[B4,A4] = butter(N,Wn4,'bandpass');
EMG_extensie_elleboog_lange_kop_filt = filtfilt(B4,A4,EMG_extensie_elleboog_lange_kop);
EMG_extensie_elleboog_laterale_kop_filt = filtfilt(B4,A4,EMG_extensie_elleboog_laterale_kop);
EMG_retroflexie_lange_kop_filt = filtfilt(B4,A4,EMG_retroflexie_lange_kop);
EMG_retroflexie_laterale_kop_filt = filtfilt(B4,A4,EMG_retroflexie_laterale_kop);

%normaliseer het EMG-signaal naar %MVC 
%formule: |x|/Xarv *100

%maak het genormaliseerde signaal
EMG_extensie_elleboog_lange_kop_norm = (abs(EMG_extensie_elleboog_lange_kop_filt)/XarvTricepsbrachii_lange_kop_filt)*100;
EMG_extensie_elleboog_laterale_kop_norm = (abs(EMG_extensie_elleboog_laterale_kop_filt)/XarvTricepsbrachii_laterale_kop_filt)*100;
EMG_retroflexie_lange_kop_norm = (abs(EMG_retroflexie_lange_kop_filt)/XarvTricepsbrachii_lange_kop_filt)*100;
EMG_retroflexie_laterale_kop_norm = (abs(EMG_retroflexie_laterale_kop_filt)/XarvTricepsbrachii_laterale_kop_filt)*100;

%plot de genormaliseerde signalen

%tijd as maken
N6 = length(EMG_extensie_elleboog_lange_kop_norm);
k6 = [0:N6-1];
t6 = k6*dt;

figure(23)
plot(t6,EMG_extensie_elleboog_lange_kop_norm)
%hold on 
%plot (t6,EMG_extensie_elleboog_laterale_kop_norm)
title('EMG van het genormaliseerde signaal tijdens extensie van de elleboog')
xlabel('tijd (s)')
ylabel('EMG (%MVC)')
legend('tricep lange kop')

N7= length(EMG_retroflexie_lange_kop_norm);
k7 = [0:N7-1];
t7 = k7*dt; 

figure(24)
plot(t7, EMG_retroflexie_lange_kop_norm)
% hold on 
% plot(t7,EMG_retroflexie_laterale_kop_norm)
title('EMG van het genormaliseerde signaal tijdens retroflexie')
xlabel('tijd (s)')
ylabel('EMG (%MVC)')
legend('tricep lange kop')

%Omhullende bepalen 
Omhullende_extensie_elleboog_lange_kop = filtfilt(B2,A2,EMG_extensie_elleboog_lange_kop_norm);
Omhullende_extensie_elleboog_laterale_kop_norm = filtfilt(B2,A2,EMG_extensie_elleboog_laterale_kop_norm);
Omhullende_retroflexie_lange_kop_norm = filtfilt(B2,A2,EMG_retroflexie_lange_kop_norm);
Omhullende_retroflexie_laterale_kop_norm = filtfilt(B2,A2,EMG_retroflexie_laterale_kop_norm);

%plot de omhullende
figure(25)
plot(t6,Omhullende_extensie_elleboog_lange_kop)
%hold on
%plot(t6,Omhullende_extensie_elleboog_laterale_kop_norm)
title('Omhullende tijdens extensie van de elleboog')
xlabel('tijd (s)')
ylabel('EMG (%MVC)')
legend('tricep lange kop')

figure(26)
plot(t7,Omhullende_retroflexie_lange_kop_norm)
%hold on
%plot(t7,Omhullende_retroflexie_laterale_kop_norm)
title('Omhullende tijdens retroflexie')
xlabel('tijd (s)')
ylabel('EMG (%MVC)')
legend('tricep lange kop')

%We zien dat de opdracht niet echt goed is gelukt. Wij hadden last van een
%kapotte kabel, waardoor wij geen goede gegevens hebben gekregen. De lange
%kop geeft namelijk geen logische signaal weer. 

%De activiteit van verschillende koppen van dezelfde spier is niet altijd
%gelijk. 
