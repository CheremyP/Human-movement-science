
%% Opdracht 1.1
% sommatie van alle z-waarde van de krachten
load('TN000002.afp')
singaal_som_z = TN000002(:,5) + TN000002(:,6) + TN000002(:,7) + TN000002(:,8);

% plot van signaal
figure(1)
plot(singaal_som_z)
xlabel('spanning(V)')
ylabel('index(i)')
legend('singaal som z')
title('Niet gecorrigeerd signaal')
%% Opdracht 1.2
% functie die corrgigeerd voor de drift
[gecorigeerdsignaal] = verwijderdrift(singaal_som_z, 1,3236,27526,30484);

% plot van signaal met drift en onder drift
figure(2)
plot(singaal_som_z,'r')
hold on, plot(gecorigeerdsignaal,'b')
plot(singaal_som_z)
xlabel('spanning(V)')
ylabel('index(i)')
legend('Niet-gecorrigeerd signaal', 'Gecorrigeerd signaal ')
title(' orgineel en gecorrigeerd signaal')
% Corrigeren is goed gegaan de y waarde van het einde is dichter naar 0 dan
% bij het niet-gecorrigeerde signaal.
% In gevallen van een exponetieel of kwadratische verband.
%% Opdracht 1.3

% laden van alle expermimenten van verschillende ranges.
load('TN000003.afp')
load('TN000004.afp')
load('TN000005.afp')


% range 1
% sommeren van alle z-waarden tot één singaal van z
singaal_som_z_1 = TN000003(:,5) + TN000003(:,6) + TN000003(:,7) + TN000003(:,8);

% drift verwijderen
cor_singaal1 = verwijderdrift(singaal_som_z_1, 1 ,2965 ,14352,15566);

% Berekening van de werekelijke krachten
massa_1 = [0 4.941 (4.941+4.982) (4.941+4.982+4.858) (4.828+4.941+4.982+4.858) (4.860+4.828+4.941+4.982+4.858)];
F_1 = massa_1.*9.81;


% waarde opstellen voor de polyfit werkelijke krachten tegen over
% corrosponderen signaal. Telkens gemiddelde genomen van de plateus
gem_waarde = [0 mean(cor_singaal1(3677:5236)) mean(cor_singaal1(6000:7400)) mean(cor_singaal1(8177:9646)) mean(cor_singaal1(10459:11414)) mean(cor_singaal1(11881:13313))];
coef_1 = polyfit(F_1, gem_waarde,1);
cal_lijn =  coef_1(1)*F_1 + coef_1(2);

% de eerste is in volt/tijd en de bovenstaande is Newton/ volt
figure(3);
plot(F_1,cal_lijn,'.',F_1,cal_lijn)
title ('kalibratielijn range 1')
xlabel ('x [kracht in N]')
ylabel ('y [computerwaarde]')

% range 2
% Berekening van de krachten
singaal_som_z_2 = TN000004(:,5) + TN000004(:,6) + TN000004(:,7) + TN000004(:,8);
massa_2 = [0 (24.981+4.828) (24.981+4.828 +9.883) (24.981+4.828 +9.883+9.9461) (24.981+4.828 +9.883+9.9461+10.918) (24.981+4.828 +9.883+9.9461+10.918+10.009)];
F_2 = massa_2.*9.81;

% drift verwijderen
cor_singaal2 = verwijderdrift(singaal_som_z_2, 1 ,2304 ,20874,21913);

% waarde opstellen voor de polyfit
gem_waarde2 = [0 mean(cor_singaal2(3260:7159)) mean(cor_singaal2(8279:10968)) mean(cor_singaal2(11860:13860)) mean(cor_singaal2(15330:16890)) mean(cor_singaal2(17772:18780))];
coef_2 = polyfit(F_2, gem_waarde2,1);
cal_lijn2 =  coef_2(1)*F_2 + coef_2(2);

% de eerste is in volt/tijd en de bovenstaande is Newton/ volt
figure(4);
plot(F_2,cal_lijn2,'.',F_2,cal_lijn2)
title ('kalibratielijn range 2')
xlabel ('x [kracht in N]')
ylabel ('y [computerwaarde]')

% range 3
% Berekening van de werkelijke krachten
singaal_som_z_3 = TN000005(:,5) + TN000005(:,6) + TN000005(:,7) + TN000005(:,8);
massa_3 = [0 (24.981+4.828) (24.981+4.828 +9.883) (24.981+4.828 +9.883+9.9461) (24.981+4.828 +9.883+9.9461+10.918) (24.981+4.828 +9.883+9.9461+10.918+10.009)];
F_3 = massa_3.*9.81;

% drift verwijderen
cor_singaal3 = verwijderdrift(singaal_som_z_3, 1 ,1860 ,15000,15600);

% waarde opstellen voor de polyfit werkelijke krachten tegen over
% corrosponderen signaal. Telkens gemiddelde genomen van de plateus
gem_waarde3 = [0 mean(cor_singaal3(2623:4229)) mean(cor_singaal3(5122:7289)) mean(cor_singaal3(7900:9200)) mean(cor_singaal3(10059:11437)) mean(cor_singaal3(12236:13480))];
coef_3 = polyfit(F_3, gem_waarde3,1);
cal_lijn3 =  coef_3(1)*F_3 + coef_3(2);


% de eerste is in volt/tijd en de bovenstaande is Newton/ volt
figure(5);
plot(F_3,cal_lijn3,'.',F_3,cal_lijn3)
title ('kalibratielijn range 3')
xlabel ('x [kracht in N]')
ylabel ('y [computerwaarde]')


%% Opdracht 2 b kalibratielijn


% gevoeligheid 1/a z richting
% a = sensiticity spec (pC/N), b = range met corresponderende pC/10 volt
% a/b = corrosponderen coef van de ranges
A1_z = 3.8/100;
A2_z = 3.8/500;
A3_z = 3.8/1000;
A3_xy = 8/100;

range_1_sen = coef_1(1,1);
range_2_sen = coef_2(1,1);
range_3_sen = coef_3(1,1);



%% Opdracht 3.a bepalen kalibratie lijn
 
load('TN000006.afp')
fs = 10;
N = length(TN000006);
dt = 1/fs;
k  = 0:N-1;
t = k'*dt;

%sommatie van krachten x-,y- en z-krachten
signaal_som_x_6 = TN000006(:,1) + TN000006(:,2);
signaal_som_y_6 = TN000006(:,3) + TN000006(:,4);
signaal_som_z_6 = TN000006(:,5) + TN000006(:,6) + TN000006(:,7) + TN000006(:,8);

% grafieken plotten van N tegenover samplefreq


% krachten berrekend door signaal te delen door de coefficenten 
F_x_6 = signaal_som_x_6/A3_xy;
F_y_6 = signaal_som_y_6/A3_xy;
F_z_6 = signaal_som_z_6/coef_1(1,1) ;

coef_6_x = polyfit(t(141:3000),F_x_6(141:3000),1);
coef_6_y = polyfit(t(141:3000),F_y_6(141:3000),1);
coef_6_z = polyfit(t(141:3000),F_z_6(141:3000),1);

figure(6)
plot(F_z_6)
hold on, plot(F_x_6)
hold on, plot(F_y_6)
xlabel('Kracht (N)')
ylabel('sample frequentie (HZ)')
title('Kracht (N) tegenover samplefrequentie')
legend('f_z','f_x','f_y')
%% Opdracht 3 b


F_z_3 = gem_waarde/coef_1(1,1);
coef_3b = polyfit(F_1,F_z_3,1);
F_z_3_gemeten = coef_3b(1,1)*F_1 + coef_3b(1,2);

% plot N tegen over N
figure(7)
plot(F_1,F_z_3)
hold on, plot(F_1,F_z_3_gemeten)
xlabel('zwaartekracht (N)')
ylabel('z-kracht (N)')

% afwijking van de FSO;
max_afwijking = max(F_1 - F_z_3_gemeten);
fso = 40/coef_3b(1,1);
fso_percentage = max_afwijking/fso*100;

%% Opdracht 3 c

singaal_som_x_1 = TN000003(:,1) + TN000003(:,2);
cor_singaal_X_1 = verwijderdrift(singaal_som_x_1, 1 ,3036 ,14140,15160);
gem_waarde_x_1 = [0 mean(cor_singaal_X_1(3677:5236)) mean(cor_singaal_X_1(6000:7400)) mean(cor_singaal_X_1(8177:9646)) mean(cor_singaal_X_1(10459:11414)) mean(cor_singaal_X_1(11881:13313))];

% Berekening van de krachten doormiddel van delen door de theoretische en
% echte richting coeffiecent bij z-component
F_x_1 = gem_waarde_x_1/A3_xy;
F_z_1 = gem_waarde/A1_z;

figure(8)
plot(F_z_1,F_x_1)
xlabel('z-component(N)')
ylabel('x-component(N)')
overspraak =  F_x_1(2:6)./F_z_1(2:6).*100;
overspraak_gem = mean(overspraak);


%% Opdracht 4 a

%Laden van het sprong signaal
load('TN000008.afp')

% Singaal sommeren en verwijdering van de drift.
singaal_som_z_8 = TN000008(:,5) + TN000008(:,6) + TN000008(:,7) + TN000008(:,8);
cor_singaal_8 = verwijderdrift(singaal_som_z_8, 1 ,95 ,1032,1084);

% Gewicht
M_z_8 = cor_singaal_8 /9.81;
M_z_8 = 1/coef_2(1,1) * M_z_8;

% 1 gewicht meting
gem_cor_singaal_8_1 = mean(M_z_8(128:197));
std_cor_singaal_8_1 = std(M_z_8(128:197));
N_8_1 = length(M_z_8(128:197));
sem_cor_singaal_8_1 = std_cor_singaal_8_1/sqrt(N_8_1);
bw_interval_1 = [gem_cor_singaal_8_1+sem_cor_singaal_8_1*1.96 gem_cor_singaal_8_1-sem_cor_singaal_8_1*1.96];

% 2 gewicht meting
gem_cor_singaal_8_2 = mean(M_z_8(335:395));
std_cor_singaal_8_2 = std(M_z_8(335:395));
N_8_2 = length(M_z_8(335:395));
sem_cor_singaal_8_2 = std_cor_singaal_8_2/sqrt(N_8_2);
bw_interval_2 = [gem_cor_singaal_8_2+sem_cor_singaal_8_2*1.96 gem_cor_singaal_8_2-sem_cor_singaal_8_2*1.96];

% 3 gewicht meting
gem_cor_singaal_8_3 = mean(M_z_8(545:602));
std_cor_singaal_8_3 = std(M_z_8(545:602));
N_8_3 = length(M_z_8(545:602));
sem_cor_singaal_8_3 = std_cor_singaal_8_3/sqrt(N_8_3);
bw_interval_3 = [gem_cor_singaal_8_3+sem_cor_singaal_8_3*1.96 gem_cor_singaal_8_3-sem_cor_singaal_8_3*1.96];



% 4 gewicht meting
gem_cor_singaal_8_4 = mean(M_z_8(745:801));
std_cor_singaal_8_4 = std(M_z_8(745:801));
N_8_4 = length(M_z_8(745:801));
sem_cor_singaal_8_4 = std_cor_singaal_8_4/sqrt(N_8_4);
bw_interval_4 = [gem_cor_singaal_8_4+sem_cor_singaal_8_4*1.96 gem_cor_singaal_8_4-sem_cor_singaal_8_4*1.96];

% 5 gewicht meting
gem_cor_singaal_8_5 = mean(M_z_8(940:992));
std_cor_singaal_8_5 = std(M_z_8(940:992));
N_8_5 = length(M_z_8(940:992));
sem_cor_singaal_8_5 = std_cor_singaal_8_5/sqrt(N_8_5);
bw_interval_5 = [gem_cor_singaal_8_5+sem_cor_singaal_8_5*1.96 gem_cor_singaal_8_5-sem_cor_singaal_8_5*1.96];

plot(M_z_8)
xlabel('massa (kg)')
ylabel('Samplefrequentie')

% Gemiddelde van de in totaal 5 metingen.
massa_cher = mean([gem_cor_singaal_8_1 gem_cor_singaal_8_2 gem_cor_singaal_8_3 gem_cor_singaal_8_4 gem_cor_singaal_8_5]);


%% Opdracht 4 b

%laden van het sprongsingaal
load('TN000010.afp')

signaal_som_x_10 = TN000010(:,1) + TN000010(:,2);
signaal_som_y_10 = TN000010(:,3) + TN000010(:,4);
signaal_som_z_10 = TN000010(:,5) + TN000010(:,6) + TN000010(:,7) + TN000010(:,8);

plot(signaal_som_z_10)

%geen filter want het kan die hoge verandering van pieken niet aan
fs = 200;
g = -9.81;
signaal_som_z_10 = verwijderdrift(signaal_som_z_10, 1 ,1526 ,4678,5716);

% f=m*a a=f/m tweemaal integreren geeft de positie
Fz = g*massa_cher;
F_z_10 = signaal_som_z_10/ coef_3(1,1);

Som_krachten = F_z_10+Fz;
%Som_krachten = max(Som_krachten,0);

% Intergreren van de versnelling -> snelheid ->
% Geen sprake van intergratiedrift omdat je de offset hebt verwijdert
a = (Som_krachten)./massa_cher;
v = cumtrapz(a(2000:4350))./fs;
r = cumtrapz(v)./fs;
spronghoogte = max(r);

%% Opdracht 4 c Bereken Centere of pressure

%laden signaal in range 2
load('TN000011.afp')

% sommaties van alle 
signaal_som_x_11 = TN000011(:,1) + TN000011(:,2);
signaal_som_y_11 = TN000011(:,3) + TN000011(:,4);
signaal_som_z_11 = TN000011(:,5) + TN000011(:,6) + TN000011(:,7) + TN000011(:,8);


signaal_som_x_11_z_drift = verwijderdrift(signaal_som_x_11, 1 ,1526 ,44930,45680);
signaal_som_y_11_z_drift = verwijderdrift(signaal_som_y_11, 1 ,1526 ,44930,45680);
signaal_som_z_11_z_drift = verwijderdrift(signaal_som_z_11, 1 ,1526 ,44930,45680);

signaal_1_z_11_z_drift = verwijderdrift(TN000011(:,5), 1 ,1526 ,44930,45680);
signaal_2_z_11_z_drift = verwijderdrift(TN000011(:,6), 1 ,1526 ,44930,45680);
signaal_3_z_11_z_drift = verwijderdrift(TN000011(:,7), 1 ,1526 ,44930,45680);
signaal_4_z_11_z_drift = verwijderdrift(TN000011(:,8), 1 ,1526 ,44930,45680);

% specificaties van platform
a = 12;
b = 20;
p_z = -6.8;

% Berekenen van de krachten
range_2coef = 8/500;

F_x_11_tot = signaal_som_x_11_z_drift/range_2coef;
F_y_11_tot = signaal_som_y_11_z_drift/range_2coef;
F_z_11_tot = signaal_som_z_11_z_drift/coef_2(1,1);
F_z_11_1 = signaal_1_z_11_z_drift/coef_2(1,1) ;
F_z_11_2 = signaal_2_z_11_z_drift/coef_2(1,1) ;
F_z_11_3 = signaal_3_z_11_z_drift/coef_2(1,1) ;
F_z_11_4 = signaal_4_z_11_z_drift/coef_2(1,1) ;


% formule pagina 186
fs= 200;
N = length(signaal_som_z_11);
dt = 1/fs;
k = 0:N-1;
t = k*dt;

P_x_cm = ((p_z*F_y_11_tot)+b*(F_z_11_1+F_z_11_2 - F_z_11_3 - F_z_11_4))./F_z_11_tot;
P_y_cm = ((p_z*F_x_11_tot)+a*(F_z_11_1 - F_z_11_2 - F_z_11_3  + F_z_11_4))./F_z_11_tot;

figure(9)
plot(t,P_x_cm)
hold on,plot(t,P_y_cm) 
title('coordinaat van het center of pressure')
xlabel('tijd[s]')
ylabel('coordinaat [cm]')
legend('x-coordinaat', 'y-coordinaat')
axis equal

% Power bereken

[Pxx_x,F_X] = pwelch(P_x_cm(3000:22460)-mean(P_x_cm(3000:22460)),[],[],[],200);

figure(10)
plot(F_X,Pxx_x)
xlabel('Frequentie(Hz)')
ylabel('Power')
title('Powerspectrum van Fx open')

[F_X_dicht,Pxx_x_dicht] = pwelch(P_x_cm(26420:43640)-mean(P_x_cm(26420:43640)),[],[],[],200);

figure(11)
plot(F_X_dicht,Pxx_x_dicht)
xlabel('Frequentie(Hz)')
ylabel('Power')
title('Powerspectrum van Fx dicht')

[Pxx_y,F_y] = pwelch(P_y_cm(3288:21720)-mean(P_y_cm(3288:21720)),[],[],[],200);

figure(12)
plot(F_y,Pxx_y)
xlabel('Frequentie(Hz)')
ylabel('Power')
title('Powerspectrum van Fy open')


[Pxx_y_dicht,F_y_dicht] = pwelch(P_y_cm(25900:43720)-mean(P_y_cm(25900:43720)),[],[],[],200);

figure(13)
plot(F_y_dicht,Pxx_y_dicht)
xlabel('Frequentie(Hz)')
ylabel('Power')
title('Powerspectrum van Fy dicht')

close all

figure(14)

plot(P_x_cm(3288:21720), P_y_cm(3288:21720))
axis equal
xlabel('X(cm)')
ylabel('Y(cm)')
title('Ogen open')


figure(15)
plot(P_x_cm(25900:43720), P_y_cm(25900:43720))
axis equal
xlabel('X(cm)')
ylabel('Y(cm)')
title('Ogen dicht')

% Ja, je ziet bij de ogen dicht veel meer schommeling dan bij de andere. 

%% Opdracht 5
close all
fs = 200;
load('TN000012.afp')

signaal_som_x_12 = TN000012(:,1) + TN000012(:,2);
signaal_som_y_12 = TN000012(:,3) + TN000012(:,4);
signaal_som_z_12 = TN000012(:,5) + TN000012(:,6) + TN000012(:,7) + TN000012(:,8);

figure(16)
plot(TN000012)

[Pxx, F] = pwelch(signaal_som_y_12 - mean(signaal_som_y_12),[],[],[],fs);

figure(17)
plot(F,Pxx)
xlabel('Frequentie(Hz)')
ylabel('Power')
title('Omwenteling fiets')
% 1,2Hz per seconde doet hij over 1 omweteling.

figure(18)
plot(TN000012(:,1)) 
hold on , plot(TN000012(:,2))
hold on, plot(signaal_som_x_12)
xlabel('Samplefrequentie')
ylabel('Spanning(V)')
title('Signalen van de x-waardes')
legend('FX 1+2','FX 3+4','FX totaal')

figure(19)
plot(TN000012(:,3)) 
hold on , plot(TN000012(:,4))
hold on, plot(signaal_som_y_12)


xlabel('Samplefrequentie')
ylabel('Spanning(V)')
title('Signalen van de y-waardes ')
legend('Fy 1+2','Fy 3+4','Fy totaal')

%Ja, want er komt een artefact die het signaal beinvloed.
close all
