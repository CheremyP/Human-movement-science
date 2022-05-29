%% Opdracht 1.a

% manueel data  geimporteerd via home en import data 
% en opslaan als mat bestand
% Bestande laden
clc,clear;

load('Groep9ABalstuit.mat')
load('Groep9ABalstuitREF.mat')

% matrix selecteren van de 
Datastuit = Groep9ABalstuit(:,2:3);
DatabalstuitREF = Groep9ABalstuitREF(:,2:3);


% markers Real world cordinaten stuit
Marker1_stuit = DatabalstuitREF(5,:);
Marker2_stuit = DatabalstuitREF(33,:);
Marker3_stuit = DatabalstuitREF(61,:);
Marker4_stuit = DatabalstuitREF(89,:);

%omreken factor stuit
y_comp_stuit = Marker1_stuit - Marker3_stuit;
x_comp_stuit = Marker2_stuit- Marker4_stuit;
lengte_real_world = 0.995;

real_world_stuit = y_comp_stuit(2) / lengte_real_world;

%omrekenen 
Real_world_singaal_suit =  Datastuit./real_world_stuit;

% offset assenstelsel beginen in bal punt
singaal_stuit = Real_world_singaal_suit - Real_world_singaal_suit(:,1);

%tijdsas aanmaken
fs = 50;
N_stuit = length(singaal_stuit);
k = 0:N_stuit-1;
dt = 1/fs;
t = k*dt;

%afgeleiden naar snelheid en vervolgens cversnelling
v_stuit = gradient(singaal_stuit(:,2),dt);


a_stuit = gradient(v_stuit,dt);



% positie snelheid en versneling waneer de bal in de lucht is werkt alleen
% zwaartekracht luchtweerstand is verwaarloosbaar. dus daar kan je val
% versnelling zien

figure(1)
subplot(3,1,1)
plot(t,singaal_stuit)
xlabel('tijd(s)')
ylabel('lengte(m)')

subplot(3,1,2)
plot(t,v_stuit)
xlabel('tijd(s)')
ylabel('snelheid(m/s)')

subplot(3,1,3)
plot(t,a_stuit)
xlabel('tijd(s)')
ylabel('versnelling(m/s^2)')


a_gem = mean(a_stuit(32:69));

close all


%% Opdracht 1 b
% manueel data  geimporteerd via home en import data 
% en opslaan als mat bestand
% Bestande laden

load('Groep9AGooi.mat')
load('Groep9AGooiREF.mat')

% matrix selecteren van de 
Datagooi = Groep9AGooi(:,2:3);
DatagooiREF = Groep9AGooiREF1(:,2:3);

% markers Real world cordinaten gooi
Marker1_gooi = Groep9AGooiREF1(1,:);
Marker2_gooi = Groep9AGooiREF1(5,:);
Marker3_gooi = Groep9AGooiREF1(9,:);
Marker4_gooi = Groep9AGooiREF1(13,:);

%omreken factor gooi
y_comp_gooi = Marker1_gooi - Marker3_gooi;
x_comp_gooi = Marker2_gooi - Marker4_gooi;

real_world_gooi = y_comp_gooi(3) / lengte_real_world;

%omrekenen 
Real_world_singaal_gooi = Datagooi./real_world_gooi;

% offset assenstelsel beginen in bal punt
singaal_gooi = Real_world_singaal_gooi - Real_world_singaal_gooi(1,1:2);




%tijdsas aanmaken
fs = 500;
N_gooi = length(singaal_gooi());
k = 0:N_gooi-1;
dt = 1/fs;
t = k*dt;


v_gooi = gradient(singaal_gooi(1:36,1),dt);
v_gooi_km_h = v_gooi*3.6;


figure(1)
subplot(2,1,1)
plot(t,singaal_gooi(:,1))
xlabel('tijd(s)')
ylabel('lengte(m)')
subplot(3,1,2)
plot(t(1:36),v_gooi_km_h)
xlabel('tijd(s)')
ylabel('snelheid(m/s)')

%gem snelheid van de gooi in horizontale richting
v_gem_gooi = mean(v_gooi_km_h);

%% opdracht 2a
clc;
clear;
[x,y,z, fs] = readndf('TN000002.ndf');
figure(1)
plot(-x,-y)
axis equal

% tijd-as
dt = 1/fs;
N = length(x);
k = [0:N-1];
t = k*dt;

%X coordinaten markers
heup_x   = -x(:,1);
knie_x   = -x(:,2);
enkel_x  = -x(:,3);
pedaal_x = -x(:,4);

figure(2)
plot(t,heup_x)
hold on, plot(t,knie_x)
hold on, plot(t,enkel_x)
hold on, plot(t,pedaal_x)
xlabel 'tijd(s)'
ylabel 'x-coordinaten'
legend ('heup_x','knie_x','enkel_x','pedaal_x')
title 'Fietsexperiment'

%y coordinaten markers
heup_y   = -y(:,1);
knie_y   = -y(:,2);
enkel_y  = -y(:,3);
pedaal_y = -y(:,4);

%filteren
N= 2; %orde
fc = 10; %afsnijfrequentie
Wn = fc/(fs/2);
[B,A] = butter(N,Wn);

% voor x twee maal in tegengestelde richting
gef_heup_x = filtfilt(B,A,heup_x);
gef_knie_x = filtfilt(B,A,knie_x);
gef_enkel_x= filtfilt(B,A,enkel_x);
gef_pedaal_x = filtfilt(B,A,pedaal_x);
figure(3)
subplot(2,1,1)
plot(t,heup_x,t,knie_x,t,enkel_x,t,pedaal_x)
xlabel 'tijd(s)'
ylabel 'x-coordinaten'
legend ('heup_x','knie_x','enkel_x','pedaal_x')
title 'Fietsexperiment'
hold on
subplot(2,1,2)
plot(t,gef_heup_x,t,gef_knie_x,t,gef_enkel_x,t,gef_pedaal_x)
xlabel 'tijd(s)'
ylabel 'x-coordinaten'
legend ('gef heup x','gef knie x','gef enkel x','gef pedaal x')
title 'Fietsexperiment gefilterd singaal'

% voor y twee maal in tegengestelde richting
gef_heup_y = filtfilt(B,A,heup_y);
gef_knie_y = filtfilt(B,A,knie_y);
gef_enkel_y = filtfilt(B,A,enkel_y);
gef_pedaal_y = filtfilt(B,A,pedaal_y);

figure(4)
subplot(2,1,1)
plot(t,gef_heup_y,t,gef_knie_y,t,gef_enkel_y,t,gef_pedaal_y)
xlabel 'tijd(s)'
ylabel 'y-coordinaten'
legend ('heup_y','knie_y','enkel_y','pedaal_y')
title 'Fietsexperiment'
hold on
subplot(2,1,2)
plot(t,gef_heup_y,t,gef_knie_y,t,gef_enkel_y,t,gef_pedaal_y)
xlabel 'tijd(s)'
ylabel 'y-coordinaten'
legend ('gef heup y','gef knie y','gef enkel y','gef pedaal y')
title 'Fietsexperiment gefilterd singaal'

%weining zin om te filteren je ziet niet veel verschil alleen als je
%inzoomt

%opdracht deel b
r_bb =  [gef_heup_x gef_heup_y] - [gef_knie_x gef_knie_y];
l_bb = sqrt(r_bb(:,1).^2 +r_bb(:,2).^2)./1000;
gem_l = mean(l_bb);
st_l = std(l_bb);
figure(5)
plot(t,l_bb)
xlabel 'tijd(s)'
ylabel 'lengte(m)'
title 'bovenbeenlengte'

%de oorzaak van de variatie klan het verschuiven van de makers op het
%broekje zijn
%opdracht deel c
max_pedaal_y = max(pedaal_y);
min_pedaal_y  = min(pedaal_y);
max_pedaal_x = max(pedaal_x);
min_pedaal_x  = min(pedaal_x);

%trap as min -max + begin vector
trapas_x = (max_pedaal_x - min_pedaal_x)/2 +min_pedaal_x ;
trapas_y = (max_pedaal_y - min_pedaal_y)/2  +min_pedaal_y;

%vector crank
r1 = gef_pedaal_x - trapas_x;
r2  = gef_pedaal_y - trapas_y;
crank_hoek = atan2(r2,r1);
crank_hoek_fil = unwrap(crank_hoek);
hoek_v = gradient(crank_hoek_fil,1/fs);

figure(99)
plot(crank_hoek_fil)

figure(6)
plot(t,hoek_v)
xlabel 'tijd(s)'
ylabel 'hoek(rad/-1s)'
title 'hoeksnelheid trapas t.o.v. pedaal '

rpm = mean(hoek_v)*60/(2*pi);
verhschil_rpm = rpm-60;

% opdracht 2 d

%Coordinaten
r_bb_x = heup_x  - knie_x;  
r_bb_y = heup_y  - knie_y;

r_ob_y = knie_y - enkel_y;
r_ob_x = knie_x - enkel_x;

% knie hoek ten opzichte van y as 0 graden
phi1 = 180 + atan2d (r_bb_x,r_bb_y);
phi2 = 180 - atan2d (r_bb_x,r_bb_y);

phi = phi2-phi1;
omega_org = gradient(phi,fs);

gef_r_bb_x = gef_heup_x  - gef_knie_x;  
gef_r_bb_y = gef_heup_y  - gef_knie_y;

gef_r_ob_y = gef_knie_y - gef_enkel_y;
gef_r_ob_x = gef_knie_x - gef_enkel_x;

% knie hoek ten opzichte van y as 0 graden
gef_phi1 = 180 + atan2d (gef_r_bb_x,gef_r_bb_y);
gef_phi2 = 180 - atan2d (gef_r_bb_x,gef_r_bb_y);

gef_phi = gef_phi2-gef_phi1;
gef_omega = gradient(gef_phi,fs);


figure(7)
plot(t,gef_phi)
hold on, plot(t,phi) 
xlabel ('tijd(s)')
ylabel 'hoek(graden)'
legend('gefilterd','oringeel')
title 'kniehoek'


figure(8)
plot(t,omega_org)
hold on, plot(t,gef_omega) 
xlabel ('tijd(s)')
ylabel 'Hoeksnelheid(graden/s^-1'
legend('orgineel','gefilterd')
title 'hoeksnelheid knie'

crank_hoek = atan2d(r2,r1);

figure(9)
plot(gef_omega,crank_hoek)
xlabel ('crankhoek(graden)')
ylabel 'Hoeksnelheid(graden/s^-1'
title 'hoeksnelheid knie t.o.v. crankhoek'

%% opdracht 3 a

[x_9,y_9,z_9, fs] = readndf('TN000009.ndf');
figure(10)
plot(-x_9,-y_9)
axis equal

% tijd-as
dt = 1/fs;
N = length(x_9);
k = [0:N-1];
t = k*dt;

%X coordinaten markers
schouder_x_9 = -x_9(:,4);
heup_x_9   = -x_9(:,1);
knie_x_9   = -x_9(:,3);
enkel_x_9  = -x_9(:,2);


figure(2)
plot(t,schouder_x_9,t,heup_x_9,t,knie_x_9,t,enkel_x_9)
xlabel 'tijd(s)'
ylabel 'x-coordinaten'
legend ('heup_x','knie_x','enkel_x','pedaal_x')
title 'Sprongexperiment'

%y coordinaten markers
schouder_y_9 = -y_9(:,4);
heup_y_9   = -y_9(:,1);
knie_y_9   = -y_9(:,3);
enkel_y_9  = -y_9(:,2);


figure(12)
plot(t,schouder_y_9,t,heup_y_9,t,knie_y_9,t,enkel_y_9)
xlabel 'tijd(s)'
ylabel 'y-coordinaten'
legend ('schouder','heup_x','knie_x','enkel_x')
title 'Sprongexperiment'
close all

schouder =  [schouder_y_9];
heup = [heup_y_9];
knie = [knie_y_9];
enkel = [enkel_y_9];

r_ob = knie - enkel;
r_bb = heup - knie;
r_romp = schouder - heup;

r_ob = r_ob*0.5/1000;
r_bb = r_bb*0.5/1000;
r_romp = r_romp*0.5/1000;


figure(13)
plot(t,r_ob,t,heup_y_9,t,r_bb,t,r_romp)
xlabel 'tijd(s)'
ylabel 'y-coordinaten'
legend ('schouder','heup','knie','enkel')
title 'Sprongexperiment'




