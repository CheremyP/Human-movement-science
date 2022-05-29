% Opdracht 1

load('MyData-01-10Mar2020at14h13m03s.mat')
Data_1 = double(DataMinute001); 
load('MyData-02-10Mar2020at14h16m56s.mat')% kalibratie_2 waarde
Data_2 = double(DataMinute001); 
load('MyData-03-10Mar2020at14h18m02s.mat') % kalibratie_3 waarde
Data_3 = double(DataMinute001);
load('MyData-04-10Mar2020at14h19m02s.mat')% kalibratie_4 waarde
Data_4 = double(DataMinute001);
load('MyData-05-10Mar2020at14h20m57s.mat')% kalibratie_5 waarde
Data_5 = double(DataMinute001);

% inputs 
% 1 = kracht
% 2 = biceps
% 3 = brachioradialis
% 4 = triceps lang
% 5 = triceps kort
% 6 = referentie

gem_1 = mean(Data_1(127:20850,1));
gem_2 = mean(Data_2(200:21130,1));
gem_3 = mean(Data_3(400:21000,1));
gem_4 = mean(Data_3(200:20500,1));
gem_5 = mean(Data_3(200:20700,1));

singaal_gem = [gem_1 gem_2 gem_3 gem_4 gem_5];

x = [0.0343 4.879+0.0344 4.882+4.879+0.0343  4.882+4.879+0.0343+4.828 4.882+4.879+0.0343+4.828+4.858];
werkelijke_F = x * 9.81;
cal_lijn = polyfit(werkelijke_F,singaal_gem,1);

x = [0:200];
y = cal_lijn(1,1)*x + cal_lijn(1,2);

figure(1)
plot(x,y)
hold on, plot(werkelijke_F(1,1), gem_1 , 'o',werkelijke_F(1,2), gem_2 , 'o', werkelijke_F(1,3), gem_3 , 'o', werkelijke_F(1,4), gem_4 , 'o' , werkelijke_F(1,5), gem_5 , 'o')
xlabel 'kracht (N)'
ylabel 'Spanning (mV)'
title('Kalibratielijn bij bekende krachten')

% Ja er is sprake van offset, dat zie je aan b en aan de grafiek die niet
% bij nul begint

F_cal = (y- cal_lijn(1,2))/cal_lijn(1,1);

% Opdracht 3 

load('MyData-06-10Mar2020at14h41m27s.mat')
Data_6 = DataMinute001;
Data_6 = double(Data_6);

N = length(Data_6(:,1));
fs = 2000;
dt = 1/fs;
k = [0:N-1];
t = k*dt;

figure(2)
plot(t,Data_6(:,1), t,Data_6(:,2), t, Data_6(:,3), t, Data_6(:,4), t, Data_6(:,5))
xlabel('tijd(s)')
ylabel 'spanning (mV)'
legend ('Kracht', 'Biceps', 'Brachioradialis', 'Triceps lang', 'Triceps kort')

% filteren van signaal
fs = 2000;
fc = 20;
N = 2;
Wn = fc/(fs/2);
[B,A] = butter(N,Wn, 'high');

gef_data_6 = filtfilt(B,A,Data_6(:,1:5));

figure(3)
plot(t,gef_data_6(:,2), t, Data_6(:,2))
xlabel('tijd(s)')
ylabel 'spanning (mV)'
legend ('gefilterd signaal biceps', 'orgineel signaal biceps')

figure(4)
plot(t,gef_data_6(:,3), t, Data_6(:,3))
xlabel('tijd(s)')
ylabel 'spanning (mV)'
legend ('gefilterd signaal brachioradialis', 'orgineel signaal brachioradialis')

figure(5)
plot(t,gef_data_6(:,4), t, Data_6(:,4))
xlabel('tijd(s)')
ylabel 'spanning (mV)'
legend ('gefilterd signaal triceps lang', 'orgineel signaal triceps lang')

figure(6)
plot(t,gef_data_6(:,5), t, Data_6(:,5))
xlabel('tijd(s)')
ylabel 'spanning (mV)'
legend ('gefilterd signaal triceps kort', 'orgineel signaal triceps kort')

% Xarv bereken tijdens Mvc
gef_signaal_biceps = gef_data_6(:,2);
N = length(gef_data_6);
fs = 2000;
N_interval = N/(fs/2); % interval van 2 seconde

x_arv_biceps = [0];
jSt = [0];
for i = 1:3
    for j = 1:N-499 % 500ms komt overeen met 500 samples
        
        x_arv_6_biceps = mean(abs(gef_signaal_biceps(j:j+499,i)));
       
        if x_arv_6_biceps>x_arv_biceps(i)
            jSt(i)= j;
            x_arv_biceps(i) = x_arv_biceps;
        end
        Begintijd(i) = t1(jSt(i));
        Eindtijd(i) = t1(jSt(i)+499);
    end
   
end

gef_signaal_brachioradialis = gef_data_6(:,3);
N = length(gef_data_6);
fs = 2000;
N_interval = N/(fs/2); % interval van 2 seconde

x_arv_brachioradialis = [0];
jSt = [0];
for i = 1:3
    for j = 1:N-499 % 500ms komt overeen met 500 samples
        
        x_arv_6_brachioradialis = mean(abs(gef_signaal_brachioradialis(j:j+499,i)));
       
        if x_arv_6_brachioradialis>x_arv_brachioradialis(i)
            jSt(i)= j;
            x_arv_brachioradialis(i) = x_arv_brachioradialis;
        end
        Begintijd(i) = t1(jSt(i));
        Eindtijd(i) = t1(jSt(i)+499);
    end
   
end


