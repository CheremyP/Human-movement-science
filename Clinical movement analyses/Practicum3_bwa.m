close all

TimeBase=[0:1:100]; % TimeBase = NormalizedTime 0-100%

%% enkelhoek

figure (1)
plot(TimeBase,NORM.JOINT.AngleAnkleRight, 'b' ,TimeBase,TC2.JOINT.AngleAnkleRight, 'r')
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([NORM.EVENT.RightTO NORM.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color','b'); % blauw

Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([TC2.EVENT.RightTO TC2.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color','r'); % rood

l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color', [ 0 0 0 ]); % black

xlabel('gait cycle [%]')
ylabel('<plflx [deg] dflx>')
title('AngleAnkleRight')
legend ('normale enkel ' , 'gekke gang')

%% kniehoek

figure (2)
plot(TimeBase,NORM.JOINT.AngleKneeRight, 'b' ,TimeBase,TC2.JOINT.AngleKneeRight, 'r')
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([NORM.EVENT.LeftIC NORM.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','b'); % blauw

Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([-TC2.EVENT.LeftIC -TC2.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','r'); % rood

l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color', [ 0 0 0 ]); % black

xlabel('gait cycle [%]')
ylabel('<ext [deg] flx>')
title('AngleKneeRight')

legend ('normale gang ' , 'gekke gang')

%% heuphoek

figure (3)
plot(TimeBase,NORM.JOINT.AngleHipRight, 'b' ,TimeBase,TC2.JOINT.AngleHipRight, 'r')
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([NORM.EVENT.LeftIC NORM.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','b'); % blauw

Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([-TC2.EVENT.LeftIC -TC2.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','r'); % rood

l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color', [ 0 0 0 ]); % black

xlabel('gait cycle [%]')
ylabel('<retro [deg] ante>')
title('AngleHipRight')
legend ('normale gang ' , 'gekke gang')

%% moment voet


figure (4)
plot(TimeBase,NORM.JOINT.SaggMomentAnkleRight , 'b ',  TimeBase , TC2.JOINT.SaggMomentAnkleRight , 'r ')
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([NORM.EVENT.LeftIC NORM.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','b'); % blauw

Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([-TC2.EVENT.LeftIC -TC2.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','r'); % rood

l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color', [ 0 0 0 ]); % black

ylabel('<dflx [Nm] plflx>')
xlabel('gait cycle [%]')

title('SaggMomentAnkleRight')

legend ('normale gang ' , 'gekke gang')


%% moment  knie

figure (5)

plot(TimeBase,NORM.JOINT.SaggMomentKneeRight, 'b' ,TimeBase,TC2.JOINT.SaggMomentKneeRight, 'r')
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([NORM.EVENT.LeftIC NORM.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','b'); % blauw

Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([-TC2.EVENT.LeftIC -TC2.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','r'); % rood

l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color', [ 0 0 0 ]); % black

ylabel('<ante [Nm] retro>')
xlabel('gait cycle [%]')

title('SaggMomentKneeRight')
legend ('normale gang ' , 'gekke gang')


%% moment heup


figure (6)

plot(TimeBase,NORM.JOINT.SaggMomentHipRight, 'b' ,TimeBase,TC2.JOINT.SaggMomentHipRight, 'r')
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([NORM.EVENT.LeftIC NORM.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','b'); % blauw

Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([-TC2.EVENT.LeftIC -TC2.EVENT.LeftIC ], [Axranges(3) Axranges(4)]);
set(l,'Color','r'); % rood

l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color', [ 0 0 0 ]); % black
ylabel('<flx [Nm] ext>')
xlabel('gait cycle [%]')

title('SaggMomentHipRight')
legend ('normale gang ' , 'gekke gang')