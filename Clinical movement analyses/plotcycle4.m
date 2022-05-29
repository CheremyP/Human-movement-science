function plotcycle4(MODEL)

%______________
% funtion file to visualize the content of the MODEL structured array
% plotcycle4(MODEL)
% input = the structured array (MODEL or other when you renamed this
%           variable)
% output = no output variable exists for this function, the figure is the
%           only output you'll get
% this version contains explicit Y-axes labeling
% authors: Jaap Harlaar & Eline Flux & Han Houdijk
% edited: 29-05-2015
% create figure
   h=figure('Name',['GaitCyclePlot MODEL  ', MODEL.FileName]);
   set(h,'DefaultAxesFontSize',8,'DefaulttextInterpreter','none');
%===============

TimeBase=[0:1:100]; % TimeBase = NormalizedTime 0-100%

h=subplot(7,3,1);
plot(TimeBase,MODEL.JOINT.AngleLowBack)
set(h,'Color',[.8 1 1]); % cyan canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<bw [deg] fw>')
title('AngleLowBack')

h=subplot(7,3,4);
plot(TimeBase,MODEL.JOINT.AngleHipRight)
set(h,'Color',[.8 1 1]); % cyan canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<retro [deg] ante>')
title('AngleHipRight')

h=subplot(7,3,7);
plot(TimeBase,MODEL.JOINT.AngleKneeRight)
set(h,'Color',[.8 1 1]); % cyan canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<ext [deg] flx>')
title('AngleKneeRight')



h=subplot(7,3,10);
plot(TimeBase,MODEL.JOINT.AngleAnkleRight)
set(h,'Color',[.8 1 1]); % cyan canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<plflx [deg] dflx>')
title('AngleAnkleRight')


h=subplot(7,3,13);
plot(TimeBase,MODEL.JOINT.AngleHipLeft)
set(h,'Color',[.5 1 1]); % cyan canvas
hold on
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<retro [deg] ante>')
title('AngleHipLeft')

h=subplot(7,3,16);
plot(TimeBase,MODEL.JOINT.AngleKneeLeft)
set(h,'Color',[.5 1 1]); % cyan canvas
hold on
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<flx [deg] ext>')
title('AngleKneeLeft')

h=subplot(7,3,19);
plot(TimeBase,MODEL.JOINT.AngleAnkleLeft)
set(h,'Color',[.5 1 1]); % cyan canvas
hold on
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<plflx [deg] dflx>')
title('AngleAnkleLeft')


h=subplot(7,3,2);
plot(TimeBase,MODEL.JOINT.SaggMomentHipRight)
set(h,'Color',[1 .9 .9]); % rose canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
ylabel('<ante [Nm] retro>')
title('SaggMomentHipRight')

h=subplot(7,3,5);
plot(TimeBase,MODEL.JOINT.SaggMomentKneeRight)
set(h,'Color',[1 .9 .9]); % rose canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<flx [Nm] ext>')
title('SaggMomentKneeRight')


h=subplot(7,3,8);
plot(TimeBase,MODEL.JOINT.SaggMomentAnkleRight)
set(h,'Color',[1 .9 .9]); % rose canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<dflx [Nm] plflx>')
title('SaggMomentAnkleRight')

h=subplot(7,3,11);
plot(TimeBase,MODEL.JOINT.FrontMomentHipRight)
set(h,'Color',[1 .8 .8]); % rose canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [Axranges(3) Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
l=line([Axranges(1) Axranges(2)],[0 0]); % zero level
set(l,'Color',[0 0 0]); % black
ylabel('<add [Nm] abd>')
title('FrontMomentHipRight')


h=subplot(7,3,3);
plot(TimeBase,MODEL.MUSCLE.EmgGMX)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgGMX')

h=subplot(7,3,6);
plot(TimeBase,MODEL.MUSCLE.EmgGMD)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgGMD')

h=subplot(7,3,9);
plot(TimeBase,MODEL.MUSCLE.EmgRF)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgRF')

h=subplot(7,3,12);
plot(TimeBase,MODEL.MUSCLE.EmgVL)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgVL')

h=subplot(7,3,15);
plot(TimeBase,MODEL.MUSCLE.EmgST)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgST')

h=subplot(7,3,18);
plot(TimeBase,MODEL.MUSCLE.EmgTA)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgTA')

h=subplot(7,3,21);
plot(TimeBase,MODEL.MUSCLE.EmgGAM)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgGAM')

h=subplot(7,3,20);
plot(TimeBase,MODEL.MUSCLE.EmgSOL)
set(h,'Color',[1 1 .8]); % yellow canvas
hold on
Axranges=axis; % Axranges=[xmin xmax ymin ymax]
axis([Axranges(1) Axranges(2) 0 Axranges(4)]); % force zero
l=line([MODEL.EVENT.RightTO MODEL.EVENT.RightTO ], [0 Axranges(4)]);
set(l,'Color',[.5 .5 .5]); % grey
ylabel('uV')
title('EmgSOL')

end


