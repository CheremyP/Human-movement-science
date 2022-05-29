function [gecorigeerdsignaal] = verwijderdrift(signaal, onbelast1start, onbelast1eind, onbelast2start, onbelast2eind)
%%
% Begin van onbelast en einde van onbelast (plateus aflezen in grafiek met de bijbehorende begin
% en  eind waardes.
x1 = [onbelast1start:onbelast1eind];
x2 = [onbelast2start:onbelast2eind];

% Corresponderende y-waarde met x ' om de matrix om te zetten. 
y = signaal([x1 x2])';

% Lengte van het signaal zodat je x van ineare formule a*x+b krijgt.
N = length(signaal);
k = [0:N-1];
% polyfit zodat je waarde a en b van lineare formule a*x+b krijgt.
c = polyfit(([x1 x2]) ,y,1);
drift = c(1)*k + c(2);
gecorigeerdsignaal = signaal - drift';