function [yloc, zloc]=createaxes (distaal,proximaal);
%
%                                                  Idsart Kingma
%
% Creert, op basis van een distale en proximale marker op een segment,
% een in de tijd varierend lokaal assenstelsel. zloc is de lokale z-as (compressie
% component), met de positieve richting van distaal naar proximaal.
% yloc is de lokale y-as (afschuif-component). Deze as is 90 graden rechtsom geroteerd ten opzichte van de z-as. 

    zloc = (proximaal-distaal);
    lengzloc = sqrt(zloc(:,1).^2 + zloc(:,2).^2);
    zloc=zloc./ [lengzloc lengzloc];
    yloc = [zloc(:,2) -zloc(:,1)];
