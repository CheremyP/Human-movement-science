function [angle] = angle2d (vectors);
%                                                Idsart Kingma
%   [angle] = angle2d (vectors);
%
%   berekent de hoek (in radialen) met de rechter horizontaal van
%   een kolom met 2-dimensionale (rij)vektoren. Er wordt gecorrigeerd
%   voor overschrijdingen de grens tussen + en - 180 graden tussen samples
%
%
%  INPUT  [Y1  Z1          OUTPUT :   [phi1
%          ........                     ..          met  m = aantal samples
%          Y1m Z1m]                    phim]
%

angle = atan2 (vectors(:,2), vectors(:,1));  % LET OP: functie atan2 wil de
                                             %  volgorde Y .. X (of Z .. Y)

angle = fixjumps(angle);





