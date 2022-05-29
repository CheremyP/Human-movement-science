function [y, z] = dattoyz(data);
%
%                                                  Idsart Kingma
%
% omzetten van datastructuur:
%
%
% van    [Y1  Z1  .... Yk  Zk             met: k = aantal markers
%         ..................                   m = aantal samples
%         Y1m Z1m .... Ykm Zkm]
%
%
%   naar 2 matrices Y en Z van format    [Y1  Y2  Yk
%                                        ...........
%                                        Y1m Y2m Ykm]



[m,n] = size (data);
y = []; z = [];
for i = 1:n/2;
  y = [y data(:,i*2-1)];
  z = [z data(:,i*2)];
end;


