function [moment] = cross2d(arm,kracht)

moment = arm(1).*kracht(:,2) - arm(2).*kracht(:,1);