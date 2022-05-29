function[odata]=afgcol(idata,fs); 
%
% [odata]=afgcol(idata,fs); 
%
% idata = input data matrix met tijdseries in kolommen
% fs = de samplefrequentie in beeldjes per seconden 
% odata = gedifferentieerde output data met tijdseries in kolommen
%
%


% c = gegeven in het programma = differentieer-constanten 
% 
  
% idata = idata'; 
 [m,n]=size(idata); 
  
 c=[-1.5 -4 4 1.5]/14*fs; 
  
 for i=1:n, 
     tussen(:,1)=[0;0;idata(1:m-2,i)]; 
     tussen(:,2)=[0;idata(1:m-1,i)]; 
     tussen(:,3)=[idata(2:m,i);0]; 
     tussen(:,4)=[idata(3:m,i);0;0]; 
  
     tussen=tussen*diag(c); 
     odata(:,i)=sum(tussen')'; 
 end 

% ONDERSTAANDE TER VERMINDERING VAN INSLINGER EFFECTEN BIJ 2e AFGELEIDE     
odata(1,:) = odata (3,:);            
odata(2,:) = odata (3,:);            
odata(m-1,:) = odata (m-2,:);
odata(m,:) = odata (m-2,:);


