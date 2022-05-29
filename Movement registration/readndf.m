function [s1,s2,s3,s4] = readndf(naam);
%function [s1,s2,s3,fr] = readndf(naam);
% read optotrak NDF data files containing reconstructed coordinates.
%
% input:	naam = filename including drive, path, and extension
% output:	s1 = xcoordinates, s2 = ycoordinates s3 = zcoordinates
% 		fr = frequency
% version:	1.0
% system:	Matlab 4.2 (and up) for Windows
% toolbox:	-
% files:	-
% outhor:	Tom Welter, 1998, Frank Zaal, 1997
%
% examples:
%		[x,y,z]=readndf('c:\C#001.ndf'); % read 3D data
%		[x,y]=readndf('c:\C#001.ndf'); % read 2D data
%		[x,y,y,freq]=readndf('c:\C#001.ndf'); % read data, get sample frequency
% hint:
%		plot3d(x,y,z);   % plot markertrajectories
%	  	plot3d(x',y',z') % plot stick figures

[fid,message]=fopen(naam);					% open de file
if fid == -1
	error(message);
end

filetype	=fread(fid,1,'char');			% header info lezen
items		=fread(fid,1,'short');
subitems	=fread(fid,1,'short');
numframes	=fread(fid,1,'long');
frequency	=fread(fid,1,'float32');
usercomment	=fread(fid,60,'char');  
rest		=fread(fid,183,'char');         % rest van header bytes niet gebruikt

aantal  =items*subitems*numframes;

data=fread(fid,aantal,'float32');  			% data lezen
data(data<-10e20) = nan;
data=reshape(data,items*subitems,numframes);
data=data';


if subitems ==1
	s1=data(:,1:1:items*subitems);
	if nargout > subitems						% als de extra parameter
		s2=frequency;							% ook gevraagd wordt
	end;
elseif subitems==2
	s1=data(:,1:2:items*subitems);
	s2=data(:,2:2:items*subitems);
	if nargout > subitems						% als de extra parameter
		s3=frequency;							% ook gevraagd wordt
	end;
elseif subitems==3
	s1=data(:,1:3:items*subitems);
	s2=data(:,2:3:items*subitems);
	s3=data(:,3:3:items*subitems);
	if nargout > subitems						% als de extra parameter
		s4=frequency;							% ook gevraagd wordt
	end;
end;

fclose(fid);