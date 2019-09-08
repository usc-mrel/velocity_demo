addpath ~/lab4


%	This is a simple matlab script that does the following:
%	  -Load a pfile header
%	  -Load the data from a pfile
%	  -Reconstruct the images from each coil and each slice (FFT)
%	  -**Replace the image data with something else:  A square of pixels
%		in one corner.
%	  -**You will change that part!
%
%	For simplicity, the 
disp('loading data');
[header,r,byteheader] = rawheadX('~/lab4/Pfile');	% Load the file we copied.
[data,header] = rawloadX('~/lab4/Pfile');	% Load the file we copied.
fip = fopen('~/lab4/Pfile2');
if (fip>=0)
  fclose(fip);
  [data2] = rawloadX('~/lab4/Pfile2');		% Load 2nd file if it exists
else
  data2 = data;		% Just copy so that script doesn't fail...
end;

% Assume data is (Nx,Ny,Ns,1,Nc)

Nc = header.ncoils;	% -- Number of coils being used
Ns = header.nslices;	% -- Number of slices being used.
scale = sqrt(header.nframes*header.frsize);	% -- Scale for FFT
for c=1:Nc	% -- Loop through coils

 %tt=sprintf('Coil %d of %d',c,Nc); disp(tt);	% -- Show status
 for s = 1:Ns	% -- Loop through slices

  %tt=sprintf('  Slice %d of %d',s,Ns); disp(tt);	% -- Show status
  im = fftshift(fft2(fftshift(squeeze(data(:,:,1,s,c)))));  % 2DFFT
  im = im/scale;	% FFT scale

  % -- Put a "frame" around the data, 10 pixels wide of value 1000. 
  im(1:10,:)=1000; im(:,1:10)=1000;
  im(end-9:end,:)=1000; im(:,end-9:end)=1000;

  % -- ** Add code to do something useful here, like calculate a 
  % -- ** flip angle (say) 

  img(:,:,1,s,c) = im;		% -- Put back into a big image array.
 end;
end;


% -- Write data to a new pfile of identical size.

disp('writing data');	
writepfileX('~/lab4/Pfile.out',byteheader,img(:,:,:,:,:),header,1);



disp('Exiting Matlab - comment out the last line to avoid this!');
%exit;	% For automated reconstruction.




