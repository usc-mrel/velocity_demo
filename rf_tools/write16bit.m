%
%function write3dshorts(fname, array)
%
%	Writes the magnitude of the array to a test file.
%

function write3dshorts(fname, array)

s = size(array);

array = round(array);	% round to integer.

% -- Check LSB for EPIC waveforms.
f = find(mod(array,2) == 1);
if (length(f)>1)
  disp('Warning - for EPIC waveforms - More than one point has LSB set.');
end;
if (mod(array(end),2)==0)
  disp('Warning - for EPIC waveforms - Last point does not have LSB set.');
end;

	
fid = fopen(fname,'w','l');
if (nargin < 2)
	fwrite(fid,s,'int16');
end;
fwrite(fid,array,'int16');
fclose(fid);



