function y = loadwave(s)

%
%  rf = loadwave('filename');
%
%   'filename' -- the name of a signa 5x waveform file
%   rf         -- waveform extracted from that file


fip = fopen(s,'r','b');
if fip == -1,
   disp(sprintf('Error opening %s for read',s));
   return;
end;

d = fread(fip,'int16');
y = d(33:length(d)-4);

