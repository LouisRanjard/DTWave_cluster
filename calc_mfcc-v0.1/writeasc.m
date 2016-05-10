function writeasc(fname, data)
% writeasc(fname, data)
%   Write out data to file fname in ascii format.
%   Each row of data is a row of output, 
%   but each row begins with an utterance number (always zero) and
%   a frame number (starts at zero).
% 2013-02-26 Dan Ellis dpwe@ee.columbia.edu

[nr, nc] = size(data);

fp = fopen(fname, 'w');

if fp <= 0
  error(['writeasc: could not write to ', fname]);
end

for i = 1:nr
  fprintf(fp, '%d %d', 0, i-1);
  fprintf(fp, ' %.6f', data(i,:));
  fprintf(fp, '\n');
end

fclose(fp);
