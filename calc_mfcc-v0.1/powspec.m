function y = powspec(x, sr, wintime, steptime, dither, preemph)
%y = powspec(x, sr, wintime, steptime, sumlin, dither)
%
% compute the powerspectrum of the input signal.
% basically outputs a power spectrogram
%
% each column represents a power spectrum for a given frame
% each row represents a frequency
%
% default values:
% sr = 8000Hz
% wintime = 25ms (200 samps)
% steptime = 10ms (80 samps)
% which means use 256 point fft
% hamming window

% for sr = 8000
%NFFT = 256;
%NOVERLAP = 120;
%SAMPRATE = 8000;
%WINDOW = hamming(200);

if nargin < 2
  sr = 8000;
end
if nargin < 3
  wintime = 0.025;
end
if nargin < 4
  steptime = 0.010;
end
if nargin < 5
  dither = 1;
end

winpts = round(wintime*sr);
steppts = round(steptime*sr);

NFFT = 2^(ceil(log(winpts)/log(2)));
% hamming back for better HTK
WINDOW = hamming(winpts)';
% even hanning is the original standard???
%WINDOW = [0,hanning(winpts-1)'];
% hanning gives much less noisy sidelobes
NOVERLAP = winpts - steppts;
SAMPRATE = sr;

% Values coming out of rasta treat samples as integers, 
% not range -1..1, hence scale up here to match (approx)
%y = abs(specgram(x,NFFT,SAMPRATE,WINDOW,NOVERLAP)).^2;

%if preemph ~= 0
%  x = filter([1 -preemph], 1, x);
%end

ncols = 1+ floor((length(x)-winpts)/steppts);
xx = x(repmat(steppts*[0:ncols-1],winpts,1)+repmat([1:winpts]',1,ncols));
% HTK: Each window is made zero-mean
xx = xx - repmat(mean(xx),winpts,1);
% pre-emphasis is per segment, as if first point is repeated
% (per HSigP.c:134)
xx = filter([1 -preemph], 1, [xx(1,:);xx]);
xx = xx(2:end,:);
xx = repmat(WINDOW',1,ncols).*xx;
% Up to this point verifies against gdb of HCopy (at least for
% first frame) out to 6sf

y = (fft(xx,NFFT));
y = y(1:(NFFT/2+1),:);

%y(1:5,1)  % good to 4sf or so; d.c. only good to <3sf

y = abs(y).^2;

%y(1:5,1)

%x(steppts+[1:5])'.*WINDOW(1:5)

%y(1:5,1)'

% imagine we had random dither that had a variance of 1 sample 
% step and a white spectrum.  That's like (in expectation, anyway)
% adding a constant value to every bin (to avoid digital zero)
if (dither)
  y = y + winpts;
end
% ignoring the hamming window, total power would be = #pts
% I think this doesn't quite make sense, but it's what rasta/powspec.c does

% that's all she wrote
