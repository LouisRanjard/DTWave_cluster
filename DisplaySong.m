function [ zoomchange ] = DisplaySong( song, sylnum, directory, TCA, timera, freqra )
% display the spectrogram of a song with the syllable boundaries marked
% Spectrogram: Hamming window of 512 samples long, 50% overlap
% TCA is The Current Axe
% Example showing one song, just knowing its filename:
% DisplaySong(song(find(strcmp({song.filename},'38_1_5.wav')==1)));
% timera(1): slider value, percent of the length of the song to display
% timera(2): zoom value, number of seconds to display
% timera(3): previous slider value
% timera(4): previous zoom value
%       example (show 100% and 200ms) : DisplaySong( song(1), 0, '', '', [100 200 0 0]); 
% freqrange(1): the minimum frequency to display
% freqrange(2): the maximum frequency to display
% sylnum: to display only a single syllable in a song segmented

if nargin>=3 && numel(directory)>0
    filename=fullfile(directory,song.filename) ;
% look for the file in the system
elseif nargin<3 || numel(directory)==0
	[status, filename] = system(['locate -n1 -r ?*/' song.filename '$']);
	if numel(strfind(filename,song.filename)==1)
		filename=filename(1:end-1) ; %remove the \n at the end
	else
		error('file not found, please give the directory path');
	end
end

if exist('timera','var')==0
    timera(1) = 1 ; % by default, show the beginning of the song
    timera(2) = 10 ; % by default, show 10 seconds
    timera(3) = 1 ; % by default, set a previous val
    timera(4) = 1 ; % by default, set a previous val
end

% watch out  Out of memory
siz =  wavread(filename,'size') ;
if siz(1)<1000000 % limit to 1 million samples but should check memory available really
    [ y , Fs ] = wavread(filename);
else
    zoomchange = [] ;
    fprintf(1,'wavfile too big for loading\n') ;
    return ;
end

% look if the file has been processed yet
if isfield(song,'SyllableS')==0
    if exist([filename '.mat'],'file') % first look in the same directory
        load([filename '.mat'],'song');
    else % look for the file in the system
        [status, filename] = system(['locate -n1 -r ?*/' [song.filename '.mat'] '$']);
        if numel(strfind(filename,[song.filename '.mat'])==1)
            filename=filename(1:round(end-1)) ; %remove the \n at the end
            load(filename,'song');
        end
    end
end
    
if isfield(song,'sequence')==0 
    song.sequence=''; song.SyllableS=[]; song.SyllableE=[];
else 
    fprintf(1,'%s %s\n',song.filename,num2str(song.sequence)) ;
end
    
% conversion from stereo to mono if required
y = y(:,1) ;
fenetre = 512 ; % window size

% show a single syllable?
if nargin>=2 && sylnum>0
    y = y(song.SyllableS(sylnum):song.SyllableE(sylnum)) ;
    fenetre = 128 ; % window size smaller for single syllables
end

% LIMIT THE TIME ACCORDING TO THE PARAMETERS
deb = max( 1 , round((length(y)-Fs*timera(2))*(timera(1)/100))) ;
fin = min( length(y) , round(deb+Fs*timera(2))) ;
% check if the previous value for zoom was equal to the total length, in
% which case there is not any change of the zoom
if fin == min( length(y) , round((timera(1)/100)*length(y)+Fs*timera(4)) )
    zoomchange=0; 
else zoomchange=1;
end
% restrict the signal to show only the zoom part
y = y(deb:fin) ;

% GET THE SPECTROGRAM
[sp,f,t] = spectrogram(y,min(fenetre,length(y)),(min(fenetre,length(y)))/2,min(fenetre,length(y)),Fs) ;
t = [0 t] ;
% LIMIT THE FREQUENCIES ACCORDING TO THE PARAMETERS
if (nargin>=6)
    L = round( (freqra(2)*length(f)) / (Fs/2) ) ;
    l = max( round( (freqra(1)*length(f)) / (Fs/2) ) , 1 ) ; % max: to avoid 0
    sp = sp( l:L, 1:end ) ;
end

%     % GET THE SPECTRAL DERIVATIVES
%     params.Fs=Fs; % sampling frequency
%     params.fpass=[500 10000]; % band of frequencies to be kept
%     params.tapers=[3 5]; % taper parameters [time bandwidth parameter, number of tapers kept (approx. 2*NW-1)]
%     params.pad=2; % pad factor for fft
%     movingwin=[0.015 0.005]; % [window winstep] i.e length of moving window and step size (factor of Fs)
%     params.trialave=1;
%     % multitaper spectrogram derivatives
%     phi=[0 pi/2];
%     [dS,t,f]=mtdspecgramc(y,movingwin,phi,params);
%     sp = squeeze(dS(2,:,:))' ;

sp = abs(sp) ; % to avoid complex values
sp(sp==0) = min(sp(sp>0)) ; % finds the minimum value of sp which greater than 0 and
sp = 20*log10(sp) ;         % affects it to sp values equal to 0 (step required to avoid log(0))

if nargin<3 || numel(TCA)==0
    F=figure( 'Name', song.filename ) ;
    ima=imagesc(sp);axis xy;
    title( num2str(song.sequence) ) ;
    scrsz = get(0,'ScreenSize'); % gets screen size [left, bottom, width, height]
    set(F,'Position',[scrsz(3)/6 scrsz(4)/2 scrsz(3)/1.5 scrsz(4)/4]) ; % sets figure position and size [left bottom width height]
else
    set(gcf,'CurrentAxes',TCA);
    ima=imagesc(sp);axis xy;
    title( num2str(song.sequence) ) ;
end

h=gca; % gets the axes of the current figure
%set(h,'TickDir','out') ; % sets the axes tick outside the figure
set(ima,'ButtonDownFcn', { @PlaySyllab, y, Fs });
colormap('bone') ;

%set( gca,'YTick',round(1:length(f)/10:length(f)) ) ;
%set( gca,'YTickLabel',round( f(round(1:length(f)/10:length(f))) ) ) ;% Hz
%set( gca,'XTick',round(1:length(t)/10:length(t)) ) ; 
%set( gca,'XTickLabel',round( t(round(1:length(t)/10:length(t)))*1000 ) ) ;% ms
%freqticks=[2000 4000 6000 8000 10000 12000 14000 16000 18000 20000 22050];
for n=1:10, freqticks(n) = round( n*Fs/20 ) ; end
freqtickidx=zeros(1,numel(freqticks));
for ft=1:numel(freqticks)
    freqtickidx(ft)=find(f<=freqticks(ft),1,'last');
end
set(gca,'ytick',freqtickidx);
set(gca,'yticklabel',freqticks); % Hz
ylabel('Frequency (Hz)');

% time ticks
timestep=round( (size(y,1)/Fs)*100 ) ; % 1/10 of the length in ms
timeticks=timestep:timestep:(t(end)*1000);
timetickidx=zeros(1,numel(timeticks));
for tt=1:numel(timeticks)
    timetickidx(tt)=find((t.*1000)<=timeticks(tt),1,'last') ;
end
set(gca,'xtick',timetickidx);
set(gca,'xticklabel',timeticks+round((deb*1000/Fs))); % ms
xlabel('Time (ms)');

% show boundaries (only if showing a song not for single syllable)
if nargin<2 || sylnum==0
    hold on ;
    SS=zeros(1,size(sp,2)) ;
    sylstart = song.SyllableS(song.SyllableS>deb)-deb+1 ;
    SS(ceil((sylstart*size(sp,2))/size(y,1))) = size(sp,1) ; 
    h = stem(SS);
    set(h,'Color','blue','Marker','none') ;
    SE = zeros(1,size(sp,2)) ; 
    sylstop = song.SyllableE(song.SyllableE>deb)-deb+1 ;
    SE(ceil((sylstop*size(sp,2))/size(y,1))) = size(sp,1) ; 
    h = stem(SE) ;
    set(h,'Color','red','Marker','none') ;
    hold off ;
    drawnow ;
end

% save the spectro as an image file
% set(gcf,'PaperPositionMode','auto');
% print('-dtiff','-zbuffer','-r100',['/home/lran022/Sounds/test_Hen&Chicken/result2006621936+SOManalysis/SOMtypeS/EncodedSong_SOM10_ClustDist/' song.filename '.tiff']) ;
% close(F);

function PlaySyllab( src, eventdata, y, Fs )
    % play a sound
    global ap ;
    if (isa(ap, 'audioplayer'))
        if (isplaying(ap))
            stop(ap);
        else
            play(ap) ;
        end
    else
        ap = audioplayer(y,Fs) ;
        playblocking(ap) ;
    end
    clear global ap;
end

end
    
