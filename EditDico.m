function [ dico ] = EditDico( dico, song, distSN, diconodes, syllab, featur, ce )
% display and allow edition of the dictionary

if nargin<7 
    ce = 0;
end

if ~exist('featur','var')
    featur='seqvect' ;
end

if (iscell(dico)) % convert the dico into an array
    dico = cell2struct(dico,'example')' ;
end


%panelColor = get(0,'DefaultUicontrolBackgroundColor');
panelColor = [0.5 0.5 0.5] ;

% create a colormap from white to black, with the first 24 white (to avoid noise)
map = zeros(64,3) ;
for n=1:24, map(n,:)=[1 1 1]; end;
for n=25:64, map(n,:)=[1-(n-24)/40 1-(n-24)/40 1-(n-24)/40]; end

%%% ------------ Callback Functions ---------------

% Figure resize function
function figResize(src,evt)
    fpos = get(dicof,'Position');
    set(botPanel1,'Position',...
        [0.05 0.05 fpos(3)*64/150 fpos(4)*27/54])
    set(botPanel2,'Position',...
        [fpos(3)*64.05/150 0.05 fpos(3)*85.9/150 fpos(4)*27/54])
    set(leftPanelty,'Position',...
        [0.05 fpos(4)*27/54 fpos(3)*32/150 fpos(4)*27/54])
    set(leftPanelex,'Position',...
        [fpos(3)*32.05/150 fpos(4)*27/54 fpos(3)*32/150 fpos(4)*27/54])
    set(centerPanel,'Position',...
        [fpos(3)*64.05/150 fpos(4)*27/54 fpos(3)*85.9/150 fpos(4)*27/54]);
end
    
% Bottom panel resize function
function botPanel1Resize(src, evt)
    rpos = get(botPanel1,'Position');
    set(listBoxall,'Position',...
        [rpos(3)*4/64 rpos(4)*4/27 rpos(3)*24/64 rpos(4)*20/27]);
    set(listBoxLabelall,'Position',...
        [rpos(3)*4/64 rpos(4)*24/27 rpos(3)*24/64 rpos(4)*2/27]);
end

% left panel type resize function
function leftPaneltyResize(src,evt)
    rpos = get(leftPanelty,'Position');
    set(listBoxty,'Position',...
        [rpos(3)*4/32 rpos(4)*4/27 rpos(3)*24/32 rpos(4)*20/27]);
    set(listBoxLabelty,'Position',...
        [rpos(3)*4/32 rpos(4)*24/27 rpos(3)*24/32 rpos(4)*2/27]);
end

% left panel example resize function
function leftPanelexResize(src,evt)
    rpos = get(leftPanelex,'Position');
    set(listBoxex,'Position',...
        [rpos(3)*4/32 rpos(4)*4/27 rpos(3)*24/32 rpos(4)*20/27]);
    set(listBoxLabelex,'Position',...
        [rpos(3)*4/32 rpos(4)*24/27 rpos(3)*24/32 rpos(4)*2/27]);
end

%%% Callback for list box note types
  function listBoxCallbackty(src,evt)
        set(0,'CurrentFigure',dicof);
        set(dicof,'CurrentAxes',a);
        set(listBoxex,'Value',1);
        listType=zeros(1,size(dico,2)) ;
        for i=1:size(dico,2)
            listType(i)=i ;
        end
        set(src,'String',listType);
        listBoxCallbackex(listBoxex,[]);
    end % listBoxCallbackty

%%% Callback for list box example
  function listBoxCallbackex(src,evt)
        % find current type selected
        types=get(listBoxty,'String');
        indexty=get(listBoxty,'Value');
        type=types(indexty,1:size(types,2)); % watchout multidimensional array
        type=str2num(type);
        if (numel(dico(type).example)>0)
            % find the example list
            listExample=[];
            for i=1:size(dico(type).example,2)
                listExample(i)=dico(type).example(i); 
            end
            % find the current example
            set(src,'String',listExample);
            examples=get(listBoxex,'String');
            indexex=get(listBoxex,'Value');
            example=examples(indexex,1:size(examples,2)); % watchout multidimensional array
            example=str2num(example);
            set(0,'CurrentFigure',dicof);
            axes(a); % Make the GUI axes current
            % create plot
            if (numel(get(0,'CurrentFigure'))==0) % does not create a new figure if no current
                return;
            end
            if numel(song)==0
                ShowS( syllab(example) ) ;
            else
                % find the index of the song which has this syllable in its sequence
                sgid = cell2mat(cellfun(@(x) sum((x-example)==0),{song.sequence},'UniformOutput',0)) ;
                DisplaySong( song(sgid>0), find(song(sgid>0).sequence==example), '', gca ) ;
                title(strrep([song(sgid>0).filename ' - ' num2str(find(song(sgid>0).sequence==example))],'_','\_')) ;
                colormap(map);
            end
            % get the other note (in all list) to calculate the comparison score
            alls=get(listBoxall,'String');
            indexall=get(listBoxall,'Value');
            all=alls(indexall,1:size(alls,2)); % watchout multidimensional array
            all=str2num(all);
            % compute comparison value and show it
%             scorval=DTWaverage( syllab(example).(featur), syllab(all).(featur), 1, 0, ce );
%             set(ScoreBox,'String',scorval);
        else % no examples for this type
            set(src,'String',[]);
            cla(a);
        end
    end % listBoxCallbackex

%%% Callback for list box all
  function listBoxCallbackall(src,evt)
        listNote = [dico.example] ;
        listNote = sort(listNote) ;
        set(src,'String',listNote);
        % get the current selected note
        notes=get(listBoxall,'String');
        index=get(listBoxall,'Value');
        note=notes(index,1:size(notes,2)); % watchout multidimensional array
        note=str2num(note);
        set(0,'CurrentFigure',dicof);
        axes(b); % Make the GUI axes current
        % create plot
        if (numel(get(0,'CurrentFigure'))==0) % does not create a new figure if no current
            return;
        end
        if numel(song)==0
            ShowS( syllab(note) ) ;
        else
            % find the index of the song which has this syllable in its sequence
            sgid = cell2mat(cellfun(@(x) sum((x-note)==0),{song.sequence},'UniformOutput',0)) ;
            DisplaySong( song(sgid>0), find(song(sgid>0).sequence==note), '', gca ) ;
            title(strrep([song(sgid>0).filename ' - ' num2str(find(song(sgid>0).sequence==note))],'_','\_')) ;
            colormap(map);
        end
        % get the dico number of this exmaple
        diconum = find(cell2mat(cellfun(@(x) numel(find(x==note))>0,{dico.example},'UniformOutput',false))>0) ;
        set(InfoBox,'String',['cluster ' num2str(diconum)]);
        % get the membership of this note to this cluster (low values == higher confidence)
        if exist('distSN','var')
            set(InfoBox,'String',[get(InfoBox,'String') ' (' num2str( distSN(note,diconodes(diconum))/sum((distSN(note,diconodes))) ) ')']) ;
        end
        % get the other note (in example list) to calculate the comparison score
        examples=get(listBoxex,'String');
        indexex=get(listBoxex,'Value');
        if (size(examples,2)>0)
            example=examples(indexex,1:size(examples,2)); % watchout multidimensional array
            example=str2num(example);
            % compute comparison value and show it
             scorval=DTWaverage( syllab(note).(featur), syllab(example).(featur), 1, 0, ce );
             set(ScoreBox,'String',scorval);
        end
  end

%%% Callback for graph1
  function Play1(src,evt)
      if (numel(get(listBoxex,'String'))==0) return; end
      if ( strcmp(get(src,'String'),'Play')==1 || strcmp(evt.Key,'return')==1 ) % button Play or Key Enter
        notes=get(listBoxex,'String');
        index=get(listBoxex,'Value');
        note=notes(index,1:size(notes,2)); % watchout multidimensional array
        note=str2num(note);
        PlaySyllab( syllab(note).signal, syllab(note).Fs );
      end
  end
%%% Callback for graph2
  function Play2(src,evt)
      if (numel(get(listBoxall,'String'))==0) return; end
      if ( strcmp(get(src,'String'),'Play')==1 || strcmp(evt.Key,'return')==1 ) % button Play or Key Enter
        notes=get(listBoxall,'String');
        index=get(listBoxall,'Value');
        note=notes(index,1:size(notes,2)); % watchout multidimensional array
        note=str2num(note);
        PlaySyllab( syllab(note).signal, syllab(note).Fs );
      end
  end

% play .wav files
function PlaySyllab( y, Fs )
    snd = audioplayer(y,Fs) ;
    playblocking(snd) ;
end

% Show a spectrogram
function ShowS( syllab )
    imagesc(syllab.(featur)) ;
    colormap(bone) ;
    axis xy off; % need to show proper axes
    %set(gca,'TickDir','out') ; % sets the axes tick outside the figure
    set(gca,'Unit','pixels') ;
    posit = get(gca,'Position') ;
    set(gca,'Position',[posit(1) posit(2) 0.2*numel(syllab.(featur)) posit(4)]) ;
    if ( isfield(syllab,'filename')==1 )
        L = length(syllab.filename) ;
        title([ '...' syllab.filename(max(1,L-40):L) ]) ;
    end
end
    
%%% Callback for popUpFeatur
function popUpFeaturCallback(src,evt)
    val = get(src,'Value') ;
    string_list = get(src,'String') ;
    featur = string_list(val,:) ;
    featur = strrep(featur,' ','') ; % remove spaces
end

%%% ------------ GUI layout ---------------

%%% Set up the figure and defaults
dicof = figure('Units','characters',...
        'Position',[30 30 150 54],...
        'Color',panelColor,...
        'HandleVisibility','callback',...
        'IntegerHandle','off',...
        'Renderer','painters',...
        'Toolbar','figure',...
        'NumberTitle','off',...
        'Name','Cluster analysis',...
        'ResizeFcn',@figResize);

%%% Create the left bottom uipanel
botPanel1 = uipanel('BorderType','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[0.05 0.05 64 27],...
    'Parent',dicof,...
    'ResizeFcn',@botPanel1Resize);

%%% Create the right bottom uipanel
botPanel2 = uipanel('BorderType','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[64.05 0.05 85.9 27],...
    'Parent',dicof);

%%% Create the first left side panel
leftPanelty = uipanel('bordertype','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[0.05 27 32 27],...
    'Parent',dicof,...
    'ResizeFcn',@leftPaneltyResize);

%%% Create the second left side panel
leftPanelex = uipanel('bordertype','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position',[32.05 27 32 27],...
    'Parent',dicof,...
    'ResizeFcn',@leftPanelexResize);

%%% Create the center panel
centerPanel = uipanel('bordertype','etchedin',...
    'BackgroundColor',panelColor,...
    'Units','characters',...
    'Position', [64.05 27 85.9 27],...
    'Parent',dicof);

%%% Add an axes to the center panel
a = axes('parent',centerPanel);
b = axes('parent',botPanel2);

%%% Add listbox and label
listBoxLabelty = uicontrol(dicof,'Style','text','Units','characters',...
        'Position',[4 24 24 2],...
        'String','Cluster',...
        'BackgroundColor',panelColor,...
        'Parent',leftPanelty);
% !!!! uicontrol listBox: If Max - Min > 1, then multiple item selection allowed If Max - Min <= 1, then not allowed 
listBoxty = uicontrol(dicof,'Style','listbox','Units','characters',...
        'Position',[4 4 24 20],...
        'BackgroundColor','white',...
        'Max',2,'Min',1,...
        'Parent',leftPanelty,...
        'Callback',@listBoxCallbackty,...
        'KeyPressFcn',@Play1);
listBoxLabelex = uicontrol(dicof,'Style','text','Units','characters',...
        'Position',[4 24 24 2],...
        'String','Object',...
        'BackgroundColor',panelColor,...
        'Parent',leftPanelex);
listBoxex = uicontrol(dicof,'Style','listbox','Units','characters',...
        'Position',[4 4 24 20],...
        'BackgroundColor','white',...
        'Max',2,'Min',1,...
        'Parent',leftPanelex,...
        'Callback',@listBoxCallbackex,...
        'KeyPressFcn',@Play1);
listBoxLabelall = uicontrol(dicof,'Style','text','Units','characters',...
        'Position',[4 24 24 2],...
        'String','Objects',...
        'BackgroundColor',panelColor,...
        'Parent',botPanel1);
listBoxall = uicontrol(dicof,'Style','listbox','Units','characters',...
        'Position',[4 4 24 20],...
        'BackgroundColor','white',...
        'Max',2,'Min',1,...
        'Parent',botPanel1,...
        'Callback',@listBoxCallbackall,...
        'KeyPressFcn',@Play2);

if exist('syllab','var') && numel(syllab)>0
    if isfield(syllab,'signal')
        popUpLabel = uicontrol(dicof,'Style','pushbutton','Units','characters',...
                'Position',[75.2 0.05 10 1.5],...
                'String','Play',...
                'BackgroundColor','white',...
                'Parent',centerPanel,...
                'Callback',@Play1);
        popUpLabel = uicontrol(dicof,'Style','pushbutton','Units','characters',...
                'Position',[75.2 0.05 10 1.5],...
                'String','Play',...
                'BackgroundColor','white',...
                'Parent',botPanel2,...
                'Callback',@Play2);
    end
    ScoreBox = uicontrol(dicof,'Style','text','Units','characters',...
            'Position',[50 23 10 1.5],...
            'String','',...
            'BackgroundColor',panelColor,...
            'Parent',botPanel1);
    fields = fieldnames(syllab) ;
    popUpFeatur = uicontrol(dicof,'Style','popupmenu','Units','characters',...
            'Position',[50 19 12 1.5],...
            'String',sprintf('|%s',fields{:}),...
            'BackgroundColor',panelColor,...
            'Parent',botPanel1,...
            'Callback',@popUpFeaturCallback);
end
InfoBox = uicontrol(dicof,'Style','text','Units','characters',...
        'Position',[4 1 50 1.5],...
        'String','',...
        'BackgroundColor',panelColor,...
        'HorizontalAlignment','left',...
        'Parent',botPanel1);
        
%%% Initialize list box
    listBoxCallbackall(listBoxall,[]);
    listBoxCallbackty(listBoxty,[]);
end 

