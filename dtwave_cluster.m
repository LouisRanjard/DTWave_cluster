% ----------------------------------------------------------------------------
%      dtwave_cluster: unsupervised sound classification
%              Copyright (C) The University of Auckland 2008
%           (Louis Ranjard, l.ranjard at auckland dot ac dot nz)
%
% dtwave_cluster is distributed under the terms of the
%                 GNU General Public License Version 2
% A copy of this license is included in the file gpl2.txt of this distribution
% -----------------------------------------------------------------------------

function [syllab1,tree1,diconodes1,classrep,DavBou] = dtwave_cluster( wavdir, varargin )
%

tic;

% arguments check, parse out the optional arguments
[htkconfigdir, verbose, outputfile,...
 wintime, hoptime, numcep, lifterexp, sumpower, preemph, dither, ...
 minfreq, maxfreq, nbands, bwidth, dcttype, fbtype, usecmp, modelorder, broaden,...
 epoch, LR, NS, NC, thexpan, gama, ce, nrepeat, depth, wHamming] = ...
process_options(varargin, 'htkconfigdir', '', 'verbose', 0, 'outputfile', fullfile(wavdir,'dtwave_cluster.csv'),...
                          'wintime', 0.3, 'hoptime', 0.1,...
                          'numcep', 13, 'lifterexp', -22, 'sumpower', 1, 'preemph', 0.97,...
	                      'dither', 0, 'minfreq', 300, 'maxfreq', 20000,...
	                      'nbands', 26, 'bwidth', 1.0, 'dcttype', 3,...
	                      'fbtype', 'htkmel', 'usecmp', 0, 'modelorder', 0, 'broaden', 0,...
                          'epoch', 4, 'LR', [0.5 0.05],...
                          'NS', [2 1], 'NC', [4 2],...
                          'thexpan', -1,'gama', 0.95,...
                          'ce', 0,'nrepeat', 5, 'depth', 0, 'wHamming', 0);

% encode wav files
if isempty(htkconfigdir)
    syllab1 = encodeWavlist(wavdir,[],[],[],1,0,'wintime', wintime, 'hoptime', hoptime,...
                            'numcep', numcep, 'lifterexp', lifterexp, 'sumpower', sumpower,...
                            'preemph', preemph, 'dither', dither, 'minfreq', minfreq,...
                            'maxfreq', maxfreq, 'nbands', nbands, 'bwidth', bwidth,... 
                            'dcttype', dcttype, 'fbtype', fbtype, 'usecmp', usecmp,...
                            'modelorder', modelorder, 'broaden', broaden);
else
    syllab1 = encodeWavlist(wavdir,htkconfigdir,[],[],1,0);
end
if thexpan<0, thexpan=round(0.05*numel(syllab1)); end

% write parameters to file?
%fid = fopen( fullfile(etreedir,'etree_parameters.txt'), 'w') ;
%fprintf(fid,'# %s\n\nnumber of epochs %f\nlearning rate %f %f\nneighborhood size %f %f\nnumber of children %f %f\nthreshold before expansion %f\nweight decay %f\n',...
%    datestr(now),epoch,LR(1),LR(2),NS(1),NS(2),NC(1),NC(2),thexpan,gama) ;
%fclose(fid) ;

fprintf(1,'\n****** dtwave_cluster: unsupervised sound classification ******\n') ;
if (numel(syllab1)==0) 
    error('No wav files found in %s',wavdir); 
end
fprintf(1,'\n%g WAV files loaded\n',numel(syllab1)) ;
if numel(NS)>1
    fprintf(1,'Parameters: %g, %g, %g, %g, %g, %g, %g, %g, %g\n', epoch,LR(1),LR(2),NS(1),NS(2),NC(1),NC(2),thexpan,gama ) ;
else
    fprintf(1,'Parameters: %g, %g, %g, %g, %g, %g, %g, %g\n', epoch,LR(1),LR(2),NS(1),NC(1),NC(2),thexpan,gama ) ;
end

classrep = zeros(size(syllab1,2),nrepeat) ;
DavBou = zeros(1,nrepeat) ;
for r=1:nrepeat
    fprintf(1,'\n--------------------- classification tree: %d/%d ---------------------\n',r,nrepeat) ;
    % create a classification ETree
    [tree1, weight1] = ETDTWrec(syllab1,epoch,LR,NS,NC,thexpan,gama,'seqvect',[],verbose,ce) ;
    [distSN1, distNN1] = distSNN(syllab1,weight1,'','seqvect') ;
    % do syllable library
    if depth>0
        %[dico, diconodes] = etreedico(tree,weight,[],[depth mean(reshape(distSN,1,[])) 0 0 0],[],distSN,distNN,[],1) ;
        [dico1, diconodes1, ~, ~ ,~, DavBou(r)] = etreedico(tree1,weight1,[],[depth 0 0 0 0],[],distSN1,distNN1,[],1) ;
    elseif depth==0
        [dico1, diconodes1, ~, ~ ,~, DavBou(r)] = etreedico(tree1,weight1,[],[depth 0 0 mean(reshape(distNN1,1,[]))/2 0],[],distSN1,distNN1,[],3) ;
    elseif depth<0 && usejava('desktop') % the desktop is available
        OptimDepth( tree1, weight1, '', distSN1, distNN1, 'seqvect' );
        object_handle = gcf;
        data = guidata(object_handle);
        %%% UI control to create a dictionary from a click on the graph
        set(gcf,'name','select the depth for the classification') ;
        set(gcf,'windowbuttondownfcn',{@getdepth,tree1,weight1,distSN1,distNN1,[],[]}) ;
        waitforbuttonpress;
        guidata(object_handle,data);
        [dico1, diconodes1, ~, ~ ,~, DavBou(r)] = etreedico(tree1,weight1,[],[depth mean(reshape(distSN1,1,[])) 0 0 0],[],distSN1,distNN1,[],1) ;
    end
    % save classification
    for n1=1:size(syllab1,2)
        for t1=1:size(dico1,2)
            if find(dico1(t1).example==n1)
                classrep(n1,r) = t1 ;
                break;
            end
        end
    end
end

% create classification output file
if exist('outputfile','var') && ~isempty(outputfile)
    fid = fopen(outputfile,'w') ;
    for nsyl=1:size(syllab1,2)
        fprintf(fid, '%s,', syllab1(nsyl).filename) ;
        c = num2str(classrep(nsyl,:),'%1d,') ;
        c = strrep(c,' ','') ;
        fprintf(fid,'%s\n',c) ;
    end
    fprintf(fid, 'Davies-bouldin,') ;
    db = num2str(DavBou,'%.2f,') ;
    db = strrep(db,' ','') ;
    fprintf(fid, '%s\n',db) ;
    fclose(fid) ;
end

% plot classification
if numel(syllab1)<100
    if wHamming>0 % use weighted Hamming distance, using 1/DaviesBouldin as weight
        whamming = zeros(size(classrep,1)) ;
        maxwhamming = sum(ones(1,size(classrep,2))./DavBou) ;
        for i=1:size(classrep,1)
           for j=1:size(classrep,1)
               whamming(i,j) = ( maxwhamming - sum((classrep(i,:)==classrep(j,:))./DavBou) ) ./ maxwhamming ;
           end
        end
        Zdist = linkage(squareform(whamming)) ;
    else
        Zdist = linkage(classrep,'average','hamming') ;
    end
    dendrogram(Zdist,numel(syllab1),'Labels',{syllab1.filename},'Orientation','left') ;
end
    
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% core functions

function [seqvect1] = norma_seqvect(seqvect,vectstruct)
% vectstruct contains for each group the set of line to consider
% for example (if MFCC_O_D_A)
% vectstruct=[1 12; 13 13; 14 25; 26 26; 27 38; 39 39; 40 51;52 52;101 112;113 124;125 136;137 148;149 160;161 172;173 184;185 196;197 208;209 220;221 232;233 244;]
    tmp=seqvect;
    if max(vectstruct)~=size(seqvect,1)
        fprintf(1,'vector size: %f %f\nnorma vector size: %f\n',size(seqvect,1),size(seqvect,2),max(max(vectstruct)));
        error('norma_seqvect(): vector size incompatible with vector structure');
    end
    if nargin>1
        % normalise subgroups of the cell according to vectstruct
        cmdtext = '';
        for n2=1:size(vectstruct,1)
            cmdtext = [cmdtext 'norma(tmp(' num2str(vectstruct(n2,1)) ':' num2str(vectstruct(n2,2)) ',:),0,1);'] ;
            % cmdtext = [cmdtext 'zscore(tmp(' num2str(vectstruct(n,1)) ':' num2str(vectstruct(n,2)) ',:));'] ;
        end
        eval(['seqvect1 = [' cmdtext '];']) ;
    else
        % normalize every coefficient (line) independently
        seqvect1=zeros(size(seqvect,1),size(seqvect,2));
        for nc=1:numel(tmp(:,1))
          seqvect1(nc,:)=norma(tmp(nc,:),0,1);
          % seqvect2(nc,:)=zscore(tmp(nc,:));
        end
    end
end

function [ normalized ]= norma (signal,mini,maxi)
% normalize a signal between mini and maxi
%if matrix, mini and maxi are calculated in 2D
    [rows, cols]=size(signal);
    if (rows==1 || cols==1)
        miniS=min(signal);
        maxiS=max(signal);
    else
        miniS=min(min(signal));
        maxiS=max(max(signal));
    end
    normalized = ((signal-miniS)*(maxi-mini) ./ (maxiS-miniS)) + mini ;
end

function [tree, weight, BMU] = ETDTWrec(syllab,epoch,LR,NS,NC,thexpan,gama,fieldnam,dico,verbose,ce)
% performs a Evolving tree clustering analysis (Samarasinghe2006, Pakkanen2004)
% syllab.'fieldnam' are the matrices to classify
% epoch is the number of Epochs
% LR(1) is the initial Learning Rate decreasing until LR(2)
% NS is the neighborhood strength, if more than one element then use a
%     bell-shaped function instead of a constant value
%     NS(1) is the maximum value obtained at middle point of the learning process
%     NS(2) is the standard deviation of the function, if ==0 then NS(2) takes (#epochs)/5
% weight contains the weigth matrix of each index of tree
% NC(1) is the initial number of children per node and NC(2) is the final number of children
% thexpan is the threshold above which a node is splitted
%       only nodes which have a cluster scatter above the mean can be divided
% gama is the weight decay on the hit counter (Pakkanen2006)
% fielnam is the name of the field in syllab
% dico is a classification to be used to improve the convergence of the clustering
%       must be decreasingly ordered according to the distance
% verbose: verbose mode? 1 or 0
% ce: use compression/exap[nsion? 1 or 0
% example ETDTWrec(syllab,5,[0.9 0.01],[4 1],[2 2],numel(syllab)*0.2,0.95,'cepvect',[],1,1)

% tree contains for each index the index of its father
% weight contains the weight matrices (cluster centres)
% BMU records the path for each sample in the ETree
% LEARNING AND DISPLAY PARAMETERS
    if nargin<11
        ce = 0 ;
        if nargin<10
            verbose = 0 ;
            if nargin<9
                dico = [] ;
                if nargin<8
                    fieldnam = 'cepvect' ;
                    if nargin<7
                        gama = 1 ;
                        if nargin<6
                            thexpan = round(numel(syllab)/10) ; % arbitrarily set the splitting threshold
                            if nargin<5
                                NC = [2 2] ; % default is binary trees
                                if nargin<4
                                    NS = [3 3] ;
                                    if nargin<3
                                        LR = [0.9 0.01] ; % minimal learning rate (Samarasinghe2006 uses a threshold of 0.01)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    % to make computation faster, record the dico index of each syllable
    if numel(dico)>0
        syldico=zeros(1,numel(syllab)) ;
        for si=1:numel(syllab)
            for di=1:numel(dico)
                if find(dico(di).example==si)
                    syldico(si)=di ;
                    break ;
                end
            end
        end
    end
    % INITIALISATION
    %%%%% mean sizes of syllab cepstrum to find the size of weight matrices
    % max 100 random samples
    idxs = ceil(numel(syllab).*rand(1,min(100,numel(syllab)))) ; % avoid extrem values
    % fprintf(1,'%d %d\n',size([syllab(idxs).(fieldnam)],1),size([syllab(idxs).(fieldnam)],2)) ;
    vals = reshape([syllab(idxs).(fieldnam)],1,[]) ;
    x = size([syllab.(fieldnam)],1) ; % frequency mean size
    y = round(size([syllab.(fieldnam)],2)/numel(syllab)) ; % time mean size
    % tree initialisation with random weights values picked from vals
    tree(1) = 0 ; % Root node
    numchildren = NC(1) ;
    for w=1:numchildren
        tree(w+1) = 1 ; % its father is the root
        rindex = 1 + round((size(vals,2)-1)*rand(x,y)) ;
        weight{w+1} = vals(rindex) ; % +1 because weight(1) is the root
        nbmu(w+1) = 0 ; % will store the number of time this node is BMU (hit counter)
    end
    clear vals ;
    % DISPLAY?
    if verbose==1
        %%% figures
        f = figure('Units','pixels','Position',[300 100 600 840]) ;
        Panel1 = uipanel('Units','pixels','Position',[20 20 270 180],'Title','Davies-Bouldin','Parent',f) ;
        Panel2 = uipanel('Units','pixels','Position',[310 20 270 180],'Title','Evolving Tree','Parent',f) ;
        Panel3 = uipanel('Units','pixels','Position',[20 220 270 180],'Title','Furthest Sample Mean Distance','Parent',f) ;
        Panel4 = uipanel('Units','pixels','Position',[310 220 270 180],'Title','Mapping Precision','Parent',f) ;
        Panel5 = uipanel('Units','pixels','Position',[20 420 270 180],'Title','Closest Sample Mean Distance','Parent',f) ;
        Panel6 = uipanel('Units','pixels','Position',[310 420 270 180],'Title','Inner Scatter Distribution','Parent',f) ;
        Panel7 = uipanel('Units','pixels','Position',[20 620 270 180],'Title','Dividing Events','Parent',f) ;
        mapprec = axes('parent',Panel4) ;
        treep = axes('parent',Panel2) ;
        daviesbouldinax = axes('parent',Panel1) ;
        innerscatterax = axes('parent',Panel6) ;
        closestsample = axes('parent',Panel5) ;
        dividingevents = axes('parent',Panel7) ;
        maxclustdist = axes('parent',Panel3) ;
    elseif verbose && usejava('desktop')
        % shows progress
        h = waitbar(0,'learning...') ;
    end
    % print the parameter values
    if verbose==1
        fprintf(1,'epoch time -- Neighborhood Size, Mapping precision, Davies-Bouldin, Tree size, Leaves number\n' ) ;
    else
        fprintf(1,'epoch time -- Neighborhood Size, Mapping precision, Tree size, Leaves number' ) ;
    end
    % some memory allocation
    BMU = zeros(numel(syllab),1,2) ; % can't know the max depth yet so just initialize with a value of 1
    p = zeros(1,epoch) ; % mapping precision
    DB = [] ; % Davies Bouldin index
    meanCS = [] ; % record the mean distance of the closest sample from each cluster centre (leaf)
    meanFS = [] ; % record the mean distance of the furthest sample from each cluster centre
    dividing = [] ; % record when there is a dividing event during learning
    divided = 0 ; % test to know when there has been dividing event
    % get total number of steps to implement a step counter
    steprecord = ceil(epoch*numel(syllab)*0.01) ;
    stepcount = 0 ;
    % use a bell shaped neighborhood size?
    if numel(NS)>1
        if NS(2)==0 % no standard deviation input, take the (total number of epoch)/5
            NS(2) = epoch/5 ;
        end
        % compute the factor to be used during updates
        nsfactor = NS(1)/(1/(NS(2)*sqrt(2*pi))) ;
    end    
    nstep = epoch * size(syllab,2) ;
    istep = 1 ;
    for e=1:epoch % each epoch
        currtime=fix(clock);
        fprintf(1,'\nepoch %g/%g %02g-%02g-%g %02g:%02g -- ',e,epoch,currtime(3),currtime(2),currtime(1),currtime(4),currtime(5));
        % randomize the order which syllables are picked up
        sindex = randperm(size(syllab,2)) ;
        P = 0 ;
        %%%% compute the neighboring size for this epoch
        if numel(NS)>1
            NeSt = max( 1 , (1/(NS(2)*sqrt(2*pi)))*exp(-(((e)-((epoch+1)/2)).^2)/(2*(NS(2))^2))*nsfactor ) ;
            NeSt = round(NeSt) ;
        else
            NeSt = NS ;
        end
        % every syllable
        for s=sindex
            %%%%% restraint due to distance from BMU and learning restraint due to time
            % time constant that brings the learning rate and the neighborhood size to a very small value with iterations
            ED = exp(-e^2/(epoch*0.75)^2) ;
            %%%%% LEARNING RATE
            Ftime = max( LR(1)*ED , LR(2) ) ; % exponentially decreases from Ftime0
            %%%%% search for BMU in the tree
            BMU(s,1,1) = 1 ; % starts at the top of the tree
            BMU(s,1,2) = 0 ;
            dep = 1 ;
            while (sum(tree==BMU(s,dep,1))>0) % this node is a father => keep on searching
                leaves = find(tree==BMU(s,dep,1)) ;
                dep = dep + 1 ;
                scor = zeros(1,max(leaves)) ;
                for L=leaves
                    [ scor(L) ] = DTWaverage( syllab(s).(fieldnam), weight{L}, 1, 0, ce ) ; % cepstrum coefficient alignment
                end
                BMU(s,dep,2) = min(scor(leaves)) ;
                BMU(s,dep,1) = find(scor==min(scor(leaves)),1) ; % find best match (Best Matching Unit), minimum distance
            end
            % re-initialised BMU(s,:,:) if it was longer before
            BMU(s,dep+1:end,1) = 0 ;
            BMU(s,dep+1:end,2) = 0 ;
            %%%%% mapping precision
            P = P + BMU(s,dep,2) ;
            %%%%% Use the samples in the same cluster to update (according to the dico)
            if numel(dico)>0 && syldico(s)>0
                cws = dico(syldico(s)).example ;
                for cwi=1:numel(dico(syldico(s)).example)
                    for L=2:numel(tree)
                        %%%%% NEIGHBORHOOD SIZE
                        Edist = TreeDist(tree,BMU(s,dep,1),L) ;
                        %Fdist = exp( (-Edist^2)/(2*(NeSt*exp(-e/(epoch*0.5)))^2) ) ; % old bug
                        Fdist = exp( (-Edist^2)/(2*(NeSt)^2) ) ;
                        H = Ftime*Fdist ;
                        %%%%% LEARNING STEP
                        %if ( H>0.01 ) % in order to save time if the factor is too small do not compute
                            [ ~, weight{L} ] = DTWaverage( syllab(cws(cwi)).(fieldnam), weight{L}, 1, H, ce ) ;
                        %end
                    end
                end
            end
            %%%%% UPDATE the tree except the root in the neighbourhood of the Best Matching Unit by pulling them closer to the input matrix
            for L=2:numel(tree)
                    %%%%% NEIGHBORHOOD SIZE
                    Edist = TreeDist(tree,BMU(s,dep,1),L) ;
                    %Fdist = exp( (-Edist^2)/(2*(NeSt*exp(-e/(epoch*0.5)))^2) ) ; % old bug
                    Fdist = exp( (-Edist^2)/(2*(NeSt)^2) ) ;
                    H = Ftime*Fdist ;
                    %%%%% LEARNING STEP
                    if ( H>0.01 ) % in order to save time if the factor is too small do not compute
                        [ ~, weight{L} ] = DTWaverage( syllab(s).(fieldnam), weight{L}, 1, H, ce ) ;
                    end
            end
            %%%%% increases nbmu 
            nbmu(BMU(s,dep,1)) = nbmu(BMU(s,dep,1))+1 ;
            %%%%% computes WITHIN CLUSTER INDEX (DaviesBouldin clustering index)
            BMU1 = BMU(:,:,1) ;
            clustidx = unique(BMU1(BMU1>0)) ;
            leaves = setdiff(1:numel(tree),unique(tree)) ; % the leaves i.e potentially clusters in the tree
            clustidx = intersect(clustidx,leaves) ; % all the nodes which are also leaves
            if numel(clustidx)>0
                S = zeros(1,max(clustidx)) ; % within cluster scatter
    %             if verbose==1
    %                 Smaxtmp = zeros(1,max(clustidx)) ; % maximum cluster distance %noneed%
    %                 D = zeros(max(clustidx),max(clustidx)) ; % between cluster separation (for DAVIES-BOULDIN) %noneed%
    %             end
                for n3=1:numel(clustidx)
                   cidx = clustidx(n3) ;
                   % within cluster scatter
                   [a, b] = find(BMU(:,:,1)==cidx) ;
                   S(cidx) = sum( BMU(a,b(1),2) ) / numel(a) ; % only consider b(1) because they're all the same depth
                   %%%%% DIVIDING STEP
                   % divide if S>mean(mapping precision) and expansion threshold reached (except for first epoch)
                   if e==1 && nbmu(cidx)>thexpan
                       divid=1;
                   elseif  nbmu(cidx)>thexpan && S(cidx)>(P(P>0)/numel(P(P>0)))
                       divid=1;
                   else
                       divid=0;
                   end
                   if divid==1
                        divided = 1 ; % record that there has been dividing for at least one cluster (for display)
                        for w=1:numchildren
                            tree = [tree cidx] ;
                            weight{numel(tree)} = (rand(size(weight{cidx},1),size(weight{cidx},2))*0.20+0.9).*weight{cidx} ; % add some random noise
                            nbmu = [nbmu 0] ;
                        end
                        % initialise the new child
                        [a, b] = find(BMU(:,:,1)==cidx) ;
                        for m=1:numel(a)
                           BMU(a(m),b(m)+1,1) = numel(tree)-floor(rand*numchildren) ; % randomly assign a child
                           BMU(a(m),b(m)+1,2) = BMU(a(m),b(m),2) ; % use the same distance 
                           % this makes DAVIES-BOULDIN worse because from one cluster with a bad inner scatter, 
                           % two are created. this is specially true at the begining of the learning process when the
                           % number of cluster is quite small.
                        end
                   end
                   % between cluster separation (for DAVIES-BOULDIN) %noneed%
    %                if verbose==1
    %                    clustidxdif = clustidx(clustidx~=cidx) ;
    %                    for n2=1:numel(clustidxdif)
    %                        cidx2 = clustidxdif(n2) ;
    %                        D(cidx,cidx2) = DTWaverage( weight{cidx}, weight{cidx2} ) ;
    %                    end
    %                end
                end
                % DISPLAY?
                if verbose==1
                    if stepcount==steprecord
                        stepcount = 0 ;
                        % COMPUTE SOME CLUSTERING STATISTICS every 0.01% update steps
                        % need to recompute all the distances because there has been updates (so not possible to use BMU)
                        [distSN, distNN] = distSNN(syllab,weight,'',fieldnam,'SNN',0) ;
                        treelayout(tree) ;
                        % depth of each node
                        depthn = zeros(1,numel(tree)) ;
                        for n4 = 1:numel(tree)
                            a=n4 ;
                            while ( tree(a)~=0 ); depthn(n4)=depthn(n4)+1; a=tree(a); end
                        end
                        leaves = setdiff(1:numel(tree),unique(tree)) ;
                        % the cluster centres are the ones of depth==limdep and the ones of inferior depth but which are leaves
                        diconodes = [ find(depthn==dep) intersect(find(depthn<dep),leaves) ] ;
                        % defines the cluster (dictionary)
                        dico = struct('example',cell(1,numel(diconodes))) ;
                        % redefine distSN and distNN with just the cluster centroids
                        distSNc = distSN(:,[1 diconodes]) ;
                        % build dico
                        for sy=1:size(distSNc,1)
                            dicoidx = find(distSNc(sy,2:end)==min(distSNc(sy,2:end)),1) ; % add 1 because we avoided the root
                            dico = AddExampleDico( dico, dicoidx, sy ) ;
                        end
                        [ cs1, fs1 ] = ClustStat(dico,diconodes,distSN) ;
                        db1 = DaviesBouldin2(dico,diconodes,distSN,distNN) ;
                        %%% PLOT THESE STATISTICS
                        meanCS = [meanCS cs1]; axes(closestsample); plot(meanCS);
                        meanFS = [meanFS fs1]; axes(maxclustdist); plot(meanFS);
                        DB = [DB db1]; axes(daviesbouldinax); plot(DB);
                        dico = [] ;
                    else
                        stepcount = stepcount+1 ;
                    end

                    % dividing events
                    if divided==1 % there has been dividing
                        divided=0;
                        dividing = [dividing 1] ;
                    else
                        dividing = [dividing 0] ;
                    end
                    axes(dividingevents);
                    plot(dividing); 
                    drawnow;
    %                 % closest sample mean distance from each cluster centre %noneed%
    %                 distSN = distSNN(syllab,weight([1 clustidx]),'',fieldnam,'SN',0) ; % watch out, distSNN() skips the first weight which is supposly the root of the ETree
    %                 meanCS = [meanCS sum(min(distSN,[],1))/numel(clustidx)] ;
    %                 axes(closestsample);
    %                 plot(meanCS);
    %                 % furthest sample mean distance from each cluster centre (need to know the clusters) %noneed%
    %                 FS = zeros(1,numel(clustidx)) ;
    %                 for cidx3=1:numel(clustidx)
    %                     [a b] = find(BMU(:,:,1)==clustidx(cidx3)) ;
    %                     FS = max(distSN(a,1+cidx3)) ;
    %                 end
    %                 meanFS = [meanFS sum(FS)/numel(clustidx)] ;
    %                 axes(maxclustdist);
    %                 plot(meanFS); 
    %                 drawnow;
                    % distribution of inner scatter
                    axes(innerscatterax);
                    hist(S(S>0),10);
    %                 % between cluster separation (for DAVIES-BOULDIN) %noneed%
    %                 R = zeros(1,max(clustidx)) ; % store the max ratio of S and D
    %                 for n=1:numel(clustidx)
    %                     cidx = clustidx(n) ;
    %                     indexdif = clustidx(clustidx~=cidx) ;
    %                     R(cidx) = max( (S(cidx)+S(indexdif))./D(cidx,indexdif) ) ;
    %                 end
    %                 % compute DAVIES-BOULDIN index %noneed%
    %                 DB(find(DB==0,1)) = sum(R)/numel(clustidx) ;
    %                 axes(daviesbouldinax); 
    %                 plot(DB(DB>0)); 
    %                 drawnow;
                elseif verbose && usejava('desktop')
                    waitbar( istep/nstep ) ;
                    istep=istep+1 ;
                end
            end
            %%%%%%%%% plot the tree with labels
            % DISPLAY?
            if verbose==1
                axes(treep) ; 
                treeplot(tree) ;
                [x,y] = treelayout(tree);
                text(x,y,num2str((1:length(tree))')); % index
                %text(x,y,num2str(nbmu(ii)')); % number of bmu
                drawnow ;
            end
        end
        % weight decay on hit counter (Pakkanen2006) after each epoch
        nbmu = round(gama.*nbmu) ;
        % expansion size decay (Louis2007)
        numchildren = max( NC(1)-(round(NC(1)/epoch)) , NC(2) ) ;
        % mapping precision
        p(e) = P/size(syllab,2) ;
        if verbose==1
            fprintf(1,'%g, %g, %g, %g, %g', NeSt, p(e), DB(end), numel(tree), numel(setdiff(1:numel(tree),unique(tree))) ) ;
        else
            fprintf(1,'%g, %g, %g, %g', NeSt, p(e), numel(tree), numel(setdiff(1:numel(tree),unique(tree))) ) ;
        end
        % DISPLAY?
        if verbose==1, axes(mapprec); plot(p); drawnow; end
        %%%%%%%%%%%%%%%% SAVING THE TEMPORARY TREE %%%%%%%%%%%%%%%%
        %save( '-v7', './tree_tmp.mat', 'tree', 'weight','e','nbmu' );
    end
    if verbose==1 && usejava('desktop')
        close(h) ;
    end
    %delete('./tree_tmp.mat') ;
    fprintf(1,'\n');
end

function [ dist, mat4 ] = DTWaverage( mat1, mat2, cow, q, ce )
% perform a dynamic time warping in order to measure a distance between matrices
% according to the algorithm of Oommen1995
% Kruskal&Liberman in Chapt4 Sankoff&Kruskal1999 (time warps...)
% combining insertion/deletion and compression/expansion
% if s_ave==1 then it computes the average between the 2 vector sequences
% q is the weight for the computation of the average (must be between 0 and 1)
%        it affects mat1 and mat2 has the weight 1-q
% cow is the weight of each coefficient for the distance computation, must be the same length than size(mat1,1)
% ce indicates whether to use compression/expansion or not
% an average sequence can be computed:
    % bak=1 expansion
    % bak=2 insertion
    % bak=3 substitution
    % bak=4 deletion
    % bak=5 compression
    % if a matrix is empty return -1
    if size(mat1,2)==0 || size(mat2,2)==0
        dist=-1;
        mat4=[];
        % error('DTWaverage(): at least of the input matrices is empty');
        return;
    end
    switched = 0;
    if size(mat1,1)~=size(mat2,1)
        if size(mat1,2)==size(mat2,2)
            mat1=mat1';
            mat2=mat2';
            switched = 1;
        else
            error('DTWaverage(): matrix sizes are not compatible');
        end
    end
    % initialise the coefficients weight by one if nothing specify or if cow==1
    if (nargin<3 || numel(cow)==0)
        % all the same weight
        cow=ones(1,size(mat1,1));
        % for 12PLPlogE tieke
        %cow=[0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.3 0.3 0.3 0.3 0.4];
    elseif cow==1
        cow=ones(1,size(mat1,1));
    end
    % if no sequence weight specify do not compute average sequence
    if (nargin<4)
        q = 0 ;
    end
    % default: do not do compression/expansion operation
    if (nargin<5)
        ce = 0 ;
    end
    %%%%%%%%%%%%%% CALL THE C FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% need to be compiled first: mex COPTIMFLAGS=-O3 [path]/DTWave.c
    [dist, mat4] = DTWave(mat1,mat2,cow,q,ce);
    % return the matrix in the same direction as input
    if switched==1, mat4 = mat4'; end
end

function [tdist] = TreeDist(tree,a,b)
    % tdist: find the number of nodes between nodes a and b in the Tree tree
    % intersection between the path back to the root for each node
    baca = zeros(1,size(tree,2)) ;
    bacb = zeros(1,size(tree,2)) ;
    idbaca = 1 ;
    idbacb = 1 ;
    while ( a~=0 )
        baca(idbaca) = a ; 
        a = tree(a);
        idbaca = idbaca + 1 ;
    end
    baca = baca(baca>0) ;
    while ( b~=0 )
        bacb(idbacb) = b ; 
        b = tree(b);
        idbacb = idbacb + 1 ;
    end
    bacb = bacb(bacb>0) ;
    if size(baca,2)>size(bacb,2)
        bacb = padarray(bacb, [0 size(baca,2)-size(bacb,2)], 'pre') ; 
        subbac = baca-bacb ;
        tdist = 2*sum(subbac~=0) - sum(bacb==0) - 1 ; % -1 because two leaves with same parent should have distance 1 (Samarasinghe2006)
    elseif size(bacb,2)>size(baca,2)
        baca = padarray(baca, [0 size(bacb,2)-size(baca,2)], 'pre') ;
        subbac = baca-bacb ;
        tdist = 2*sum(subbac~=0) - sum(baca==0) - 1 ;
    else
        subbac = baca-bacb ;
        tdist = 2*sum(subbac~=0) - 1 ;
    end
end

function [distSN, distNN, distSS] = distSNN(syllab,weight,reper,fieldnam,moded,verbose,ce)
% compute distance between a dataset and a tree
    if nargin<7
        ce = 0;
        if nargin<6
            if ~exist('verbose','var') || isempty(verbose), verbose=0 ; end
            if nargin<5
                moded='SNN' ;
                if nargin<4
                    fieldnam='cepvect' ;
                    if nargin<3
                        reper='' ;
                    end
                end
            end
        end
    end
    if verbose && usejava('desktop')  % shows progress
        h = waitbar(0,'computing distances...') ;
    end
    if numel(strfind(moded,'SN'))>0
        distSN=zeros(numel(syllab),numel(weight));
        nstep = numel(syllab) * numel(weight) ;
        istep = 1 ;
        if verbose==1, fprintf(1,'\ndistSN(%g): ',numel(syllab)) ; end
        for s=1:numel(syllab)
            if verbose==1, fprintf(1,'%g ',s ) ; end
            for w=2:numel(weight)
                if verbose==1 && w==numel(weight) && s==numel(syllab), fprintf(1,'\n') ; end
                if verbose && usejava('desktop')
                    waitbar( istep/nstep ) ;
                    istep=istep+1 ;
                end
                distSN(s,w) = DTWaverage(syllab(s).(fieldnam),weight{w}, 1, 0, ce) ;
            end
        end
        if numel(reper)>0
            save('-v7',fullfile(reper,'distSN.mat'),'distSN');
        end
    end
    if verbose && usejava('desktop')
        close(h) ;
        h = waitbar(0,'compute interneuron distances...') ;
    end
    if numel(strfind(moded,'NN'))>0
        distNN=zeros(numel(weight),numel(weight));
        nstep = ( numel(weight) * (numel(weight)-1) )/2 ;
        istep = 1 ;
        if verbose==1, fprintf(1,'\ndistNN(%g): ',numel(weight)) ; end
        for w=2:numel(weight)
            if verbose==1
                fprintf(1,'%g ',w) ;
                if w==numel(weight), fprintf(1,'\n') ; end
            end
            for w2=(w+1):numel(weight)
                if verbose && usejava('desktop')
                    waitbar( istep/nstep ) ;
                    istep=istep+1 ;
                end
                distNN(w,w2) = DTWaverage(weight{w},weight{w2}, 1, 0, ce) ;
                distNN(w2,w) = distNN(w,w2) ;
            end
        end
        if numel(reper)>0
            save('-v7',fullfile(reper,'distNN.mat'),'distNN');
        end
    end
    if nargout>3
        distSS=zeros(numel(syllab),numel(syllab));
        nstep = ( numel(syllab) * (numel(syllab)-1) )/2 ;
        istep = 1 ;
        for s=1:numel(syllab)
            for w=(s+1):numel(syllab)
                if verbose && usejava('desktop')
                    waitbar( istep/nstep ) ;
                    istep=istep+1 ;
                end
                distSS(s,w) = DTWaverage(syllab(s).(fieldnam),syllab(w).(fieldnam), 1, 0, ce) ;
                distSS(w,s) = distSS(s,w) ;
            end
        end
    end
    if verbose && usejava('desktop')
        close(h) ;
    end
end

function [ dico ] = AddExampleDico( dico, type_index, syllab_index )
% adds a new example in the dictionary
% returns the dictionary and the index of the newly added syllable
%
% type_index: syllable type to which an example has to be added
% syllab_index: syllable index in the array syllab of the new example
    dico(type_index).example = [dico(type_index).example syllab_index] ;
end

function [ CS, FS ] = ClustStat(dico,diconodes,distSN)
% computes several clustering statistics
% CS: closest sample of each centroid
% FS: furthest sample from each centroid
    % check if needed to compute distance matrices
    if numel(distSN)==0
        [distSN] = distSNN(syllab,weight,'',fieldnam,'SN',0) ;
    end
    CS = 0;
    FS = 0;
    nfullc = 0 ; % number of non empty clusters
    for nclust=1:numel(diconodes)
        % closest sample of each centroid
        CS = CS + min(distSN(:,diconodes(nclust)),[],1) ;
        % furthest sample from each centroid
        if numel(dico(nclust).example)>0
            nfullc = nfullc+1 ;
            FS = FS + max(distSN(dico(nclust).example,diconodes(nclust)),[],1) ;
        end
    end
    % compute the mean overall clusters
    CS = CS/nfullc ;
    FS = FS/nfullc ;
end

function [indice] = DaviesBouldin2(dico,nodes,distSN,distNN)
% Davies-Bouldin validity index calculation
% need a dictionary (dico), the nodes (it gives for each dico index the corresponding node in the distance matrices) 
% and the distance matrices between each node (cluster centroid, distNN) and each syllable versus each node (distSN)
    num = size(dico,2) ;
    DBtmp = 0 ;
    for clst=1:num % each cluster
        % retrieves the index of this class for the distance matrix
        c1 = nodes(clst) ;
        temporaire = zeros(1,size(distSN,1)) ; % for speeed, prealloc total number of syllables
        % Inner distance of class c
        Sc = sum(distSN(dico(clst).example,c1)) / numel(dico(clst).example) ;
        for clst2=1:num % against each class
            c2 = nodes(clst2) ;
            if ( c2==c1 )
                temporaire(c2) = 0 ;
                continue ; 
            end
            % inner distance of class c2
            Sc2 = sum(distSN(dico(clst2).example,c2)) / numel(dico(clst2).example) ;
            % distance between c and c2
            dcc2 = distNN( c1, c2 ) ; % if this value is very low (two cluster centroids are very close, then the DB indice will be very large
            if (dcc2==0)
                temporaire(c2) = 0 ; % the same syllable occurs twice in syllab array
            else
                temporaire(c2) = (Sc+Sc2)/dcc2 ;
            end
        end
        DBtmp = DBtmp + max(temporaire) ;
    end
    indice = DBtmp/num ;
end

function [syllab] = encodeWavlist( wavdir, analysisdir, coeffdir, normavec, transpos, savesig, varargin )
% encode the wav files in wavdir according to the coding parameters of the configuration files in analysisdir
% if there is more than one configuration file, each one is used and the vector sequences are merged
% the final vector sequence is saved in coeffdir with .vect extension instead of .wav
% need VOICEBOX for writing and reading HTK file format
    if ~exist('normavec','var') || isempty(normavec), normavec=[]; end
    if ~exist('coeffdir','var') || isempty(coeffdir), coeffdir=wavdir; end
    if ~exist('transpos','var') || isempty(transpos), transpos=0; end
    if ~exist('savesig','var')  || isempty(savesig), savesig=0; end
    if ~exist('analysisdir','var') || isempty(analysisdir), analysisdir=''; end
    [wintime1, hoptime1, numcep1, lifterexp1, sumpower1, preemph1, dither1, ...
     minfreq1, maxfreq1, nbands1, bwidth1, dcttype1, fbtype1, usecmp1, modelorder1, broaden1] = ...
             process_options(varargin, 'wintime', 0.025, 'hoptime', 0.010, ...
                  'numcep', 13, 'lifterexp', -22, 'sumpower', 1, 'preemph', 0.97, ...
                  'dither', 0, 'minfreq', 300, 'maxfreq', 20000, ...
                  'nbands', 26, 'bwidth', 1.0, 'dcttype', 3, ...
                  'fbtype', 'htkmel', 'usecmp', 0, 'modelorder', 0, 'broaden', 0);
    if nargout>0
        syllab = struct('seqvect',{},'filename','') ;
    end
    wavfiles = dir(fullfile(wavdir,'*.wav')) ;
    % check wav files found
    if isempty(wavfiles)==1 
        fprintf(1,'No wavfiles found in %s\n',wavdir); return
    end
    for sylf=1:numel(wavfiles)
        seqvect = [] ;
        sylfilename0 = fullfile(wavdir,wavfiles(sylf).name) ;
        syllab(sylf).filename = wavfiles(sylf).name ;
        if ~isempty(analysisdir) % use HTK:HCopy to encode files
            conffiles = dir(fullfile(analysisdir,'*')) ;
            for cff=1:numel(conffiles)
                if isdir(fullfile(analysisdir,conffiles(cff).name)), continue; end % avoid directories
                tmp=regexprep(wavfiles(sylf).name(end:-1:1),'vaw.','kth.','once');tmp=tmp(end:-1:1); % allows to replace just once, the last one
                sylfilename1 = fullfile(coeffdir,tmp) ;
                system(['HCopy -A -C ' fullfile(analysisdir,conffiles(cff).name) ' ' sylfilename0 ' ' sylfilename1]) ;
                [coeffseq,fp] = readhtk(sylfilename1) ;
                seqvect = [seqvect coeffseq]; % matrix transposed with readhtk (compared to readmfcc)
                delete(sylfilename1);
            end
        else % use calc_mfcc
            [seqvect] = calc_mfcc(sylfilename0, '', 'wintime', wintime1, 'hoptime', hoptime1,...
                                    'numcep', numcep1, 'lifterexp', lifterexp1, 'sumpower', sumpower1,...
                                    'preemph', preemph1, 'dither', dither1, 'minfreq', minfreq1,...
                                    'maxfreq', maxfreq1, 'nbands', nbands1, 'bwidth', bwidth1,... 
                                    'dcttype', dcttype1, 'fbtype', fbtype1, 'usecmp', usecmp1,...
                                    'modelorder', modelorder1, 'broaden', broaden1) ;
            fp = hoptime1;
        end
        % fprintf(1,'%i %i\n',size(seqvect,1),size(seqvect,2)) ;
        if numel(normavec)>0
            % normalise from0 to 1 according to a vector giving the structure, e.g. [1 12; 13 13; 14 25]
            % seqvect = norma_seqvect(seqvect',[1 12; 13 13; 14 25]) ; % NEED TO TRANSPOSE THE MATRIX
            seqvect = norma_seqvect(seqvect',normavec) ; % NEED TO TRANSPOSE THE MATRIX
            seqvect = seqvect' ;
        end
        tc = 9 ; % always use USER data format for HTK, required for reading the data later
        if nargout==0
            tmp=regexprep(wavfiles(sylf).name(end:-1:1),'vaw.','kth.','once');tmp=tmp(end:-1:1); % allows to replace just once, the last one
            writehtk(fullfile(coeffdir,tmp),seqvect,fp,tc); % NEED TO TRANSPOSE THE MATRIX BACK
        else
            if transpos==1
                syllab(sylf).seqvect = seqvect' ;
            else
                syllab(sylf).seqvect = seqvect ;
            end
        end
        if savesig==1
             [syllab(sylf).signal, syllab(sylf).Fs] = wavread(sylfilename0) ;
        end
    end
end

function [dico, diconodes, distSN, distNN, distSbmuT, DaBo] = etreedico(...
                                    tree,weight,syllab,lim,distSS,distSN,distNN,distSbmuT,typecl,fieldnam,ce)
% build a dico for a specific depth
% also compute the distance matrices between syllab and tree nodes if needed
% if typecl==1 then find all the nodes at the specified depth and then get find the best match
% if typecl==2 then the clustering is done going from the root down the leaves in the tree (like learning, stop at the specified depth)
% if typecl==3 then one cluster for each node and merge close nodes using a topology restriction (limdistNN2)
% if typecl==4 then one cluster for each node and classify distNN using dendrograms
% if typecl==5 then cluster the distSN(x,:) vectors using kmeans
% else a syllable is added in the cluster containing its best match over all nodes
% check also EtreePath() and bmudico()
% example [dico diconodes distSN distNN distSbmuT]=etreedico(tree,weight,[],[0 0.5 0.05 0 0],[],distSN,distNN,[],4);
% example [dico diconodes distSN distNN distSbmuT]=etreedico(tree,weight,[],[0 mean(reshape(distSN,1,[])) mean(reshape(distSN,1,[]))/2 0 0],[],distSN,distNN,[],4);
% example [dico6 diconodes6 distSN6 distNN6 distSbmuT6]=etreedico(tree,weight,[],[6 mean(reshape(distSN,1,[])) 0 0 0],[],distSN,distNN,[],1);
    % no need of syllab if method different from 2
    if typecl==2 && numel(syllab)==0
        error('need to provide a syllab array for this method');
    elseif (numel(distSN)==0 || numel(distNN)==0) && numel(syllab)==0
        error('need to provide a syllab array for this method');
    else
        syllab=[];
    end
    % default values for the limits
    if numel(lim)==0
        lim=[0 mean(reshape(distSN,1,[])) mean(reshape(distNN,1,[]))/100 0 0]; % /100 : no node merging
    end
    % lim array contains different limit and threshold for the classification
    limdep=lim(1);
    limdistSN=lim(2);
    limdistNN=lim(3);
    limdistNN2=lim(4);
    limdistIC=lim(5);
    if limdep==0 % no depth required, get the max depth of the tree
        [~,~,limdep] = treelayout(tree) ;
        clear x y s ; % we just need the depth of the tree
    end
    if limdistSN==0 % no limit for distance between sample and centroid
        limdistSN = max(reshape(distSN,1,[]))+1 ;
    end
    if nargin<11
        ce = 0;
        if nargin<10
            fieldnam='cepvect' ;
        end
    end
    % default value for Davies Bouldin indice
    DaBo = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTE DISTANCE MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if numel(distSN)==0 || numel(distNN)==0
        %%%%% NORMALIZATION (need to convert to cell to perform functions in one step) like it is done in EvolTreeDTW()
        % with this step the others fields of the structure "syllab" are deleted
        if ( max(reshape(syllab(1).(fieldnam),1,numel(syllab(1).(fieldnam))))~=1 || min(reshape(syllab(1).(fieldnam),1,numel(syllab(1).(fieldnam))))~=0 )
          syllab = cell2struct( cellfun(@(x) norma(x,0,1),{syllab.(fieldnam)}, 'UniformOutput', false), fieldnam, 1 )' ;
        end

        % each syllable against each node
        if numel(distSN)==0 && numel(distSbmuT)==0
            distSN=zeros(numel(syllab),numel(weight));
            fprintf(1,'\ndistSN(%g): ',numel(syllab)) ;
            for s=1:numel(syllab)
               fprintf(1,'%g ',s ) ;
               for w=2:numel(weight)
                   % distSN(s,w) = DTWoommen(syllab(s).(fieldnam),weight{w}) ;
                   distSN(s,w) = DTWaverage(syllab(s).(fieldnam),weight{w}, 1, 0, ce) ;
               end
            end
            % ANOTHER SYLLABLE PAIRWISE DISTANCE 
            % the pairwise syllable distance is the number of edges between BMU of each syllable in the etree
            distSbmuT=zeros(numel(syllab));
            fprintf(1,'\ndistSbmuT(%g): ',numel(syllab)) ;
            for s1=1:numel(syllab)
                fprintf(1,'%g ',s1 ) ;
                for s2=s1+1:numel(syllab)
                    distSbmuT(s1,s2)=TreeDist(tree,find(distSN(s1,2:end)==min(distSN(s1,2:end))),find(distSN(s2,2:end)==min(distSN(s2,2:end))));
                    distSbmuT(s2,s1)=distSbmuT(s1,s2);
                end
            end
        end
        % each node against each node
        if  numel(distNN)==0 && limdistNN~=0
            distNN=zeros(numel(weight),numel(weight));
            fprintf(1,'\ndistNN(%g): ',numel(weight)) ;
            for w=2:numel(weight)
               fprintf(1,'%g ',w) ;
               for w2=w+1:numel(weight)
                   % distNN(w,w2) = DTWoommen(weight{w},weight{w2}) ;
                   distNN(w,w2) = DTWaverage(weight{w},weight{w2}, 1, 0, ce) ;
                   distNN(w2,w) = distNN(w,w2) ;
               end
            end
        end
    end
    %%%%% DEPTH of each node
    depthn = zeros(1,numel(tree)) ;
    for n = 1:numel(tree)
        a=n ;
        while ( tree(a)~=0 ); depthn(n)=depthn(n)+1; a=tree(a); end
    end
    %%%%% FIND NODES FOR THIS DEPTH
    % find the nodes of the required depth and the ones of inferior depth
    % which are leaves
    leaves = setdiff(1:numel(tree),unique(tree)) ;
    % keep track of the correspondance between the dico index and the whole tree index
    % i.e which tree node corresponds to which dico index?
    % the cluster centres are the ones of depth==limdep and the ones of inferior depth but which are leaves
    diconodes = [ find(depthn==limdep) intersect(find(depthn<limdep),leaves) ] ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILL THE DICO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply a threshold to allow any adding in the dico
    if typecl==1
        %%%%% FIRST APPROACH: find best match for each syllable among all the leaves of specified depth
        dico = struct('example',cell(1,numel(diconodes))) ;
        for s=1:size(distSN(:,diconodes),1)
            if ( min(distSN(s,diconodes))<limdistSN )
                dico = AddExampleDico( dico, find(distSN(s,diconodes)==min(distSN(s,diconodes)),1), s ) ;
            end
        end
    elseif typecl==2
        %%%%% SECOND APPROACH: go through the tree for each syllable and stop when depth is reached (the node is in diconodes),
        %%%%% like the learning process do
        dico = struct('example',cell(1,numel(diconodes))) ;
        BMU = zeros(numel(syllab),limdep,2) ;
        for s=1:numel(syllab)
            BMU(s,1,1) = 1 ;
            dep = 1 ;
            while ( sum(tree==BMU(s,dep,1))>0 && numel(intersect(BMU(s,end,1),diconodes))==0 ) % this node is a father => keep on searching
                feuilles = find(tree==BMU(s,dep,1)) ;
                dep = dep + 1 ;
                scor = zeros(1,max(feuilles)) ;
                for L=feuilles
                    if numel(distSN)==0
                        % NORMALIZATION of the weights (need to deal with cells to faster calculations)
                        cellfun(@(x) norma(x,0,1),{weight{feuilles}}, 'UniformOutput', false) ;
                        [ scor(L) ] = DTWaverage( syllab(s).(fieldnam), weight{L}, 1, 0, ce ) ; % cepstrum coefficient alignment
                    else
                        scor(L) = distSN(s,L);
                    end
                end
                BMU(s,dep,2) = min(scor(feuilles)) ;
                BMU(s,dep,1) = find(scor==min(scor(feuilles)),1) ; % find best match (Best Matching Unit), minimum distance
            end
            if BMU(s,dep,2)<limdistSN
                dico = AddExampleDico( dico, find(diconodes==BMU(s,dep,1)), s ) ;
            end
        end
    elseif typecl==3 || typecl==4
        %%%%% THIRD/FOURTH APPROACH: one cluster for each node, clusters will be merged later
        diconodes = 1:numel(weight) ; % useful for Davies bouldin index calculation later on
        dico = struct('example',cell(1,numel(weight))) ;
        for s=1:size(distSN,1)
            if ( min(distSN(s,2:end))<limdistSN ) % avoid the root
                dicoidx = find(distSN(s,2:end)==min(distSN(s,2:end)),1) + 1 ; % add 1 because we avoided the root
                dico = AddExampleDico( dico, dicoidx, s ) ;
            end
        end
    elseif typecl==5
        %%%%% FIFTH APPROACH: cluster distSN vectors with kmeans
        diconodes = []; % diconodes is not required
        nc = 20 ; % number of clusters
        dico = struct('example',nc) ;
        idx=kmeans(distSN,nc,'distance','cosine'); % empirically, it is the distance giving the best silhouette (mean(silh)) results
        %[silh,h]=silhouette(distSN,idx);
        for dci = 1:nc
            dico(dci).example=find(idx==dci) ;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REMOVE EMPTY CLUSTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d=1 ;
    while ( d<=numel(dico) )
        if ( numel(dico(d).example)==0 ) 
            dico = DeleteTypeDico( dico, d ) ;
            if (d<numel(dico));  diconodes(d:end-1)=diconodes(d+1:end) ; end
            diconodes = diconodes(1:end-1) ;
        else d = d+1 ;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MERGE NODES IN THE ETREE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if typecl==4
        %%%%% FOURTH APPROACH: classify distNN using dendrograms
        distNNsq = squareform(distNN(diconodes,diconodes));
        Z = linkage(distNNsq,'average'); %figure;dendrogram(Z,0,'colorthreshold',limdistNN);
        T = cluster(Z,'cutoff',limdistNN,'criterion','distance');
        dc1 = 1 ;
        while dc1<size(diconodes,2)
            dc2 = dc1+1 ;
            while dc2<size(diconodes,2)
                if T(dc1)==T(dc2)
                    dico(dc1).example = [ dico(dc1).example dico(dc2).example ] ;
                    dico = DeleteTypeDico( dico, dc2 ) ;
                    if (dc2<numel(dico)); diconodes(dc2:end-1)=diconodes(dc2+1:end) ; end
                    diconodes = diconodes(1:end-1) ;
                    if (dc2<numel(T)); T(dc2:end-1)=T(dc2+1:end) ; end
                    T = T(1:end-1) ;
                    dc2 = dc2-1 ;
                end
                dc2 = dc2+1 ;
            end
            dc1 = dc1+1 ;
        end
    elseif typecl==3
        %%%%%%%%%% REFINE THE DICO BY MERGING CLOSED NODES WITH A TOPOLOGY RESTRICTION
        merged = 0 ; % to know if there had been a merging done or not
        if ( limdistNN2~=0 )
           tidx = 1 ;
            while tidx<numel(tree)
                sons = find(tree==tidx) ;
                dicidx1 = find(diconodes==tidx) ;
                if numel(dicidx1)==1 
                    if numel(dico(dicidx1).example)>0 % only merge is there is at least one example
                        for sidx = 1:numel(sons)
                            if distNN(sons(sidx),tidx)<limdistNN2
                                % merging the cluster in the dictionary
                                dicidx2 = find(diconodes==sons(sidx)) ;
                                if numel(dicidx2)==1
                                    if numel(dico(dicidx2).example)>0 % only merge is there is at least one example
                                        dico(dicidx1).example = [ dico(dicidx1).example dico(dicidx2).example ] ;
                                        dico = DeleteTypeDico( dico, dicidx2 ) ;
                                        if (dicidx2<numel(dico)); diconodes(dicidx2:end-1)=diconodes(dicidx2+1:end) ; end
                                        diconodes = diconodes(1:end-1) ;
                                        % update the tree
                                        [tree, weight, distNN] = deleteNode(tree, weight, distNN, sons(sidx), tidx, ce) ;
                                        merged = 1 ;
                                    end
                                end
                            end
                        end
                    end
                end
                if merged==0; tidx = tidx+1 ;
                else merged = 0 ; end
                %treeplot(tree) ; drawnow ;
            end
        end
    elseif typecl==1 || typecl==2 || typecl==3
        %%%%%%%%%% REFINE THE DICO BY MERGING CLOSED NODES
        d1 = 1 ;
        while d1<size(diconodes,2)
            d2 = d1+1 ;
            while d2<size(diconodes,2)
                if ( distNN(diconodes(d1),diconodes(d2))<limdistNN )
                    %fprintf(1,'merge %g and %g\n',d1,d2);
                    % merge d1 and d2
                    dico(d1).example = [ dico(d1).example dico(d2).example ] ;
                    % delete d2 in the dico
                    dico = DeleteTypeDico( dico, d2 ) ;
                    % delete in diconodes, just keep the first centroid
                    % !!!!!!!!!!!! THIS IS BIAISING THE CLUSTERING INDICES !!!!!!!!!!
                    % the nodes are not centroid anymore, centroid should be merged too
                    if (d2<numel(dico)); diconodes(d2:end-1)=diconodes(d2+1:end) ; end
                    diconodes = diconodes(1:end-1) ;
                    d2 = d2-1 ;
                end
                d2 = d2+1 ;
            end
            d1 = d1+1 ;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING STATISTICS COMPUTATION AND PRINTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% compute the INNER DISTANCE of each cluster and remove the cluster if too high
    if ( limdistSN~=0 && limdistIC~=0 )
        InnerDist = zeros(1,numel(dico)) ;
        for clst=1:numel(dico)
            cnd = diconodes(clst) ;
            InnerDist(clst) = sum(distSN(dico(clst).example,cnd)) / numel(dico(clst).example) ;
        end
        for d=(find(InnerDist>=limdistIC))
            dico = DeleteTypeDico( dico, d ) ;
            if (d<numel(dico)); diconodes(d:end-1)=diconodes(d+1:end) ; end
            diconodes = diconodes(1:end-1) ;
        end
    end
    %%%%%%%%%%%% compute CLUSTERING INDICE
    if (strcmp(language_of,'matlab')), fprintf(1,'%g clusters\n',sum(arrayfun(@(x) length(x.example),dico)>0)) ; end % NOT FOR OCTAVE
    if ( numel(distSN)>1 && numel(distNN)>1 && numel(diconodes)>0 )
        DaBo = DaviesBouldin2(dico,diconodes,distSN,distNN);
        fprintf(1,'Davies-Bouldin indice: %10.2f\n', DaBo) ;
    end
    if ( numel(distSS)>1 )
        res=zeros(1,sum(cellfun(@(x) numel(x),{dico.example}))) ;
        for d=1:numel(dico), res(dico(d).example)=d ; end
        fprintf(1,'Dunn indice: %g\n',DunnIndice(distSS,res)) ;
    end
    if ( numel(distSN)>1 )
        fprintf(1,'Classified syllable %g\n',sum(cellfun(@(x) numel(x),{dico.example}))/size(distSN,1)) ;
    end
end

function [tree, weight, distNN] = deleteNode(tree, weight, distNN, noded, noder, ce)
    % delete tree node noded (noder is his father and gets its sons) : set its sons to the noder and set no father to noded
    tree(tree==noded) = noder ;
    tree(noded) = 0 ;
    % compute the average of the two nodes weight and empty the matrix of noded
    [ ~, weight{noder} ] = DTWaverage( weight{noded}, weight{noder}, 1, 0.5, ce ) ;
    weight{noded} = [] ;
    % update distNN noder with his father
    if tree(noder)>0
        distNN(noder,tree(noder)) = DTWaverage( weight{tree(noder)}, weight{noder}, 1, 0, ce) ;
        distNN(tree(noder),noder) = distNN(noder,tree(noder)) ;
    end
    % update distNN noder with his sons
    fils = find(tree==noder) ;
    for fidx = 1:numel(fils)
        if fils(fidx)~=noded
            distNN(noder,fils(fidx)) = DTWaverage( weight{fils(fidx)}, weight{noder}, 1, 0, ce) ;
            distNN(fils(fidx),noder) = distNN(noder,fils(fidx)) ;
        end
    end
    % update distNN noder and noded which not have to be checked again
    distNN(noder,noded) = max(reshape(distNN,1,[])) ;
end

function [ dico ] = DeleteTypeDico( dico, t_index )
% deletes a type in the dictionary
%
% dico: dictionary to updated
% t_index: type to be deleted
    for n=t_index:numel(dico)-1
        dico(n) = dico(n+1) ;
    end
    dico = dico(1:end-1) ;
end

function [dunn] = DunnIndice (dmat,res)
% computes Dunn cluster validity indice
% dmat is a distance matrix
% res is a matrix defining the groups, size(res,1)=size(dmat,1) the value in res define the group which the sample belongs to
% res(x)==0 then x is not classified
    % find the biggest diameter (biggest intra-cluster distance)
    diam = 0 ;
    for clust=min(res(res>0)):max(res(res>0))
        if ( numel(find(res==clust))==0 ); continue; end
        % the distances between the element of this cluster
        dmatt = dmat(res==clust,res==clust) ;
        diamt = max(max(dmatt)) ;
        if (diamt(1)>diam) 
            diam = diamt(1) ; 
        end
    end
    % find the minimum of the minimum of inter-cluster distance divided by the biggest diameter
    dunn = [] ;
    for clust=min(res(res>0)):max(res(res>0))
        if ( numel(find(res==clust))==0 ); continue; end
        tmp1 = [] ;
        for clust2=min(res(res>0)):max(res(res>0))
            if ( numel(find(res==clust2))==0 || clust2==clust ); continue; end
            dmatt = dmat(res==clust,res==clust2) ;
            distmt = min(min(dmatt)) ;
            % compute the Dunn indice
            tmp1t = distmt(1) / diam ;
            tmp1 = [tmp1 tmp1t] ; 
        end
        dunn = [dunn min(tmp1)] ;
    end
    if (numel(dunn)>0)
        dunn = min(dunn) ;
    else
        dunn = 0 ;
    end 
end

function [ depth ] = getdepth(~,~,tree,weight,distSN,distNN,song,pathstr)
% get the mouse position in graph units
  mouse = get(gca,'currentpoint') ;
  % calculate the point nearest to the mouse click (depth)
  %[val, depth] = min( sqrt((get(objet,'xdata')-mouse(1,1)).^2 ...
  %   +(get(objet,'ydata')-mouse(1,2)).^2 ));
  depth = round(mouse(1,1)) ;
  fprintf(1,'->- depth selected: %i\n',round(depth)) ;
  % ask for confirmation
  if ~isempty(song)
      button = questdlg(['Create a library using depth ' num2str(depth) '?'],'library') ;
      if strcmp(button,'Yes')
         makelibrary(depth,tree,weight,distSN,distNN,song,pathstr);
      end
  end
end

function [ DB, CS, FS, NC ] = OptimDepth( tree, weight, syllab, distSN, distNN, fieldnam, plot3D )
% check different clustering statistics going through different depth of an evolving tree
% DB: Davies-Bouldin clustering index
% CS: closest sample of each centroid
% FS: furthest sample from each centroid
% NC: number of centroids for each depth
% depth: allow to chose the depth from the plot of DEPTH versus FS/NC (if called from )
    if ~exist('plot3D','var') || isempty(plot3D); plot3D=0; else plot3D=1; end
    % check if needed to compute distance matrices
    if numel(distSN)==0 && numel(distNN)==0
            [distSN, distNN] = distSNN(syllab,weight,'',fieldnam,'SNN',0) ;
    elseif numel(distSN)==0 && numel(distNN)>0
            [distSN] = distSNN(syllab,weight,'',fieldnam,'SN',0) ;
    elseif numel(distSN)>0 && numel(distNN)==0
            [distNN] = distSNN(syllab,weight,'',fieldnam,'NN',0) ;
    end
    [~,~,maxdep] = treelayout(tree) ;
    DB = zeros(1,maxdep) ;
    CS = zeros(1,maxdep) ;
    FS = zeros(1,maxdep) ;
    NC = zeros(1,maxdep) ;
    %%%%% DEPTH of each node
    depthn = zeros(1,numel(tree)) ;
    for n = 1:numel(tree)
        a=n ;
        while ( tree(a)~=0 ); depthn(n)=depthn(n)+1; a=tree(a); end
    end
    % record the leaves of the tree
    leaves = setdiff(1:numel(tree),unique(tree)) ;
    for dep=1:maxdep
        % the cluster centres are the ones of depth==limdep and the ones of inferior depth but which are leaves
        diconodes = [ find(depthn==dep) intersect(find(depthn<dep),leaves) ] ;
        % defines the cluster (dictionary)
        dico = struct('example',cell(1,numel(diconodes))) ;
        % redefine distSN and distNN with just the cluster centroids
        distSNc = distSN(:,[1 diconodes]) ;
        % build dico
        for s=1:size(distSNc,1)
            dicoidx = find(distSNc(s,2:end)==min(distSNc(s,2:end)),1) ; % add 1 because we avoided the root
            dico = AddExampleDico( dico, dicoidx, s ) ;
        end
        [ cs1, fs1 ] = ClustStat(dico,diconodes,distSN) ;
        CS(dep) = cs1 ;
        FS(dep) = fs1 ;
        NC(dep) = numel(diconodes) ;
        % Davies-Bouldin
        DB(dep) = DaviesBouldin2(dico,diconodes,distSN,distNN) ;
    end
    if plot3D==1  % plot FS 3D
        figure;
        plot3(1:maxdep,FS,NC);
        grid on;axis square;
        xlabel('Classification tree depth');
        ylabel('Mean largest intra cluster distance');
        zlabel('Number of clusters');
    end
    % plot 2D
    % a figure already exists with a varpanel?
    varpanel = getappdata(gcf,'varpanel') ;
    %set(varpanel,'Units','pixels');
    if numel(varpanel)>0
        axes('Parent',varpanel,'Units','pixels','Position',[100 60 400 400]) ;
        padre = varpanel ;
    else
        padre = gcf ;
    end
    set(gca,'Units','pixels');
    H1 = line(1:maxdep,DB,'Color','red','Marker','+','LineStyle','--');
    ax1 = gca;
    ax2 = axes('Parent',padre,'Units','pixels','Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','left','Color','none');
    line(1:maxdep,FS,'Parent',ax2,'Marker','.');
    ax3 = axes('Parent',padre,'Units','pixels','Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','right','Color','none');
    line(1:maxdep,NC,'Parent',ax3,'Marker','x');
    set(ax1,'YLim',[0.95 1.05].*(get(ax1,'YLim')),'Ytick',[] );
    set(ax2,'YLim',[0.95 1.05].*(get(ax2,'YLim')) );
    set(ax3,'YLim',[0.95 1.05].*(get(ax3,'YLim')) );
    set(ax1,'XLim',[0.95 1.05].*(get(ax1,'XLim')) );
    set(ax2,'XLim',[0.95 1.05].*(get(ax2,'XLim')) );
    set(ax3,'XLim',[0.95 1.05].*(get(ax3,'XLim')) );
    set(get(ax3,'YLabel'),'String','Number of clusters (x)','FontSize',14);
    set(get(ax2,'YLabel'),'String','Mean largest intra-cluster distance (.)','FontSize',14);
    set(get(ax1,'YLabel'),'String','Davies-Bouldin cluster index (--)');
    set(get(ax1,'XLabel'),'String','Classification tree depth','FontSize',14);
    legend(ax1,'show','Davies-Bouldin Index','Location','NorthWest'); legend(ax1,'Boxoff'); % don't know why need to do this in 2 commands
    % remove the Davies-Bouldin index line
    legend(ax1,'hide'); set(H1,'LineStyle','none','Marker','none');
end

function [language] = language_of()
    dir_list  = lower(path);
    n_octave  = strfind(dir_list, 'octave');
    n_octave  = size(n_octave, 2);
    n_matlab  = strfind(dir_list, 'matlab');
    n_matlab  = size(n_matlab, 2);
    if n_octave > n_matlab
            language = 'octave';
    elseif n_matlab > n_octave
            language = 'matlab';
    else
            error('language: can not determine if octave or matlab.');
    end
    return
end


end

        
