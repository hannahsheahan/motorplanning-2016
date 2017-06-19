%% Motor Planning, Not Execution, Separates Motor Memories
% Sheahan, Franklin & Wolpert. (2016), Neuron 92(4)
% - This is an analysis script for assessing learning in 4 main groups of
%   subjects (n=6 per group), and a control group (n=4).

% Group 1: full follow-through
% Group 2: no follow-through
% Group 3: execution only
% Group 4: planning only
% Group 5: control

% Notes:  script doesn't include kinematic stats
% Author: Hannah Sheahan, sheahan.hannah@gmail.com
% Date:   11/12/2015 (original)
%         28/05/2017 (cleaned up script)
%--------------------------------------------------------------------------

%% settings
clc; clear all; close all;
tic
FontSize = 14;
set(0, 'DefaultAxesFontSize',FontSize); clear FontSize;
set(0,'DefaultFigureWindowStyle','docked');
FLAGS.normMPE = 1;              % normalise MPE by group peak speed
FLAGS.plotextra = 0;            % plot additional figures (subject-by-subject analysis, hand paths, target appearance time distributions)

%% load subject data
ngroups = 5;
for group = 1:ngroups  % load in by group because 'planning only' and 'control' groups have different total no. trials
    
    switch group
        case 1
            load Fullfollowthrough_combined;
        case 2
            load Nofollowthrough_combined;
        case 3
            load Executiononly_combined;
        case 4
            load Planningonly_combined;
        case 5
            load Control_combined;
    end
    %% define general variables
    
    % Frame, trial, subject, group numbers
    trial  = D.ResampTrialNumber';
    ntrial = max(trial);                        %number of trials for each subject
    N      = D.Trials;                          %total number of trials across all subjects
    ndata  = size(D.FrameData.RobotPosition,2); %number of samples (data points) in longest trial
    nsubj  = max(D.subj);
    
    % Field types
    S.FIELD_NULL     = 0 +1;
    S.FIELD_VISCOUS  = 1 +1;
    S.FIELD_CHANNEL  = 2 +1;
    S.FIELD_PMOVE    = 3 +1;
    
    % State Names
    S.STATE_GO       = 5;
    S.STATE_MOVEWAIT = 6;
    S.STATE_MOVING0  = 7;
    S.STATE_MOVING1  = 8;
    S.STATE_FINISH   = 9;
    
    % binary field-type matrix
    fieldlist = unique(D.FieldType);
    field = zeros(N,length(fieldlist)+1);
    for i=1:length(fieldlist)
        FN = fieldlist(i)+1;
        field(:,FN) = D.FieldType==fieldlist(i);
    end
    
    % Experiment phases
    if group<4
        phaselist = cumsum([5 150 3]);
    else
        phaselist = cumsum([3 75 2]);
    end
    Baseline = D.PhaseIndex<=phaselist(1);  indB = find(Baseline==1);
    exposurephase = D.PhaseIndex>phaselist(1) & D.PhaseIndex<=phaselist(2);  indE = find(exposurephase==1);
    Post     = D.PhaseIndex>phaselist(2) & D.PhaseIndex<=phaselist(3);  indP = find(Post==1);
    
    clear i FN fieldlist phaselist Baseline Exposure Post
    usefulthings = ws2struct();             % save some general access variables
    frames = D.FrameData.Frames;
    state = D.FrameData.State;
    trialtime = D.FrameData.TrialTime;
    
    %----------------------------
    %% Kinematic variable processing
    
    % Reformat kinematic variables
    %   - pad kinematic variables with NaNs in the excess frames, ready for resampling
    posx   = squeeze(D.FrameData.RobotPosition(:,:,1));
    posy   = squeeze(D.FrameData.RobotPosition(:,:,2));
    velx   = squeeze(D.FrameData.RobotVelocity(:,:,1));
    vely   = squeeze(D.FrameData.RobotVelocity(:,:,2));
    forcex = squeeze(D.FrameData.RobotForces(:,:,1));
    forcey = squeeze(D.FrameData.RobotForces(:,:,2));
    
    ind1 = D.FrameData.Frames;
    ind2 = ndata-ind1;
    
    counter = repmat(D.fdir, [1,ntrial])';  % indicate which subjects are counterbalanced for field direction
    counter = counter(:);
    fdir = -sign(D.TargetAngle);            % use this to signal each force-field direction
    fdir = fdir.*counter;                   % multiply all *-1 if field association is counterbalanced
    usefulthings.fdir = fdir;
    
    %pad with Nans (since each trial lasts a different number of frames)
    padnan = @(x,ind1,ind2,N) cell2mat( arrayfun(@(k)([x(k,1:ind1(k)) NaN(1,ind2(k))]), 1:N, 'UniformOutput', false)');
    
    posx   = padnan(posx,ind1,ind2,N); 
    posy   = padnan(posy,ind1,ind2,N);
    velx   = padnan(velx,ind1,ind2,N);
    vely   = padnan(vely,ind1,ind2,N);
    forcex = padnan(forcex,ind1,ind2,N);
    forcey = padnan(forcey,ind1,ind2,N);
    
    posfullx = posx;  % store the full trajectory position data for plotting later
    posfully = posy;  % store the full trajectory position data for plotting later
    
    % Remove data after viapoint
    ind2 = findfirst(D.FrameData.State >= S.STATE_MOVING1,2);
    ind1 = ones(size(ind2));
    pad  = ndata -(ind2-ind1) -1;
    
    chopnan = @(x,ind1,ind2,N) cell2mat(arrayfun(@(k)([x(k,ind1(k):ind2(k)) NaN(1,pad(k))]), 1:N, 'UniformOutput', false)');
    
    posx = chopnan(posx,ind1,ind2,N);
    posy = chopnan(posy,ind1,ind2,N);
    velx = chopnan(velx,ind1,ind2,N);
    vely = chopnan(vely,ind1,ind2,N);
    forcex = chopnan(forcex,ind1,ind2,N);
    forcey = chopnan(forcey,ind1,ind2,N);
    pathlength = nansum(sqrt(diff(posy,1,2).^2+diff(posx,1,2).^2),2);
    
    %----------------------------
    %% Resample kinematic variables (to 1000 samples per trial)
    nsamp = 1000;
    
    % keep some full trajectory data
    original.state = state;
    original.trialtime = trialtime;
    original.posy = posfully;
    original.posx = posfullx;
    original.lngth = size(state,2);
    
    % Resample data for first section of movement to via point
    ind1 = findfirst(state >= S.STATE_MOVING0,2);
    ind2 = findfirst(state >= S.STATE_MOVING1,2);
    
    resample = @(x,frames) interp1(1:length(x),x,linspace(1,length(x),frames));
    kinematic_resample = @(x,ind1,ind2,N) cell2mat( arrayfun(@(k) resample(x(k,ind1(k):ind2(k)),nsamp), 1:N, 'UniformOutput', false)');
    
    F.state     = kinematic_resample(state,ind1,ind2,N);
    F.trialtime = kinematic_resample(trialtime,ind1,ind2,N);
    posx      = kinematic_resample(posx,ind1,ind2,N);
    posy      = kinematic_resample(posy,ind1,ind2,N);
    velx      = kinematic_resample(velx,ind1,ind2,N);
    vely      = kinematic_resample(vely,ind1,ind2,N);
    forcex    = kinematic_resample(forcex,ind1,ind2,N);
    forcey    = kinematic_resample(forcey,ind1,ind2,N);
    
    % save a copy of variables before rotating (for later trajectory plotting)
    F.rpos = [reshape(posx',[1,nsamp,N]);reshape(posy',[1,nsamp,N])];
    F.rvel = [reshape(velx',[1,nsamp,N]);reshape(vely',[1,nsamp,N])];
    F.rforce = [reshape(forcex',[1,nsamp,N]);reshape(forcey',[1,nsamp,N])];
    original.pos = [reshape(original.posx',[1,original.lngth,N]);reshape(original.posy',[1,original.lngth,N])];
    
    clear  posx posy velx vely forcex forcey ind1 ind2 state trialtime
    
    %% Rotate kinematics from each start position
    startangle = -D.HomeAngle.*(pi/180);
    startpos   = repmat(reshape(D.StartPosition(:,1:2)', [2,1,N]),[1,nsamp,1]);
    original.startpos = repmat(reshape(D.StartPosition(:,1:2)', [2,1,N]),[1,original.lngth,1]);
    
    rotate = @(x,theta,N) cell2mat(( arrayfun(@(k)(reshape([cos(theta(k)), -sin(theta(k)); sin(theta(k)), cos(theta(k))]' * squeeze(x(:,:,k)),[2,1,size(x,2)])), 1:N, 'UniformOutput', false)));
    
    original.pos = rotate((original.pos-original.startpos),startangle,N);
    pos   = rotate((F.rpos-startpos),startangle,N);
    vel   = rotate(F.rvel,startangle,N);
    force = rotate(F.rforce,startangle,N);
    
    F.posx = squeeze(pos(1,:,:));
    F.posy = squeeze(pos(2,:,:));
    F.velx = squeeze(vel(1,:,:));
    F.vely = squeeze(vel(2,:,:));
    F.forcex = squeeze(force(1,:,:));
    F.forcey = squeeze(force(2,:,:));
    original.posx = squeeze(original.pos(1,:,:));
    original.posy = squeeze(original.pos(2,:,:));
    
    %% Calculate MPE and adaptation
    
    % calculate mpe
    [test,ind] = max(abs(F.posx'));
    ind = sub2ind(size(F.posx),1:N,ind);
    mpe = F.posx(ind)'.*(fdir==1) - F.posx(ind)'.*(fdir==-1);     % accounts for counterbalance subjects and switching field directions
    
    % Calculate adaptation (% complete)
    ph1 = find(field(:,S.FIELD_VISCOUS)==1);
    fieldconstant = D.FieldConstants(ph1(1),1);                   % same for all subjects
    perffx = fieldconstant.*repmat(fdir,[1,1000]).*F.vely;        % ideal force. NB: we can only measure forces in the x-direction
    
    % find +-200ms window either side of time of max speed (NB: often takes entire length of trial, but limits very long trials)
    speed = sqrt((F.velx.^2)+(F.vely.^2));
    [~,vmaxi] = max(speed,[],2);                                  % index of max speed
    tsmpl = 0.2;
    lvmaxi = sub2ind(size(F.trialtime),[1:N]',vmaxi);
    tvmax  = F.trialtime(lvmaxi);
    t_low   = tvmax-tsmpl;                                        % limits of time interval
    t_up    = tvmax+tsmpl;
    for i=1:N
        smpl_L(i,1)  = length(find(F.trialtime(i,:) < tvmax(i) & F.trialtime(i,:) > t_low(i)));
        smpl_U(i,1)  = length(find(F.trialtime(i,:) > tvmax(i) & F.trialtime(i,:) < t_up(i)));
    end
    p1 = vmaxi-smpl_L;
    p2 = vmaxi+smpl_U;
    
    % regress actual force onto ideal force (zero intercept) for adaptation measure. (lambda/anonymous function here takes <half computing time of 'for' loop)
    adaptregress = @(x,y,p1,p2,N) cell2mat( arrayfun(@(k)(100.*regress(y(k,p1(k):p2(k))',x(k,p1(k):p2(k))') ), 1:N, 'UniformOutput', false));
    adaptation = adaptregress(perffx,F.forcex,p1,p2,N);
    
    %------------------------
    %% Calculate other useful checks e.g. planning time, misstrials, duration...
    
    % percentage of misstrials per subject
    for i=1:nsubj
        ind = [ntrial.*(i-1)+1,ntrial.*i];
        misstrialrate(i) = max(D.MissTrials(ind(1):ind(2))) ./ (max(D.MissTrials(ind(1):ind(2)))+ntrial);
    end
    
    % percentage of overshot trials per subject (for planning only group)
    %(NB: overshot = if subject travels 3cm+ (missthresh) in the rotated +y direction after the viapoint)
    missthresh = 3;
    targetdist = 12;
    overdistance = original.posy - repmat(targetdist,size(original.posy));
    
    ch = (original.state >= S.STATE_MOVING1);                     % after hitting the via-point
    overshot = sign(sum(ch & (overdistance >= missthresh),2));
    overshot(field(:,S.FIELD_CHANNEL)==1)=0;                      % ignore count on channel trials (subjects supposed to follow through)
    nincludedtrials = sum(field(1:ntrial,S.FIELD_CHANNEL)~=1);    
    for i=1:nsubj
        ind = [ntrial.*(i-1)+1,ntrial.*i];
        overshootrate(i) = sum(overshot(ind(1):ind(2))) ./ nincludedtrials;
    end
    
    % pre-movement time (from when target first appears)
    ph  = find((field(:,S.FIELD_CHANNEL)==1));
    appeardist = [0,0,10,6,0];
    appeardist = appeardist(group);
    
    % find duration and timing of target appearances
    iRX = sub2ind(size(original.state),[1:N],findfirst(original.state' >= S.STATE_MOVING0));
    iGo = sub2ind(size(original.state),[1:N],findfirst(original.state' >= S.STATE_GO));
    iFin = sub2ind(size(original.state),[1:N],findfirst((original.state >= S.STATE_GO),2,1,'last')');
    rxtime = original.trialtime(iRX);                                           % reaction time (incl. 300ms delay period)
    
    if group~=5
        ind = sub2ind(size(original.posy),[1:N],findfirst(original.posy' >= appeardist));
        appeartime = original.trialtime(ind) - rxtime;                          % time of target appear after movement onset
        appeartimeaftergo = original.trialtime(ind) - original.trialtime(iGo);  % time of target appear after go cue
        % nan appearance times for channel trials (will have appearance time of 0ms)
        appeartime(ph) = nan;
        appeartimeaftergo(ph) = nan;
        
    else
        appeartime = D.TargetAppearTime;
        % NB: if there's a misstrial, the actual appeartime measure becomes
        % inaccurate, so use the desired appear time for these trials (an approximation that doesn't consider ~11ms runtime delay; applies to a total n=4 trials across subjects)
        ind = findfirst(F.trialtime >= repmat(appeartime,[1,size(F.trialtime,2)]),2);
        appeartime(ind==0) = D.DesiredAppearTime(ind==0);                   % note: control onset times measured from start of delay period, not from movement onset
        appeartimeaftergo  = appeartime' - original.trialtime(iGo);         % will be negative when target appears before go cue
    end
    duration = original.trialtime(iFin) - original.trialtime(iGo);          % duration of trial
    
    % find dwell time at via-point
    viapos = repmat([0;targetdist],[1,size(original.posx)]);
    viaradius = 1.25;
    posfromvia = (original.pos - viapos);
    distvia = sqrt(squeeze(posfromvia(1,:,:)).^2 + squeeze(posfromvia(2,:,:)).^2);
    ind1 = sub2ind(size(original.trialtime),[1:N]',findfirst((distvia<viaradius),2));           % subject enters via point
    ind2 = sub2ind(size(original.trialtime),[1:N]',findfirst((distvia<viaradius),2,1,'last'));  % subject leaves via point
    dwell = original.trialtime(ind2)-original.trialtime(ind1);
    peakspeed  = speed(lvmaxi);
    
    % perpendicular error at halfway to via-point (6cm)
    halfwayy = 6;
    ind = sub2ind(size(F.posx),[1:N]',findfirst(F.posy>=halfwayy,2));
    halfwayperr = F.posx(ind).*fdir;
    
    timing.rxtime = rxtime;
    timing.appeartime = appeartime;
    timing.appeartimeaftergo = appeartimeaftergo;
    timing.peakspeed = peakspeed;
    timing.dwell = dwell;
    timing.duration = duration;
    timing.halfwayperr = halfwayperr;
    timing.overshootrate = overshootrate;
    timing.misstrialrate = misstrialrate;
    %------------------------
    
    %% save the data we care about according to group name
    MPE{group} = mpe;
    Adaptation{group} = adaptation;
    FrameData{group} = F;
    Timing{group} = timing;
    Experiment{group} = usefulthings;
    
    clearvars -except MPE FLAGS Adaptation FrameData Timing Experiment ngroups S
end

%-----------------------------------------------------------
%% Plot our data across subjects and groups
% (and do a little data smoothing)
fh = [];
nmpe_scalefactor = 52.4706;
colours = ColourSelect('ColourSelect.jpg',5);
% create plotting patch for exposure phase
P.Vertices = [ 4.5 -100; 155.5 -100; 155.5 100; 4.5 100];
P.Faces = [1 2 3 4];
P.FaceColor = [.3, .3, .3];
P.FaceAlpha = 0.08;
P.EdgeColor = 'white';
P.LineWidth = 0.1;

for group = 1:ngroups
    mpe = MPE{group};
    adaptation = Adaptation{group};
    F = FrameData{group};
    timing = Timing{group};
    usefulthings = Experiment{group};
    S = usefulthings.S;
    N = usefulthings.N;
    D = usefulthings.D;
    field  = usefulthings.field;
    ntrial = usefulthings.ntrial;
    nsubj  = usefulthings.nsubj;
    fdir   = usefulthings.fdir;
    indE   = usefulthings.indE;
    %-----
    
    %% Plot MPE for each subject, and group average
    statsblocks = 4;
    % normalise MPE by mean peak speed across all groups
    if FLAGS.normMPE
        mpe = mpe./timing.peakspeed;
        mpe = mpe.*nmpe_scalefactor;
    end
    pre=5;  exp=150; post=3;
    exposurephase = zeros(N,1); exposurephase(indE) = 1;
    indexp = find((field(:,S.FIELD_CHANNEL)~=1) & exposurephase);
    
    % smooth MPE by block (c=8)
    smooth_factor = 2;
    ind = find(field(:,S.FIELD_CHANNEL)~=1);                      % null and exposure trials give mpe measures
    c = 8;                                                        % non-channel trials per block
    mpeblock = mean(reshape(mpe(ind), c, length(mpe(ind))/c),1); 
    mpeblock = reshape(mpeblock,length(mpeblock)/nsubj,nsubj);
    if group==4                                                   % if it's the planning only condition, get rid of some extra data on the ends
        mpeblock = mpeblock(2:end-1,:);
    end
    if group<5
        stats.firstmpeblock{group} = mpeblock(pre+1:pre+statsblocks,:);
        stats.finalmpeblock{group} = mpeblock(pre+exp+1-statsblocks:pre+exp,:);
    end
    % smooth MPE by 2 blocks in exposure phase (c=16)
    c = c*smooth_factor;
    mpeblockexp = mean(reshape(mpe(indexp), c, length(mpe(indexp))/c),1);
    mpeblockexp = reshape(mpeblockexp,length(mpeblockexp)/nsubj,nsubj);
    
    % plot mpe per subject (smoothed per block only)
    if FLAGS.plotextra
        figure();
        for i=1:nsubj
            subplot(1,nsubj,i);
            patch(P); hold all;
            plot(linspace(1,ntrial,size(mpeblock,1)),mpeblock(:,i),'Color',colours(group,:)); hold on;
            plot([0 ntrial],[0,0],'k');
            axis([0 ntrial -1.5 4]);
            if i==1
                xlabel('Trials');
                ylabel('NMPE (cm)');
            end
        end
        strng = sprintf('NMPE per subject, group %d', group);
        suptitle(strng);
    end
    
    % average across subjects
    mpeblock_mean = mean(mpeblock,2);
    mpeblock_se = std(mpeblock,0,2) ./ sqrt(size(mpeblock,2));
    mpeblockexp_mean = mean(mpeblockexp,2);
    mpeblockexp_se = std(mpeblockexp,0,2) ./ sqrt(size(mpeblockexp,2));
    
    % plot average across all subjects in group
    if group~=5
        figure(1000);
        subplot(1,2,1);
        shadeplot(1:pre, mpeblock_mean(1:pre), mpeblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
        shadeplot(linspace(pre+1,exp+pre,length(mpeblockexp_mean)), mpeblockexp_mean, mpeblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
        shadeplot(exp+pre+1:exp+pre+post, mpeblock_mean(end-post+1:end), mpeblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % pre-exposure
        plot([0 pre+exp+post],[0 0],'k');
        ylabel('NMPE (cm)');
        xlabel('Block');
        axis([0 pre+exp+post+20 -1.5 4]);
        
        % plot the average after effects in a subpanel
        aftereffects_mean = mean(mean(mpeblock(end-post+1:end,:),1));
        aftereffects_se = std(mean(mpeblock(end-post+1:end,:),1)) ./ sqrt(nsubj);
        errorbar(pre+exp+post+3+group*3, aftereffects_mean,aftereffects_se,'k'); hold on;
        plot(pre+exp+post+3+group*3, aftereffects_mean, 'o','MarkerSize',7,'MarkerEdgeColor',colours(group,:), 'MarkerFaceColor',colours(group,:)); hold on;
    end
    %% Plot adaptation for each subject on the same figure (1 per group)
    
    % smooth adaptation by block
    if group<4
        ind = find(field(:,S.FIELD_CHANNEL)==1);                            % channel trials give adaptation measures
        indexp = find((field(:,S.FIELD_CHANNEL)==1) & exposurephase);
    elseif group==4
        inclchannel = (round(D.HomeAngle)==-180) | (round(D.HomeAngle==0));
        ind = find((field(:,S.FIELD_CHANNEL)==1) & inclchannel);         
        indexp = find((field(:,S.FIELD_CHANNEL)==1) & (exposurephase & inclchannel));
    elseif group==5
        inclchannel = (round(D.HomeAngle)==-180) | (round(D.HomeAngle==0));
        ind = find(((D.DesiredAppearTimeIndex==-1) & (field(:,S.FIELD_CHANNEL)==1)) & inclchannel);   % channel trials with targets shown from start of trial
        indexp = find(((D.DesiredAppearTimeIndex==-1) & (field(:,S.FIELD_CHANNEL)==1)) & (exposurephase & inclchannel));
    end
    
    c = 2;                                                                             % channel trials per block
    adaptationblock = mean(reshape(adaptation(ind), c, length(adaptation(ind))/c),1);  % smooth by block
    adaptationblock = reshape(adaptationblock,length(adaptationblock)/nsubj,nsubj);
    
    if group==4           % in the planning only condition, get rid of some extra data on the ends so group plotting is comparable
        adaptationblock = adaptationblock(2:end-1,:);
    end
    
    if group<5
        stats.firstadaptblock{group} = adaptationblock(pre+1:pre+statsblocks,:);
        stats.finaladaptblock{group} = adaptationblock(pre+exp+1-statsblocks:pre+exp,:);
    end
    % smooth adaptation by 2 blocks in exposure phase
    c = c*smooth_factor;
    adaptationblockexp = mean(reshape(adaptation(indexp), c, length(adaptation(indexp))/c),1);
    adaptationblockexp = reshape(adaptationblockexp,length(adaptationblockexp)/nsubj,nsubj);
    
    % plot by subject in group
    if FLAGS.plotextra
        figure();
        for i=1:nsubj
            subplot(1,nsubj,i);
            patch(P); hold all;
            plot(linspace(1,ntrial,size(adaptationblock,1)),adaptationblock(:,i),'Color',colours(group,:)); hold on;
            plot([0 ntrial],[0,0],'k');
            axis([0 ntrial -20 100]);
            if i==1
                xlabel('Trials');
                ylabel('Adapatation (%)');
            end
        end
        strng = sprintf('Adaptation per subject, group %d', group);
        suptitle(strng);
    end
    
    % average across subjects
    adaptblock_mean = mean(adaptationblock,2);
    adaptblock_se = std(adaptationblock,0,2) ./ sqrt(size(adaptationblock,2));
    adaptblockexp_mean = mean(adaptationblockexp,2);
    adaptblockexp_se = std(adaptationblockexp,0,2) ./ sqrt(size(adaptationblockexp,2));
    
    % plot average across all subjects in group
    if group~=5
        figure(1000);
        subplot(1,2,2);
        shadeplot(1:pre, adaptblock_mean(1:pre), adaptblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
        shadeplot(linspace(pre+1,exp+pre,length(adaptblockexp_mean)), adaptblockexp_mean, adaptblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
        h = shadeplot(exp+pre+1:exp+pre+post, adaptblock_mean(end-post+1:end), adaptblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % post-exposure
        plot([0 pre+exp+post],[0 0],'k');
        ylabel('NMPE (cm)');
        xlabel('Block');
        axis([0 pre+exp+post -20 60]);
        fh = [fh,h];
    end
    
    %% Plot some special figures for the control experiment
    % NMPE, adaptation, + moving average of adaptation with onset time
    if group==5
        figure(1100);
        subplot(1,3,1);   %mpe
        pre=5;  exp=150; post=3;
        shadeplot(1:pre, mpeblock_mean(1:pre), mpeblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
        shadeplot(linspace(pre+1,exp+pre,length(mpeblockexp_mean)), mpeblockexp_mean, mpeblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
        shadeplot(exp+pre+1:exp+pre+post, mpeblock_mean(end-post+1:end), mpeblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % pre-exposure
        plot([0 pre+exp+post],[0 0],'k'); hold on;
        ylabel('NMPE (cm)');
        xlabel('Block');
        axis([0 pre+exp+post -1.5 4]);
        
        subplot(1,3,2);   % adaptation on channel trials appearing from start
        shadeplot(1:pre, adaptblock_mean(1:pre), adaptblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
        shadeplot(linspace(pre+1,exp+pre,length(adaptblockexp_mean)), adaptblockexp_mean, adaptblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
        shadeplot(exp+pre+1:exp+pre+post, adaptblock_mean(end-post+1:end), adaptblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % pre-exposure
        plot([0 pre+exp+post],[0 0],'k'); hold on;
        ylabel('Adaptation (%)');
        xlabel('Block');
        axis([0 pre+exp+post -20 60]);
        
        subplot(1,3,3);   % adaptation on random appearing channel trials vs. target appearance time
        secondhalfexp = (D.PhaseIndex > size(mpeblockexp,1)/2+3) & (exposurephase);
        ind = find(((D.DesiredAppearTimeIndex~=-1) & (field(:,S.FIELD_CHANNEL)==1)) & (secondhalfexp));           % random appear time channel trials in exposure phase
        adaptrandappear = adaptation(ind);
        randappeartime  = timing.appeartime(ind);
        tlow = 0.1;
        thigh = 0.6;
        ts=linspace(tlow,thigh,100);
        twindow = 0.15;
        
        adaptrandappear = reshape(adaptrandappear,length(adaptrandappear)/nsubj,nsubj);
        randappeartime = reshape(randappeartime,length(randappeartime)/nsubj,nsubj);
        
        for sub=1:nsubj
            x = adaptrandappear(:,sub);
            t = randappeartime(:,sub);
            tmp = rangemean(ts,t,x,twindow,0);
            mx(sub,:)=tmp.mx;
            my(sub,:)=tmp.my;
        end
        
        shadeplot(mean(mx),mean(my),stderr(my),'-',colours(group,:),0.2); hold on;
        plot([tlow thigh],[0 0],'k'); hold on;
        plot([0.3 0.3],[-20 50],'--k'); hold on;
        text(0.24,55,'Cue to move')
        ylabel('Adaptation (%)');
        xlabel('Target appearance time (s)');
        axis([tlow thigh -20 60]);
        
        if FLAGS.plotextra
            % Plot the distribution of onset times for all trials
            figure();
            subplot(1,2,1);
            h = histogram(1000*timing.appeartime, 30);
            h.FaceColor = colours(group,:);
            xlabel('Appear time (ms)');
            title('All trials');
            
            subplot(1,2,2);
            % NB: actual appear times will not be perfectly uniformly sampled across the 0 -
            % 700ms window. This is because if subject moved 10cm into reach before
            % target was scheduled to appear, target would appear anyway to ensure
            % it was always visible before reaching the via-point. So some of the
            % larger scheduled 'appearance times' (close to 700ms) were pushed
            % slightly earlier in practice.
            ph = find(D.DesiredAppearTime ~= 0);
            h = histogram(1000*timing.appeartime(ph),30);
            h.FaceColor = colours(group,:);
            xlabel('Appear time (ms)');
            title('Random-onset trials only');
            suptitle('Control experiment: distribution of target appearance times');
        end
    end
    
    %% Plot the trajectories for each group, at different experiment phases
    % Note: this takes a long time to execute
    if FLAGS.plotextra
        if group~=5                     % ignore trajectories in the control group
            if group~=4
                pre=[4,5];
                earlyexp=[1,2];
                lateexp=[149,150];
                post=[6,7];
            else
                pre=[5,6];
                earlyexp=[1,2];
                lateexp=[149,150];
                post=[7,8];
            end
            
            % for each starting position, extract left and right trajectories for 2-block epochs in each phase
            figure();
            homeangles = unique(round(D.HomeAngle));
            for angle=1:length(homeangles)
                i = homeangles(angle);
                
                n1 = find(((field(:,S.FIELD_NULL)==1) & (fdir==1))  & (round(D.HomeAngle)==i)); %leftward target    0 degrees
                n2 = find(((field(:,S.FIELD_NULL)==1) & (fdir==-1))  & (round(D.HomeAngle)==i)); %rightward target
                n1 = reshape(n1,length(n1)/nsubj,nsubj);
                n2 = reshape(n2,length(n2)/nsubj,nsubj);
                pren1 = n1(pre,:);
                postn1 = n1(post,:);    
                pren2 = n2(pre,:);
                postn2 = n2(post,:);
                
                e1 = find(((field(:,S.FIELD_VISCOUS)==1) & (fdir==1)) & (round(D.HomeAngle)==i)); %leftward target    0 degrees
                e2 = find(((field(:,S.FIELD_VISCOUS)==1) & (fdir==-1)) & (round(D.HomeAngle)==i)); %rightward target
                e1 = reshape(e1,length(e1)/nsubj,nsubj);
                e2 = reshape(e2,length(e2)/nsubj,nsubj);
                earlye1 = e1(earlyexp,:);
                latee1 = e1(lateexp,:);
                earlye2 = e2(earlyexp,:);
                latee2 = e2(lateexp,:);
                
                % average the trajectories across blocks and subjects
                x = mean(reshape(squeeze(F.rpos(1,:,pren1(:)))',[2,nsubj,1000]),1);   % pre-exposure, left
                y = mean(reshape(squeeze(F.rpos(2,:,pren1(:)))',[2,nsubj,1000]),1);
                trajectory_x.pre1(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.pre1(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                x = mean(reshape(squeeze(F.rpos(1,:,pren2(:)))',[2,nsubj,1000]),1);   % pre-exposure, right
                y = mean(reshape(squeeze(F.rpos(2,:,pren2(:)))',[2,nsubj,1000]),1);
                trajectory_x.pre2(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.pre2(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                x = mean(reshape(squeeze(F.rpos(1,:,earlye1(:)))',[2,nsubj,1000]),1);   % early-exposure, left
                y = mean(reshape(squeeze(F.rpos(2,:,earlye1(:)))',[2,nsubj,1000]),1);
                trajectory_x.early1(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.early1(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                x = mean(reshape(squeeze(F.rpos(1,:,earlye2(:)))',[2,nsubj,1000]),1);   % early-exposure, right
                y = mean(reshape(squeeze(F.rpos(2,:,earlye2(:)))',[2,nsubj,1000]),1);
                trajectory_x.early2(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.early2(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                x = mean(reshape(squeeze(F.rpos(1,:,latee1(:)))',[2,nsubj,1000]),1);   % late-exposure, left
                y = mean(reshape(squeeze(F.rpos(2,:,latee1(:)))',[2,nsubj,1000]),1);
                trajectory_x.late1(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.late1(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                x = mean(reshape(squeeze(F.rpos(1,:,latee2(:)))',[2,nsubj,1000]),1);   % late-exposure, right
                y = mean(reshape(squeeze(F.rpos(2,:,latee2(:)))',[2,nsubj,1000]),1);
                trajectory_x.late2(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.late2(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                x = mean(reshape(squeeze(F.rpos(1,:,postn1(:)))',[2,nsubj,1000]),1);   % post-exposure, left
                y = mean(reshape(squeeze(F.rpos(2,:,postn1(:)))',[2,nsubj,1000]),1);
                trajectory_x.post1(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.post1(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                x = mean(reshape(squeeze(F.rpos(1,:,postn2(:)))',[2,nsubj,1000]),1);   % post-exposure, right
                y = mean(reshape(squeeze(F.rpos(2,:,postn2(:)))',[2,nsubj,1000]),1);
                trajectory_x.post2(:,:) = [squeeze(mean(x)),squeeze(std(x)./sqrt(nsubj))];
                trajectory_y.post2(:,:) = [squeeze(mean(y)),squeeze(std(y)./sqrt(nsubj))];
                
                % plot the group average for different points in the experiment
                subplot(1,4,1);   % pre-exposure
                % matlab downloaded 'circles' function to draw the circles
                % in the background
                
                shadedTrajectory(trajectory_x.pre1', trajectory_y.pre1', [0.8 0.2 0.2], 0.1); hold on;
                shadedTrajectory(trajectory_x.pre2', trajectory_y.pre2', [0.2 0.2 0.8], 0.1); hold on;
                title('Pre-exposure');
                xlabel('Position (cm)');
                ylabel('Position (cm)');
                axis equal;
                
                subplot(1,4,2);   % early-exposure
                shadedTrajectory(trajectory_x.early1', trajectory_y.early1', [0.8 0.2 0.2], 0.1); hold on;
                shadedTrajectory(trajectory_x.early2', trajectory_y.early2', [0.2 0.2 0.8], 0.1); hold on;
                title('Early exposure');
                axis equal;
                
                subplot(1,4,3);   % late-exposure
                shadedTrajectory(trajectory_x.late1', trajectory_y.late1', [0.8 0.2 0.2], 0.1); hold on;
                shadedTrajectory(trajectory_x.late2', trajectory_y.late2', [0.2 0.2 0.8], 0.1); hold on;
                title('Late exposure');
                axis equal;
                
                subplot(1,4,4);   % post-exposure
                shadedTrajectory(trajectory_x.post1', trajectory_y.post1', [0.8 0.2 0.2], 0.1); hold on;
                shadedTrajectory(trajectory_x.post2', trajectory_y.post2', [0.2 0.2 0.8], 0.1); hold on;
                title('Post-exposure');
                axis equal;
                clear trajectory_x trajectory_y
            end
            strng = sprintf('Average hand paths, group %d ',group);
            suptitle(strng);
        end
    end
    
    clearvars -except Adaptation P MPE stats ngroups nsubj ntrial N S FLAGS Timing Experiment colours FrameData nmpe_scalefactor fh A M
end
%%
figure(1000);                       % plot some background shaded areas
subplot(1,2,1)
pre=5;  exp=150; post=3;
Pg = P;
Pg.Vertices = [ pre+.5 -100; exp+pre+.5 -100; exp+pre+.5 100; pre+.5 100];
patch(Pg); hold on;
subplot(1,2,2)
patch(P); hold on;
legend([fh(1)',fh(2)',fh(3)',fh(4)'],{'Full follow-through','No follow-through','Execution only','Planning only'});

figure(1100);
subplot(1,3,1);
patch(Pg); hold on;
subplot(1,3,2);
patch(P); hold on;

clear pre exp post fh P Pg ntrial N
%------------------------------
%% Statistical Analysis
% Used first 4 blocks (= 64 exposure trials, 8 channel trials) and last 4
% blocks of data for epochs, and do repeated measures anova.

for group=1:4                    % within-group learning

    [~, statstable.mpe{group}] = anova_rm([squeeze(mean(stats.firstmpeblock{group}))', squeeze(mean(stats.finalmpeblock{group}))'],'off');
    [~, statstable.adaptation{group}] = anova_rm([squeeze(mean(stats.firstadaptblock{group}))', squeeze(mean(stats.finaladaptblock{group}))'],'off');
    
    if group == 1
        mpestats.group1 = [squeeze(mean(stats.firstmpeblock{group}))', squeeze(mean(stats.finalmpeblock{group}))'];
        adaptstats.group1 = [squeeze(mean(stats.firstadaptblock{group}))', squeeze(mean(stats.finaladaptblock{group}))'];
    elseif group == 4
        mpestats.group4 = [squeeze(mean(stats.firstmpeblock{group}))', squeeze(mean(stats.finalmpeblock{group}))'];
        adaptstats.group4 = [squeeze(mean(stats.firstadaptblock{group}))', squeeze(mean(stats.finaladaptblock{group}))'];
    end
end

% compare learning between full follow-through and planning only groups
[~, statstable.comparegroup1_4.mpe] = anova_rm({mpestats.group1, mpestats.group4},'off');  
[~, statstable.comparegroup1_4.adaptation] = anova_rm({adaptstats.group1, adaptstats.group4},'off'); 

% Can do some additional kinematic stats analysis, if desired.
toc
