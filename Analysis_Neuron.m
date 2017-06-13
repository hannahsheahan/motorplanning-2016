
%% Motor Planning, Not Execution, Separates Motor Memories
% Sheahan, Franklin & Wolpert. (2016), Neuron
% - Analysis script for assessing learning in 4 groups of subjects (n=6 per group)

% Group 1: full follow-through
% Group 2: no follow-through
% Group 3: execution only
% Group 4: planning only

% Author: Hannah Sheahan, sheahan.hannah@gmail.com
% Date:  11/12/2015 (original)
%        28/05/2017 (cleaned up script)
%--------------------------------------------------------------------------

%% settings
clc; clear all; close all;
tic
FontSize = 14;
set(0, 'DefaultAxesFontSize',FontSize); clear FontSize;
set(0,'DefaultFigureWindowStyle','docked');

%% load subject data
ngroups = 4;
%ngroups = 1;
for group = 1:ngroups  % load in group by group because 'planning only' has different trial numbers (phasessize: 6 150 4 vs 5 150 3)
    
    switch group
        case 1
            load Fullfollowthrough_combined;
        case 2
            load Nofollowthrough_combined;
        case 3
            load Executiononly_combined;
        case 4
            load Planningonly_combined;
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
    
    % state Names
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
    if group < 4
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
    
    % Reformat kinematic variables ***Turn into function later padNan
    %   - padd kinematic variables with NaNs in the excess frames, ready for resampling
    %   - can be used to flip x/y kinematics by the field direction for counterbalance subjects)
    posx   = squeeze(D.FrameData.RobotPosition(:,:,1));
    posy   = squeeze(D.FrameData.RobotPosition(:,:,2));
    velx   = squeeze(D.FrameData.RobotVelocity(:,:,1));
    vely   = squeeze(D.FrameData.RobotVelocity(:,:,2));
    forcex = squeeze(D.FrameData.RobotForces(:,:,1));
    forcey = squeeze(D.FrameData.RobotForces(:,:,2));
    
    %pad with Nans (since each trial lasts a different number of frames)
    ind1 = D.FrameData.Frames;
    ind2 = ndata-ind1;
    
    flipdir = ones(N,1);                    % flipping before you rotate will mess up so don't do it
    flipdir = flipdir(:);
    counter = repmat(D.fdir, [1,ntrial])';  % indicate which subjects are counterbalanced for field direction
    counter = counter(:);
    fdir = -sign(D.TargetAngle);            % use this to signal each force-field direction (useful for MPE and adaptation calculations)
    fdir = fdir.*counter;                   % multiply all *-1 if field association is counterbalanced
    usefulthings.fdir = fdir;
    
    %This multiplies the x/y variable by flipdir (NB: we have to rotate later, so don't flip yet!) and pads the unused samples with nans from ind1 to ind2
    padnan = @(x,r,ind1,ind2,N) cell2mat( arrayfun(@(k)(r(k)*[x(k,1:ind1(k)) NaN(1,ind2(k))]), 1:N, 'UniformOutput', false)');
    
    posx   = padnan(posx,flipdir,ind1,ind2,N);
    posy   = padnan(posy,flipdir,ind1,ind2,N);
    velx   = padnan(velx,flipdir,ind1,ind2,N);
    vely   = padnan(vely,flipdir,ind1,ind2,N);
    forcex = padnan(forcex,flipdir,ind1,ind2,N);
    forcey = padnan(forcey,flipdir,ind1,ind2,N);
    
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
    mpe = F.posx(ind)'.*(fdir==1) - F.posx(ind)'.*(fdir==-1);                           % accounts for counterbalance subjects and switching field directions
    
    %mpe = max(F.posx')'.*(fdir==1) - min(F.posx')'.*(fdir==-1);                           % accounts for counterbalance subjects and switching field directions
    
    % Calculate adaptation (% complete)
    ph1 = find(field(:,S.FIELD_VISCOUS)==1);
    fieldconstant = D.FieldConstants(ph1(1),1);          % same for all subjects
    perffx = fieldconstant.*repmat(fdir,[1,1000]).*F.vely;  % ideal force. NB: we can only measure forces in the x-direction
    
    % find ?200ms window either side of time of max speed (NB: often takes entire length of trial, but limits very long trials)
    speed = sqrt((F.velx.^2)+(F.vely.^2));
    [~,vmaxi] = max(speed,[],2);                     % index of max speed
    tsmpl = 0.2;
    lvmaxi = sub2ind(size(F.trialtime),[1:N]',vmaxi);
    tvmax  = F.trialtime(lvmaxi);
    t_low   = tvmax-tsmpl;                          % limits of time interval
    t_up    = tvmax+tsmpl;
    for i=1:N
        smpl_L(i,1)  = length(find(F.trialtime(i,:) < tvmax(i) & F.trialtime(i,:) > t_low(i)));
        smpl_U(i,1)  = length(find(F.trialtime(i,:) > tvmax(i) & F.trialtime(i,:) < t_up(i)));
    end
    p1 = vmaxi-smpl_L;
    p2 = vmaxi+smpl_U;
    
    % regress actual force onto ideal force (zero intercept) for adaptation measure. (NB: lambda/anonymous function here takes <half computing time of 'for' loop)
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
    
    ch = (original.state >= S.STATE_MOVING1);       % after hitting the via-point
    overshot = sign(sum(ch & (overdistance >= missthresh),2));
    overshot(field(:,S.FIELD_CHANNEL)==1)=0;        % ignore count on channel trials (subjects supposed to follow through)
    nincludedtrials = sum(field(1:ntrial,S.FIELD_CHANNEL)~=1);    % same for each subject in a group
    for i=1:nsubj
        ind = [ntrial.*(i-1)+1,ntrial.*i];
        overshootrate(i) = sum(overshot(ind(1):ind(2))) ./ nincludedtrials;
    end
    
    % pre-movement time (from when target first appears)
    ph  = find((field(:,S.FIELD_CHANNEL)==1));
    appeardist = [0,0,10,6];
    appeardist = appeardist(group);
    
    % find duration and timing of target appearances
    ind = sub2ind(size(original.posy),[1:N],findfirst(original.posy' >= appeardist));
    iRX = sub2ind(size(original.state),[1:N],findfirst(original.state' >= S.STATE_MOVING0));
    iGo = sub2ind(size(original.state),[1:N],findfirst(original.state' >= S.STATE_GO));
    iFin = sub2ind(size(original.state),[1:N],findfirst((original.state >= S.STATE_GO),2,1,'last')');
    rxtime = original.trialtime(iRX);                 % reaction time (incl. 300ms delay period)
    appeartime = original.trialtime(ind) - rxtime;    % time of target appear after movement onset
    appeartimeaftergo = original.trialtime(ind) - original.trialtime(iGo);  % time of target appear after go cue
    duration = original.trialtime(iFin) - original.trialtime(iGo);          % duration of trial
    
    % nan appearance times and rxtimes for channel trials
    rxtime(ph) = nan;
    appeartime(ph) = nan;
    appeartimeaftergo(ph) = nan;
    
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
    timing.misstrialrate = misstrialrate;  % the miss trials and overshoot rate are kind of timing-dependent, so put them here.
    %------------------------
    
    %% save the data we care about according to group name
    MPE(group,:) = mpe;
    Adaptation(group,:) = adaptation;
    FrameData{group} = F;
    Timing{group} = timing;
    Experiment{group} = usefulthings;
    
    clearvars -except MPE Adaptation FrameData Timing Experiment ngroups S
end

%% Now all that's left to do is to plot
% (and a little data smoothing)
fh = [];
FLAGS.normMPE = 1;
FLAGS.plotextra = 1;
nmpe_scalefactor = 52.4706;
colours = ColourSelect('ColourSelect.jpg',4);
for group = 1:ngroups
    mpe = MPE(group,:)';
    adaptation = Adaptation(group,:)';
    F = FrameData{group};
    timing = Timing{group};
    usefulthings = Experiment{group};
    v2struct(usefulthings);
    %-----
    
    % create plotting patch for exposure phase
    P.Vertices = [ 4.5 -100; 155.5 -100; 155.5 100; 4.5 100];
    P.Faces = [1 2 3 4];
    P.FaceColor = [.3, .3, .3];
    P.FaceAlpha = 0.08;
    P.EdgeColor = 'white';
    P.LineWidth = 0.1;
    
    
    %% Plot MPE for each subject, and group average
    
    % normalise MPE by mean peak speed across all groups
    if FLAGS.normMPE
        mpe = mpe./timing.peakspeed;
        mpe = mpe.*nmpe_scalefactor;
    end
    
    exposurephase = zeros(N,1); exposurephase(indE) = 1;
    indexp = find((field(:,S.FIELD_CHANNEL)~=1) & exposurephase);
    
    % smooth MPE by block
    smooth_factor = 2;
    ind = find(field(:,S.FIELD_CHANNEL)~=1);                      % null and exposure trials give mpe measures
    c = 8;                                                        % non-channel trials per block
    mpeblock = mean(reshape(mpe(ind), c, length(mpe(ind))/c),1);  % smooth by block (c=8)
    mpeblock = reshape(mpeblock,length(mpeblock)/nsubj,nsubj);
    
    % smooth MPE by 2 blocks
    c = c*smooth_factor;
    mpeblockexp = mean(reshape(mpe(indexp), c, length(mpe(indexp))/c),1); % smooth every 2 blocks (c=16) in exposure phase
    mpeblockexp = reshape(mpeblockexp,length(mpeblockexp)/nsubj,nsubj);
    
    % plot mpe per subject (smoothed per block only)
    figure();
    for i=1:nsubj
        subplot(1,6,i);
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
    
    if group==4           % if it's the planning only condition, get rid of some extra data on the ends
        mpeblock = mpeblock(2:end-1,:);
    end
    
    % average across subjects
    mpeblock_mean = mean(mpeblock,2);
    mpeblock_se = std(mpeblock,0,2) ./ sqrt(size(mpeblock,2));
    mpeblockexp_mean = mean(mpeblockexp,2);
    mpeblockexp_se = std(mpeblockexp,0,2) ./ sqrt(size(mpeblockexp,2));
    
    % plot average across all subjects in group
    figure(1000);
    subplot(2,1,1);
    pre=5;  exp=150; post=3;
    shadeplot(1:pre, mpeblock_mean(1:pre), mpeblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
    shadeplot(linspace(pre+1,exp+pre,length(mpeblockexp_mean)), mpeblockexp_mean, mpeblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    shadeplot(exp+pre+1:exp+pre+post, mpeblock_mean(end-post+1:end), mpeblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % pre-exposure
    plot([0 pre+exp+post],[0 0],'k');
    ylabel('NMPE (cm)');
    xlabel('Block');
    Pg = P;
    Pg.Vertices = [ pre+.5 -100; exp+pre+.5 -100; exp+pre+.5 100; pre+.5 100];
    if group ==1
        patch(Pg); hold on;
    end
    axis([0 pre+exp+post -1.5 4]);
    
    
    % add in a close-up of the NMPE at the end of the MPE plot
    
    
    
    %% Plot adaptation for each subject on the same figure (1 per group)
    
    % smooth adaptation by block
    indexp = find((field(:,S.FIELD_CHANNEL)==1) & exposurephase);
    ind = find(field(:,S.FIELD_CHANNEL)==1);                      % null and exposure trials give mpe measures
    if group <4
        c = 2;      % non-channel trials per block
    else
        c = 4;
    end
    adaptationblock = mean(reshape(adaptation(ind), c, length(adaptation(ind))/c),1);  % smooth by block (c=8)
    adaptationblock = reshape(adaptationblock,length(adaptationblock)/nsubj,nsubj);
    
    % smooth adaptation by 2 blocks
    c = c*smooth_factor;
    adaptationblockexp = mean(reshape(adaptation(indexp), c, length(adaptation(indexp))/c),1); % smooth every 2 blocks (c=16) in exposure phase
    adaptationblockexp = reshape(adaptationblockexp,length(adaptationblockexp)/nsubj,nsubj);
    
    % plot by subject in group
    figure();
    for i=1:nsubj
        subplot(1,6,i);
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
    
    if group==4           % if it's the planning only condition, get rid of some extra data on the ends
        adaptationblock = adaptationblock(2:end-1,:);
    end
    
    % average across subjects
    adaptblock_mean = mean(adaptationblock,2);
    adaptblock_se = std(adaptationblock,0,2) ./ sqrt(size(adaptationblock,2));
    adaptblockexp_mean = mean(adaptationblockexp,2);
    adaptblockexp_se = std(adaptationblockexp,0,2) ./ sqrt(size(adaptationblockexp,2));
    
    % plot average across all subjects in group
    figure(1000);
    subplot(2,1,2);
    shadeplot(1:pre, adaptblock_mean(1:pre), adaptblock_se(1:pre),'-',colours(group,:),0.3); hold all; % pre-exposure
    shadeplot(linspace(pre+1,exp+pre,length(adaptblockexp_mean)), adaptblockexp_mean, adaptblockexp_se,'-',colours(group,:),0.3); hold all;  % exposure
    h = shadeplot(exp+pre+1:exp+pre+post, adaptblock_mean(end-post+1:end), adaptblock_se(end-post+1:end),'-',colours(group,:),0.3); hold all;  % pre-exposure
    plot([0 pre+exp+post],[0 0],'k');
    ylabel('NMPE (cm)');
    xlabel('Block');
    Pg = P;
    Pg.Vertices = [ pre+.5 -100; exp+pre+.5 -100; exp+pre+.5 100; pre+.5 100];
    if group ==1
        patch(Pg); hold on;
    end
    axis([0 pre+exp+post -20 60]);
    fh = [fh,h];
    
    %% Plot the trajectories for each group, at different experiment phases
    % Note: this takes a long time to execute
    if FLAGS.plotextra
        
        if group<4
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
            
            n1 = find(((field(:,S.FIELD_VISCOUS)~=1) & (fdir==1))  & (round(D.HomeAngle)==i)); %leftward target    0 degrees
            n2 = find(((field(:,S.FIELD_VISCOUS)~=1) & (fdir==-1))  & (round(D.HomeAngle)==i)); %rightward target
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
            shadedTrajectory(trajectory_x.pre1', trajectory_y.pre1', [0.8 0.2 0.2], 0.1); hold on;
            shadedTrajectory(trajectory_x.pre2', trajectory_y.pre2', [0.2 0.2 0.8], 0.1); hold on;
            axis equal;
            
            subplot(1,4,2);   % early-exposure
            shadedTrajectory(trajectory_x.early1', trajectory_y.early1', [0.8 0.2 0.2], 0.1); hold on;
            shadedTrajectory(trajectory_x.early2', trajectory_y.early2', [0.2 0.2 0.8], 0.1); hold on;
            axis equal;
            
            subplot(1,4,3);   % late-exposure
            shadedTrajectory(trajectory_x.late1', trajectory_y.late1', [0.8 0.2 0.2], 0.1); hold on;
            shadedTrajectory(trajectory_x.late2', trajectory_y.late2', [0.2 0.2 0.8], 0.1); hold on;
            axis equal;
            
            subplot(1,4,4);   % post-exposure
            shadedTrajectory(trajectory_x.post1', trajectory_y.post1', [0.8 0.2 0.2], 0.1); hold on;
            shadedTrajectory(trajectory_x.post2', trajectory_y.post2', [0.2 0.2 0.8], 0.1); hold on;
            axis equal;
            clear trajectory_x trajectory_y
        end
        strng = sprintf('Average hand paths, group %d ',group);
        suptitle(strng);
    end
end
%%
toc


%figure(1000);
%legend([fh(1)',fh(2)',fh(3)',fh(4)'],{'Full follow-through','No follow-through','Execution only','Planning only'});


% Run stats on anything you want to run stats on (don't bother doing this
% again...?)




%{
    
    %% Statistical Analysis
  
    
    
    %% Learning stats. analysis within each group
    
    % Original (Howard et al) used first 4 blocks (= 64 exposure trials, 8
    % channel trials) and last 4 blocks of data for epochs in their repeated
    % measures ANOVA (epoch =2, levels=2 for experimental group).
    
    % Compare within each condition the first and last epoch for MPE.
    
    for CondNum=1:ConditionNum
        
        MPE = MPESUBS{CondNum};
        Adaptation = ADAPTATIONBOTHBLOCK{CondNum};
        
        SubTotal = length(SUBJECTLIST.Condition{CondNum});
        
        if (CondNum ==1)
            for SubNum=1:SubTotal
                
                % Extract trial numbers for the first 64 (8 blocks) exposure trials
                % and the first 8 averaged channel trials (8 blocks) within the exposure
                % phase. Dont used the data smoothed over 2 blocks for
                % this, just smoothed over one block so its accurately
                % taking the right sections of information.
                
                Targets    = Paradigm.Targets(:,:,SubNum);
                Field      = Paradigm.Field(:,:,SubNum);
                
                meanMPEfirst = mean(MPE(6:13,SubNum));
                meanMPElast = mean(MPE(end-10:end-3,SubNum));
                
                meanADAPTfirst = mean(Adaptation(11:26,SubNum));
                meanADAPTlast = mean(Adaptation(end-21:end-6,SubNum));
                
                % Group into format for anova_rm.m analysis
                FullMPE(SubNum,:) = [meanMPEfirst, meanMPElast];
                FullADAPT(SubNum,:) = [meanADAPTfirst, meanADAPTlast];
                
            end
            
            MPECON(CondNum,:,:) = FullMPE;
            MPEADAPT(CondNum,:,:) = FullADAPT;
            
            % Full FT ranova_rm.m
            [p.LearningFull.MPE, Stattable.LearningFull.MPE] = anova_rm(FullMPE);  %p-value = 0.0013
            [p.LearningFull.Adapt, Stattable.LearningFull.Adapt] = anova_rm(FullADAPT);  %p-value = 0.0057 therefore learning within FullFT.
            
         end
            
        % No FT
        
          if (CondNum ==2)
            for SubNum=1:SubTotal
                
                % Extract trial numbers for the first 64 (8 blocks) exposure trials
                % and the first 8 averaged channel trials (8 blocks) within the exposure
                % phase. Dont used the data smoothed over 2 blocks for
                % this, just smoothed over one block so its accurately
                % taking the right sections of information.
                
                
                meanMPEfirst = mean(MPE(6:13,SubNum));
                meanMPElast = mean(MPE(end-10:end-3,SubNum));
                
                meanADAPTfirst = mean(Adaptation(11:26,SubNum));
                meanADAPTlast = mean(Adaptation(end-21:end-6,SubNum));
                
                % Group into format for anova_rm.m analysis
                NoFTMPE(SubNum,:) = [meanMPEfirst, meanMPElast];
                NoFTADAPT(SubNum,:) = [meanADAPTfirst, meanADAPTlast];
                
            end
            
            MPECON(CondNum,:,:) = NoFTMPE;
            MPEADAPT(CondNum,:,:) = NoFTADAPT;
            
            % No FT anova_rm.m
            [p.LearningNoFT.MPE, Stattable.LearningNoFT.MPE] = anova_rm(NoFTMPE);  %p-value = 0.3925
            [p.LearningNoFT.Adapt, Stattable.LearningNoFT.Adapt] = anova_rm(NoFTADAPT);  %p-value = 0.7878 therefore no learning in No FT.
            
         end
        
        
        % Appear
        
          if (CondNum ==3)
            for SubNum=1:SubTotal
                
                % Extract trial numbers for the first 64 (8 blocks) exposure trials
                % and the first 8 averaged channel trials (8 blocks) within the exposure
                % phase. Dont used the data smoothed over 2 blocks for
                % this, just smoothed over one block so its accurately
                % taking the right sections of information.
                
                meanMPEfirst = mean(MPE(6:13,SubNum));
                meanMPElast = mean(MPE(end-10:end-3,SubNum));
                
                meanADAPTfirst = mean(Adaptation(11:26,SubNum));
                meanADAPTlast = mean(Adaptation(end-21:end-6,SubNum));
                
                % Group into format for anova_rm.m analysis
                AppearMPE(SubNum,:) = [meanMPEfirst, meanMPElast];
                AppearADAPT(SubNum,:) = [meanADAPTfirst, meanADAPTlast];
                
            end
            
            MPECON(CondNum,:,:) = AppearMPE;
            MPEADAPT(CondNum,:,:) = AppearADAPT;
            
            % Appear anova_rm.m
            [p.LearningAppear.MPE, Stattable.LearningAppear.MPE] = anova_rm(AppearMPE);  %p-value = 0.0386, some significant decrease in MPE, but no force learning
            [p.LearningAppear.Adapt, Stattable.LearningAppear.Adapt] = anova_rm(AppearADAPT);  %p-value = 0.6382 no learning.
           
            
          end
         
          
        % Disappear (Note different block lengths for disappear condition)
          if (CondNum == 4)
            for SubNum=1:SubTotal
                
                % Extract trial numbers for the first 64 (8 blocks) exposure trials
                % and the first 8 averaged channel trials (8 blocks) within the exposure
                % phase. Dont used the data smoothed over 2 blocks for
                % this, just smoothed over one block so its accurately
                % taking the right sections of information.
                
                meanMPEfirst = mean(MPE(6:13,SubNum));    % this is the same format as other groups since the first and last 8 MPE null trials haven't been saved here
                meanMPElast = mean(MPE(end-10:end-3,SubNum));
                
                meanADAPTfirst = mean(Adaptation(13:28,SubNum));
                meanADAPTlast = mean(Adaptation(end-23:end-8,SubNum));
                
                % Group into format for anova_rm.m analysis
                DisappearMPE(SubNum,:) = [meanMPEfirst, meanMPElast];
                DisappearADAPT(SubNum,:) = [meanADAPTfirst, meanADAPTlast];
                
                
            end
            
            MPECON(CondNum,:,:) = DisappearMPE;
            MPEADAPT(CondNum,:,:) = DisappearADAPT;

            % Disappear anova_rm.m
            [p.LearningDisappear.MPE, Stattable.LearningDisappear.MPE] = anova_rm(DisappearMPE);  %p-value = 0.00039
            [p.LearningDisappear.Adapt, Stattable.LearningDisappear.Adapt] = anova_rm(DisappearADAPT);  %p-value = 0.00006
            
          end
          
          diffMPE = MPECON(CondNum,:,2) - MPECON(CondNum,:,1);
          diffADAPT = MPEADAPT(CondNum,:,2) - MPEADAPT(CondNum,:,1);
          
          DiffMPEBar_m(CondNum) = mean(diffMPE);   %late exposure (last 4 blocks) minus early exposure (first four blocks)
          DiffADAPTBar_m(CondNum) = mean(diffADAPT);
          
          DiffMPEBar_se(CondNum) = std(diffMPE) ./ (sqrt(length(diffMPE)));
          DiffADAPTBar_se(CondNum) = std(diffADAPT) ./ (sqrt(length(diffADAPT)));

    end
    
      %--------------
    
    %% Between group learning stats. analysis
    
    % Contrast adaptation measures across groups to check for interaction
    % between epoch and group (look at previous paper for detail about what
    % this means...)
    
  
    % This is the only important one
    
       % Full FT vs. Disappear Condition  > should be the same
    [p.LearningBetweenFullFT_Disappear.MPE, Stattable.LearningBetweenFullFT_Disappear.MPE] = anova_rm({FullMPE, DisappearMPE});  %p-value group = 0.365 ; p-value interaction = 0.0594
    [p.LearningBetweenFullFT_Disappear.Adapt, Stattable.LearningBetweenFullFT_Disappear.Adapt] = anova_rm({FullADAPT, DisappearADAPT});  %p-value group = 0.931 (most important); p-value interaction = 0.516
  
  
%}

