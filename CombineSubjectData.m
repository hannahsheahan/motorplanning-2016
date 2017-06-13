%% FUNCTION: Combines subject data for each group into a single datafile
% INPUTS:  N/A
% OUTPUTS: N/A
% NOTES:    N/A
% ISSUES:   N/A
% REFS:     Wolpert analysis method
% AUTHOR:   Hannah Sheahan, sheahan.hannah@gmail.com
%--------------------------------------------------------------------------
%% settings
clear all; close all; clc;

%% specify subject data
ngroups   = 1;
nsubjects = 6; % per group
datName   = {
   
    % Full follow-through
    'JL-control4d-fullft-24Mar2016-1405'; 
    'MC-control4d-fullft-24Mar2016';
    'DG-control-fullft4d-25mar2016-1003';
    'SM-control-fullft4d-counter-25Mar2016-1216';
    'AM-control4d-fullft-counter-25Mar2016-1429';
    'JY-control4d-fullft-counter-25Mar2016-1645';
    
    %{
    % No follow-through
    'JDL-control-noft4d-16Mar2016-1143';
    'ZW-control-noft4d-16Mar2016-1412';
    'DB-controlnoft4d-17Mar2016-0921';
    'SM-control-noft-counter4d-17Mar2016-1412'; 
    'SWW-control-noft4d-counter-18Mar2016-0914';
    'HE-control-noft-counter4d-18Mar2016-1427';
    
    % Execution only
    'BS-appear4d-21Mar2016-1147';
    'JB-appear4d-21Mar2016-1630';
    'HK-appear4d-22Mar2016-1525';
    'YP-appear4d-counter-23Mar2016-1136';
    'TVP-appear4d-counter-23Mar2016-1411';
    'JS-appear4d-counter-23Mar2016-1710';
    
    % Planning only
    'CP-disappear4dnew-01Apr2016-0956';
    'AA-disappear4dnew-02April2016-1013';
    'JB-disappear4dnew-02April2016-1421';
    'JM-disappear4d-counter-03April2016-1041';
    'LK-disappear4dnew-counter-03April2016-1451';
    'TB-disappear4dnew-counter-04April2016-1005';
    %}
    };

%% specify parameters you want to save
grouplist = 1:ngroups;  % number subjects by group
group  = repmat(grouplist,nsubjects,1);
group  = group(:);
nfiles = size(datName,1);
fdir   = [ones(length(group)/2,1); -1.*ones(length(group)/2,1)];            % mark counterbalance subjects
D = {};

%loop over files to load data and remove elements you don't want
for k=1:nfiles
    sprintf('Subject file: %d...',k);
    ww=['../01 Import Data/data/' datName{k}]; %filenames with .dat/.mat stripped off
    tmp=DATAFILE_Load(ww);
    tmp=DATAFILE_Select(tmp,1:2:tmp.Trials); %only store outward movements
    
    % Remove the fields we don't want
    % Regular fields:
    fields = {'Files'...
        'MovementDurationTooSlowToVia'...
        'MovementDurationTooSlow'...
        'MovementDurationTooFast'...
        'ViaSpeedThreshold'...
        'ViaTimeOutTime'...
        'ViaToleranceTime'...
        'PMoveEndPosition'...
        'PMoveStartPosition'...
        };
    for j=1:length(fields)
        if isfield(tmp,fields(j))
            tmp = rmfield( tmp,fields(j));
        end
    end
    % Framedata fields:
    fields = {'GraphicsVerticalRetraceNextOffsetTime',...
        'GraphicsSwapBuffersToVerticalRetraceTime',...
        'GraphicsSwapBuffersCount',...
        'PMoveStatePosition',...
        'PMoveStateRampValue',...
        'CursorPosition',...
        'HandleTorques',...
        'HandleForces',...
        'ForcesFunctionPeriod',...
        'ForcesFunctionLatency',...
        'StateGraphics',...
        'GraphicsVerticalRetraceNextOnsetTime' ...
        'PMoveStateTime' ...
        'PMoveState' ...
        };
    for j=1:length(fields)
        tmp.FrameData = rmfield( tmp.FrameData,fields{j});
    end
    
    % add the group and subject info to the datafile
    tmp.group=group(k)*ones(size(tmp.TrialNumber))';
    tmp.subj=k*ones(size(tmp.TrialNumber))';
    D = DATAFILE_Append(D,tmp);
    clear tmp;
end

D.datName=datName;
D.fdir=fdir;

% remove the z-dimension from the kinetic/kinematic data (keep files small)
D.FrameData.RobotForces(:,:,3)=[];
D.FrameData.RobotPosition(:,:,3)=[];
D.FrameData.RobotVelocity(:,:,3)=[];

%%
% save it all in a single datafile (much faster to load later!)
save -v7.3 Fullfollowthrough_combined D

