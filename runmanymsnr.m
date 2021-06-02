function [pksteps,stats] = runmanymsnr(sigsin,samprate,trpks,segsh,mtchoset)
%	The input record   
%
% Inputs: sigsin - m x n array of signals in double format, each row corresponding to one channel;
%                  the polarity of the signal is preferably selected so that the heart beats
%                  are characterized by dominant upward/positive peaks 
%                  the signals will be low pass filtered to 45 kHz; they should
%                  not be prefiltered in a manner that results in oversmoothing
%         trpks - time of true peaks (optional) in ms
%         samprate - sampling rate in ms
%         segsh - 5 second segment to run; if left empty, all segments will be run
%         mtchoset - if trpks is input, the time offset between the 
%                    peak time in trpks and the time of the detected peak
%                    by this function; determined by experimentation
%
% Outputs: pksteps - peak times in ms 
%          stats - various statistics
%
%
% Other m-files required: domultisnr
% Subfunctions: 
% Packages required: signal, nnet
% NOTE:  The getcandpeaks function, called by the chain runmanymsnr/domultimsnr, runs a slightly modified version
% version of the Octave findpeaks function (Version 1.4.1 of the signal package).  In particular,
% a bug in the peak width assignment statement was fixed, and the parabola fitted to the vertex was
% normalized according to the height to allow a single peak width paramter to govern signals with
% a different amplitudes: pp(1) = (ind-xm).^2 \ (data(ind)/H-1)

%
% 
% Author: Bruce Hopenfeld
% May 2021
% The code is covered by US patent 9402557, pending patent applications
% and copyright; license to use the code is granted solely for research purposes
%------------- BEGIN CODE --------------

disp("The code is covered by US patent 9402557, pending patent applications")
disp("and copyright; license to use the code is granted solely for research ")
disp("according to the terms of the license available on the GitHub license")
disp("page associated with this software.")
 
if ~exist('trpks','var')
  trpks=[];
end

scorenodrlast=0;
tmslast=[];
nseg=5*256;
rrsteps=[];
pksteps=[];
rrlast=[];
stats=[];

numchan=length(sigsin(:,1));

[b,a]=butter(5,45/(samprate/2));

for kk=1:numchan
  sigsf(kk,:)=resample(filtfilt(b,a,sigsin(kk,:)),256,samprate); 
end

if ~isempty(segsh)
  segproc=segsh;
else
  segproc=[1:floor(length(sigsf(1,:))/1280)];
end

pklast=[];
for iictr = 1:length(segproc) 
    ii=segproc(iictr);
    svr=round(nseg*(ii-1)+[1:nseg]);
    sigsr=sigsf(:,svr); %get segments
    tmoset=5000*(ii-1);
    trpksa=trpks(find(trpks > tmoset & trpks <= tmoset+5000));
    trpksu=trpksa-tmoset;
    
    for kk=1:numchan  %compute second difference signals to assess correlation
      ddsigs(kk,:)=sigsr(kk,3:end)-2*sigsr(kk,2:end-1)+sigsr(kk,1:end-2); %compute second difference to get correlation
    end
    
    cormat=tril(corr(ddsigs'),-1)';  
                                     
    [rc,cc]=find(abs(cormat) > 0.7);%add together channels that are highly correlated
                                     %implemented to handle no more than a single pair
                                     %of highly correlated channels
    if ~isempty(rc) & 1==0
      amp2=std(sigsr(rc(1),:));
      amp3=std(sigsr(cc(1),:));
      mn2=mean(sigsr(rc(1),:));
      mn3=mean(sigsr(cc(1),:));
      sigsr(rc(1),:)=(amp2*(sigsr(rc(1),:)-mn2)+amp3*(sigsr(cc(1),:))-mn3)/(amp2+amp3);
      sigsr(cc(1),:)=[];
    end

    tmsallref=[];
    statsa=domultisnr(sigsr,trpksu,256,mtchoset);
    datuse=statsa.dat;
    tmsall=datuse(:,3:end);
    scoreso=datuse(:,2);
    rrsa=datuse(:,1);
    stats{ii}=statsa;

    if isempty(datuse)
      rrsteps(ii)=NaN;
      continue
    end
    mind=0;
    if ~isempty(tmslast)
     [scores,rrso]=meshsegs(tmsall,tmslast,rrsa,scoreslast,scoreso); %enhance scores for sequences that mesh
     [mxscr,mind]=max(scores);
     if ~isempty(trpks)
        %nmtch=getmatches(tmslast(mind,:),trpkslast);
        trpkslast=trpksu;  
     end
    end
    if mind > 0 & ii > 1
      pkstepsl(ii-1,1:length(tmslast(mind,:)))=tmslast(mind,:);
      rrsteps(ii-1)=rrso(mind);
    elseif ii > 1
      rrsteps(ii-1)=NaN;
      pksteps(ii-1)=NaN;
    end
tmslast=tmsall; %save sequence times and scores to mesh with next segment
scoreslast=scoreso;
    
end
  
endfunction
   
