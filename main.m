function [pksteps,stats,statsr,statst] = main(sigsin,samprate,trpks,segsh,mtchoset)
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
statsr=[];
statst=[];
scoresaugcur=[];

numchan=length(sigsin(:,1));

[b,a]=butter(5,45/(samprate/2));

for kk=1:numchan
  sigsb=filtfilt(b,a,sigsin(kk,:));
  sigsf(kk,:)=resample(sigsb,256,samprate); 
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

    for kk=1:numchan  %compute second difference signals to assess correlation
      ddsigs(kk,:)=sigsr(kk,3:end)-2*sigsr(kk,2:end-1)+sigsr(kk,1:end-2); %compute second difference to get correlation
    end
    
    cormat=tril(abs(corr(ddsigs')),-1)';                                
    [rc,cc]=find(cormat > 0.6);%add together channels that are highly correlated
                                     %implemented to handle no more than a single pair
                                     %of highly correlated channels
    if ~isempty(rc) 
      amp2=std(sigsr(rc(1),:));
      amp3=std(sigsr(cc(1),:));
      mn2=mean(sigsr(rc(1),:));
      mn3=mean(sigsr(cc(1),:));
      sigsr(rc(1),:)=(amp2*(sigsr(rc(1),:)-mn2)+amp3*(sigsr(cc(1),:))-mn3)/(amp2+amp3);
      sigsr(cc(1),:)=[];
      cormat(rc(1),cc(1))=0;
    end
    coradj=1-min(1,3*max(0,max(cormat(:))-0.3)); 
    
    if ~isempty(trpks)
      [trpksu,statst] = procref(trpks,tmoset,statst,ii);
    else
      trpksu=[];
    end

    statsa=domultisnr(sigsr,trpksu,256,mtchoset,coradj);
    stats{ii}=statsa;
    datuse=statsa.dat;
    %round(datuse)
    tmsall=datuse(:,5:end);
    scoreso=datuse(:,2);
    rrsa=datuse(:,1);
    
    if isempty(datuse)
      rrsteps(ii)=NaN;
      continue
    end
    mind=0;
    if ~isempty(tmslast)
      mold=max(scoreslast);
      if isempty(scoresaugcur)
          [scores,scoresaugcur]=meshsegs(tmsall,tmslast,rrsa,rrslast,scoreslast,scoreso,[]); %enhance scores for sequences that mesh
          %[tmslast scores scoreslast]

      else
          [scores,scoresaugcur]=meshsegs(tmsall,tmslast,rrsa,rrslast,scoreslast,scoreso,scoresaugcur); %enhance scores for sequences that mesh
          %[tmslast scores scoreslast]
      end 
     [mxscr,mind]=max(scores);

     if ~isempty(trpks)
        rrt=mean(trpkslast(2:end)-trpkslast(1:end-1));
        [nmtchs,dp,ngd]=getmatches(tmslast,trpkslast,mtchoset);
        if ~isempty(nmtchs)
          nmtch=nmtchs(mind);
          statsv=[ii-1 nmtch ngd(mind) dp(mind) max(nmtchs)  length(trpkslast) mxscr mold rrslast(mind)];
          %tmslast(mind,:)
        else
          statsv=[ii-1 0 0 0 0  length(trpkslast) mxscr mold rrslast(mind)] ;
        end
        
        statsr(ii-1,:)=statsv;
        %[rrslast(:) scoreslast(:) scores(:) nmtchs(:)]
        
    else
        %[ii-1 mxscr mold] 
        %tmslast(mind,:)
            
    end
    
    else
         nmtchs=getmatches(tmsall,trpksu,mtchoset);
         %[tmsall nmtchs(:) scoreso rrsa]

    end 
    trpkslast=trpksu;

    if mind > 0 & ii > 1
      pkstepsl(ii-1,1:length(tmslast(mind,:)))=tmslast(mind,:);
      rrsteps(ii-1)=rrslast(mind);
 
    elseif ii > 1
      rrsteps(ii-1)=NaN;
      pksteps(ii-1)=NaN;
    end
tmslast=tmsall; %save sequence times and scores to mesh with next segment
rrslast=rrsa;
scoreslast=scoreso;
    
end
  
endfunction

function [trpksu,statst] = procref(trpks,tmoset,statst,indin)

trpksa=trpks(find(trpks > tmoset & trpks <= tmoset+5000));
trpksu=trpksa-tmoset;   
tmsallref=[];
rrstr=trpksu(2:end)-trpksu(1:end-1);
rrstr=mean(rrstr);
[trpksu,rrstr]=permsearch(double(trpksu),500,0.2,3);
regsc=[];
skips=[];
if ~isempty(rrstr)
   for ka=1:length(trpksu(:,1))
      [regsc(ka),rrstr(ka),skips(ka)]=getmaxscore(trpksu(ka,:),rrstr(ka));
   end
   cfs=[-.6 0.002004];
   nseq=sum(trpksu > 0,2);
   lhs=double(regsc(:))./(nseq(:)-2); 
   smat=[skips(:) lhs]';
   scr=cfs*smat;
   [dum,mind]=max(scr);
   trpksu=trpksu(mind,:);
   trpksu=trpksu(find(~isnan(trpksu) & trpksu >0));
   statst(indin,1:length(trpksu)+1)=[scr(mind) trpksu];
else
   statst(indin,1:length(trpksu)+1)=[NaN trpksu];
end
   
endfunction