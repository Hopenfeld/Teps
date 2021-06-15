
function [scoresn,nsa] = meshsegs(tmsnxt,tmscur,rrsnxt,rrscur,scores,scoresnxt)
  %if the last part of a sequence from the previous segment (tmslast) temporally 
  %meshes with the start of a sequence from the current segment (tmscur)
  %add a fraction of the previous sequence's score to the present sequence's score
  
  scoresn=scores;
  rrso=zeros(1,length(tmscur(:,1)));
  nsa=zeros(length(tmsnxt(:,1)),1);

    for jj=1:length(tmscur(:,1))
      ns=0;
      tmb=tmscur(jj,find(tmscur(jj,:)));
      tmb=tmb(find(tmb>2500)); 
      if length(tmb) < 2
        continue
      end
      
        rrn=rrscur(jj);
        
      for ii=1:length(tmsnxt(:,1)) 
        tma=tmsnxt(ii,find(tmsnxt(ii,:)))+5000;
        tma=tma(find(tma<7500));
        if length(tma) < 2
          continue
        end
        seql=[tmb tma];

        rrd=rrsnxt(ii);
        slst=scoresnxt(ii);
        if abs(rrn - rrd) > 200
          continue
        end
        rrm=mean([rrn rrd]);
        [regsc,rrl,skips,seql] = getmaxscore(seql,rrm,1);
        regsc=regsc/(length(seql)-2);
        ovap=~isempty(find(seql<=5000)) & ~isempty(find(seql >5000));  
        
        if  (regsc > -700 & ovap) 
          
          sfa=min(1,(rrm+double(regsc))/350);
          sfb=max(0,(40-abs(rrn-rrd))/40);
          ns=max(ns,int16(max(0,slst)*sfa));
          ns=max(ns,int16(max(0,slst)*sfb));
          
          nsa(ii)=max(nsa(ii),int16(max(0,scores(jj))*sfa));
          nsa(ii)=max(nsa(ii),int16(max(0,scores(jj))*sfb));

        end

      end
      scoresn(jj)=scoresn(jj)+ns;

  end
  
  
endfunction



