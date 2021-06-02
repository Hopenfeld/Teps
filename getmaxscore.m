function [lhood,rrl,skips,tmsin] = getmaxscore(tmsin,rrl,msh)
load multichanparams stdarr rrmult

lhood=[];
if exist('msh','var')
  mflg=1;
else
  mflg=0;
end


if rrl==0
  return 
end
tmsin=tmsin(find(~isnan(tmsin) & tmsin > 0));
rrs=tmsin(2:end)-tmsin(1:end-1);
rrs=rrs(:)';
if mflg
    [rrl,stdv,ns,tmsin] = fixrr(rrl,rrs,tmsin);
    rrs=tmsin(2:end)-tmsin(1:end-1);
    rrs=rrs(:)';
  else
  [rrl,stdv,ns] = fixrr(rrl,rrs);
end
rruse=rrl;
tml=rrs/rruse;
dus=max(1,round(tml));
skips=sum(dus)-length(dus);

if length(find(dus==0)>0) || rruse == 0 || rruse < 300
   return
end
    
rrsd=rrs./dus;
drrsd=abs(rrsd(2:end)-rrsd(1:end-1));
drr=double(abs(rrs(2:end).*dus(1:end-1)-rrs(1:end-1).*dus(2:end)));
dusu=dus-1;
ddus=dusu(1:end-1)+dusu(2:end);
dususe=max(dus(1:end-1),dus(2:end));
stuse=stdarr(dususe)*rruse*rrmult;
drr(find(drr==0))=1;
s0p=sum(-(drr./stuse).^2/2);
s0p1=sum(0.5*log(2*pi)+log(stuse));
lhood=int16(100*(s0p-s0p1));
%-(drr./stuse).^2/2-(0.5*log(2*pi)+log(stuse))
if rrl < 1000
  tol=60;
else
  tol=150;
end

skips=skips+floor(tmsin(1)/((rrl+tol)));
skips=skips+floor((5000-tmsin(end))/((rrl+tol)));

  

  
   
   

    
 