% -------------------------------------------------------
%
%    hactseg.m
%
%    Ver. 1.0.0
%
%    Created:       Jorge Sanchez Arciniegas      (02.07.2020)
%    Last modified: Jorge Sanchez Arciniegas      (18.02.2021)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology (KIT)
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2021 - All rights reserved.
%
% ------------------------------------------------------
function [accp,timelngth] = Hactseg(bi,smplfrq,a)

bifn = normalize(bi,'center');
th = hilbert(bifn);
phasen = atan(imag(th)./real(th));

method = "env";

if method == "rad"
    middlex = (min(real(th)) + max(real(th)))/2;
    middley = (min(imag(th)) + max(imag(th)))/2;
    distnt = sqrt((middley-imag(th)).^2+(middlex-real(th)).^2);
    nstd = std(distnt)/mean(distnt);
    k=0.5;
    rm = mean(distnt)+k*std(distnt);
elseif method == "env"
    distnt = abs(th);
    if mean(distnt)<median(distnt)
        rm = median(distnt)+1*std(distnt);
    else
        rm = mean(distnt)+1*std(distnt);
    end
else
    disp('Choose a valid method!')
end

phase=0;
N = 1;
acc=zeros(length(distnt)-(N-1),1);
for i=1:length(distnt)-(N-1)
    if distnt(i)>rm
        if phase<=pi*2
            phase = phase + abs(phasen(i+(N-1)));
            acc(i)=1;
        elseif phase>pi*2
            acc(i)=0;
            phase=0;
        end
    end
end

accp = movmax(acc,[6 a]); %Smoothing detected segments Valor origianl 6 -- 20
if length(accp)>5
    accp(end-5:end) = 0;
    accp(1:5) = 0;
end

cnt=0;
k=1;
sgmtlngth=zeros(1,1000);
for i=2:length(accp)
    if accp(i)-accp(i-1)~=accp(i)
        cnt=cnt+1;
    else
        sgmtlngth(k)=cnt;
        cnt=0;
    end
    k=k+1;
end
sgmtlngth(sgmtlngth==0) = [];
timelngth=sgmtlngth.*1000/smplfrq;

%step=1;

indx = find(diff(accp));
try
    indx = reshape(indx.',2,[]);
catch 
    indx = reshape([ 1; indx].',2,[]);
end
indx(3,:) = diff(indx).*1000/smplfrq;
idxtmlgth=find(indx(3,:)<10);

if min(timelngth)<10 % Mininmum physiological EGM duration
    timelngth(idxtmlgth) = [];
    indx(:,idxtmlgth) = [];
    accp = zeros(1,size(bi,1));
    
    for i=1:size(indx,2)
        accp(indx(1,i):indx(2,i)) = 1;
    end
end

end