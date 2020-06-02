function opval=loglike(y,yhat)

if nargin==0
    disp('no data passed to the function loglike, using default values')
    y=[1 1 0 1 1 0 0 1 1 1 0 0 1 1 0];
    disp('no data passed to the function loglike, using default values')
    yhat=[0.899 0.540 0.317 0.561 0.698 0.457 0.234 0.899 0.561 0.764 0.457 0.561 0.698 0.457 0.899]; %corrected from Tab and Fid
    %yhat=0.6.*ones(1,15);
elseif nargin==2
    %use parameters passed
else
    disp('wrong number of parameters passed to nargin')
end

llvec=(y.*log(yhat) +(1-y).*log(1-yhat));
opval=sum(llvec);

