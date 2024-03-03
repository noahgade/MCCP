function band=bandw1(X_s)


[ms,ns]=size(X_s);

Ds=zeros(ns,ns); % D is the pairwise distance matrix
for i=1:ns
    temp=bsxfun(@minus,X_s(:,(i+1):ns),X_s(:,i));
    Ds(i,(i+1):ns)=sum(temp.^2,1);
end
% use rule of thumb to determine the bandwidth
band=median(Ds(Ds~=0));

end