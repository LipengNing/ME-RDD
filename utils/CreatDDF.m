function  p = CreatDDF(spec)
w = spec.w; % weight of different components
center = spec.center;
sig = spec.sig;

N = length(w);
D = spec.D;
R = spec.R;
delta_D = D(2)- D(1);
delta_R = R(2) - R(1);
Nd = length(D);
Nr = length(R);
p = zeros(Nd,Nr,N);


for i = 1:N
    cD = center(i,1);
    cR = center(i,2);
    XD = (D(:)-cD)*ones(1,Nr);
    XR = ones(Nd,1)*(R(:)'-cR);
    X = [XD(:) XR(:)];
    if(ismatrix(sig))
        X = sum(X.^2,2);
        X = X./(2*sig(i));
    else
        M = inv(sig(:,:,i));
        X = sum((X*M).*X,2)/2;
    end
    ppi = exp(-X);
    ppi = ppi/sum(ppi(:))/(delta_D*delta_R);
    ppi = reshape(ppi,Nd,Nr)*w(i);
    p(:,:,i) = ppi;
end
p = sum(p,3);    
end

