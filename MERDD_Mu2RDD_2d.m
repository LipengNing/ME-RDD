function p = MERDD_Mu2RDD_2d(Lambda,spec)
% convert Lambda to DDF 

delta_b = spec.delta_b;
delta_t = spec.delta_t;
D0=spec.D0;
r0 = spec.r0;
Nd = spec.Nd;
Nr = spec.Nr;
Tmin = spec.Tmin;
D = linspace(0,D0,Nd);
r = linspace(0,r0,Nr);

switch(spec.method)
   case  'Delta1'
    p = eval_ddf_Delta1(Lambda,delta_b,delta_t,D,r,Tmin);
   case 'Delta1KL'
    p = eval_ddf_Delta1KL(Lambda,delta_b,delta_t,D,r,Tmin);
   case 'Delta2'
    p = eval_ddf_Delta2(Lambda,delta_b,delta_t,D,r,Tmin);
   case 'Delta2KL'
    p = eval_ddf_Delta2(Lambda,delta_b,delta_t,D,r,Tmin);
    otherwise
    disp('Unknown method.')

end

end

function p = eval_ddf_Delta1(lambda,delta_b,delta_t,D,r,Tmin)
    [Kb,Kt] = size(lambda);
    b = [0:Kb-1]'*delta_b;
    t = Tmin+[0:Kt-1]'*delta_t;
    first_term = 0;
    for k = 1:Kb
        for l = 1:Kt
            Mkl = -lambda(k,l)*exp(-b(k)*D')*exp(-t(l)*r);
            first_term = first_term + Mkl;
        end
    end

    p= exp(first_term-1);
end





function p = eval_ddf_Delta1KL(lambda,delta_b,delta_t,D,r,Tmin)
[Kb,Kt] = size(lambda);
b = [0:Kb-1]'*delta_b;
t = Tmin+[0:Kt-1]'*delta_t;
first_term = 0;
for k = 1:Kb
    for l = 1:Kt
        Mkl = -lambda(k,l)*exp(-b(k)*D')*exp(-t(l)*r);
        first_term = first_term + Mkl;
    end
end

prior_term = exp(.5*D(:))*exp(2*Tmin*r);
p= exp(first_term-1).*prior_term;


end


function p = eval_ddf_Delta2(lambda,delta_b,delta_t,D,r,Tmin)
[Kb,Kt] = size(lambda);
b = [0:Kb-1]'*delta_b;
t = [0:Kt-1]'*delta_t;
first_term = 0;
for k = 1:Kb
    for l = 1:Kt
        Mkl = -lambda(k,l)*exp(-b(k)*D')*exp(-t(l)*r);
        first_term = first_term + Mkl;
    end
end

p= exp(first_term);
p = p.*(ones(length(D),1)*exp(Tmin*r));
end


function p = eval_ddf_Delta2KL(lambda,delta_b,delta_t,D,r,Tmin)
[Kb,Kt] = size(lambda);
b = [0:Kb-1]'*delta_b;
t = [0:Kt-1]'*delta_t;
first_term = 0;
for k = 1:Kb
    for l = 1:Kt
        Mkl = -lambda(k,l)*exp(-b(k)*D')*exp(-t(l)*r);
        first_term = first_term + Mkl;
    end
end

%sec_term = ones(length(D),1)*exp(5*Tmin*r);
sec_term = exp(.5*D(:))*exp(2*Tmin*r);
p= exp(first_term).*sec_term;
p = p.*(ones(length(D),1)*exp(Tmin*r));
end





