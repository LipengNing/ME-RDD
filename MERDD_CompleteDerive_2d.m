function [dd,M] = MERDD_CompleteDerive_2d(Mu)
% compute the complete derivative of dwi signal
[Kb,Kt] = size(Mu);
M = [];

for k = 1:Kb
    for l = 1:Kt
        M_kl = [];
        for k_order = 0:k-1
            for l_order = 0:l-1
                if(max(k_order,l_order)>0)
                    ck = DiffCoef(k_order);
                    cl = DiffCoef(l_order);
                    c = ck(:)*cl(:)';
                    X = zeros(Kb,Kt);
                    X(k-k_order:k,l-l_order:l) = c;
                    M_kl = [M_kl; X(:)'];
                end
            end
        end
        M = [M;M_kl];
    end
end

dd = M*Mu(:);
end



function c = DiffCoef(order)
c = 1;
for i = 1:order
    c = conv(c,[1 -1]);
end
end