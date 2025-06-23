%Good nodes set method
function [GD] = Goodnode(M,N)
tmp1 = [1:M]'*ones(1,N);
Ind = [1:N];
prime1 = primes(100*N);
[p,q]=find(prime1>=(2*N+3));
tmp2 = (2*pi.*Ind)/prime1(1,q(1));
tmp2 = 2*cos(tmp2);
tmp2 = ones(M,1)*tmp2;
GD = tmp1.*tmp2;
GD = mod(GD,1);
%% For debuging
% plot(GD(:,1),GD(:,2),'*');
end

