function [xhat, errbar, alpha0] = VBML(Phi,y,netSize)
%% initialize
[m,n] = size(Phi);
tol = 1e-8; %1e-10
iter = 10000;    % 10000
dis = 2;
a0 = 1e-6;
b0 = 1e-6;
c0 = 0;
d0 = 0;
j0 = 1e-6;
h0 = 1e-6;

K = 2;
F = 2*K/(2*K+1);
rho = K/(2*K+1);
f = K/(F/(1-F)-rho/(1-rho));
e = rho/(1-rho)*f;
e0 = e; f0 = f;
e_new = zeros(n,1);
f_new = zeros(n,1);


for i=1:m
    y(i) = y(i)/norm(Phi(i,:));
    Phi(i,:) = Phi(i,:)/norm(Phi(i,:));
end
%%
Pty = Phi'*y;
mu = ones(n,1)*std(y)/1.5;
z = ones(n,1);
Sigma = zeros(n,n);
res = inf*ones(iter + dis,1);
alpha0_temp = [];

sigma2 =  std(y)^2/1e2;
PHI2 = sum(Phi.^2)';            %对每一列求和，再变为一列
ratio = (Pty.^2)./PHI2;
[maxr,index] = max(ratio);
alpha = PHI2(index)/(maxr-sigma2);

for t = 1:iter
    % update lambda
    j = j0;
    h = h0 + alpha;
    lambda = j./h;
   
    % update alpha
    c = c0+1/2; d = lambda + ( mu.^2 + diag(Sigma) ) /2;
    alpha = c./d;
%     alpha = -0.5./lambda + sqrt(0.25./(lambda^2) + (mu.^2 + diag(Sigma)./lambda) );

    % update alpha0
    a = a0+m/2;
    zzt = z*z'+diag(z)*diag(1-z);
    wwt = mu*mu'+ Sigma;
    b = b0 + (y'*y - 2*y'*Phi*(mu.*z) + ones(1,n)*(zzt.*wwt.*(Phi'*Phi))*ones(n,1))/2;
    alpha0 = a/b;
    alpha0_temp = [alpha0_temp;alpha0];
    % update mu and Sigma
    BB = alpha+alpha0*z.*(1-z).*diag(Phi'*Phi); D = diag(1./BB);
    Phin = Phi*diag(z);
    Sigma_new = 1/alpha0*eye(m)+Phin*D*Phin';
    Sigma = D - D*Phin'*(Sigma_new\(Phin))*D;
    mu = alpha0*Sigma*diag(z)*Pty;
    x = mu;
    restt = y - Phi*x;
    
    %% 
    yj = repmat(restt,[1,n])+Phi.*repmat(x,[1,m])';
    
%     for i = 1:n
% %         fprintf( '%d \n',i);
%         neigb = neig(i, netSize);
%         z_f = z(i) + z(neigb);
%         e_new(i) = e0 + z_f;
%         f_new(i) = f0 + 2 - z_f;
%     end
% 
%     ln_pi = psi( e_new ) - psi(e_new + f_new);
%     ln_pi_1 = psi( f_new ) - psi(e_new + f_new);
%     % update z
%     q_0 = exp(ln_pi_1);
%     temp = -alpha0*0.5*((mu.^2+diag(Sigma)).*diag(Phi'*Phi)...
%         -2*mu.*diag(Phi'*yj));
%     temp(temp>=1e2) = 1e2;
%     q_1 = exp (ln_pi) .* exp (temp);
%     z = q_1./(q_0+q_1);

    %% 
%     
%     if options.dis
%         pause(0.1);
%         figure(1000)
%         subplot(3,1,1),plot(z);axis([0,n,0,1])
%         subplot(3,1,2),plot(temp);
%         subplot(3,1,3),plot(res);
%     end

    Sigmay = alpha0^(-1)*eye(m,m)+(Phi*diag(z))*diag(alpha.^(-1))*(Phi)';
    dSigmay = det(Sigmay);
    if dSigmay == 0
        break;
    end
    res(t+dis) = log(dSigmay)+y'*inv(Sigmay)*y;
    
    if   abs(res(t+dis)-mean(res(t:t+dis-1)))/abs(mean(res(t:t+dis-1))) < tol
        break;
    end
end
xhat = mu.*z;
errbar = sqrt(diag(Sigma).*z);
fprintf('Cluster VB using laplace prior, iterations: %d\n',t);

end

function neigb = neig(i,size)
    fir = floor(i/size);
    sec = i - fir*size;
    if sec ~=0
        sec = sec -1;
        neigb = fir + sec * size + 1;
        return;
    end
    if sec == 0
        neigb = fir;
    end
end