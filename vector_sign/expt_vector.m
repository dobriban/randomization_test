% evaluatee rand test (two-sample test)
cd('C:\Dropbox\Projects\Randomization test\exp\vector_sign')
%% normal and t distribution
%% the variable tdof sets the degrees of freedom
p = 100;
num_mu = 20;
mu = linspace(0, 1, num_mu);
n1 = 200;
n2 = 200;

num_rep = 1000;
K_arr=[19,99];
names = ["Deterministic", "Randomization K=19", "Randomization K=99"];
rej = zeros(length(names), num_mu, num_rep, 2);

alpha= 0.05;
q = ((1-alpha)^(1/p)+1)/2;
t = norminv(q);
rng(2);

%% simulation
for k_ind = 1:2
    K = K_arr(k_ind);
    index = min(ceil((1-alpha)*(K+1)), K);
    for i=1:num_mu
        disp(i);
        S1 = zeros(p,1);
        S1(1)= mu(i);
        
        for j = 1:num_rep
            N1 = randn(n1,p);
            X1 = ones(n1,1)*(S1')+ N1;
            N2 = randn(n2,p);
            X2 =  N2;
            Del = (mean(X2,1)-mean(X1,1))/sqrt(1/n1+1/n2);
            T = max(abs(Del));
            
            %Deterministic
            rej(1,i,j,k_ind)=(T>t);
            
            %Randomized
            gT = zeros(K,1);
            for k=1:K
                M1 = 2*binornd(1,ones(n1,1)/2)-1;
                M2 = 2*binornd(1,ones(n2,1)/2)-1;
                gX1 = diag(M1)*X1;
                gX2 = diag(M2)*X2;
                gX = (mean(gX2,1)-mean(gX1,1))/sqrt(1/n1+1/n2);
                gT(k) = max(abs(gX));
            end
            x  = sort(gT);
            thresh = x(index);
            if T> thresh
                rej(2,i,j,k_ind)=1;
            end
        end
    end
end

%%
pow = mean(rej(:, :, :, :),3);

%% plot
figure, hold on;
mark = {':', '-.', '-'};
rng(2);
set(gca,'fontsize',18)

plot(mu, pow(1, :, 1), 'lineWidth', 3, 'color',rand(1,3), 'DisplayName', names(1), 'linestyle', mark{1});
plot(mu, pow(2, :, 1), 'lineWidth', 3, 'color',rand(1,3),'DisplayName', names(2), 'linestyle', mark{2});
plot(mu, pow(2, :, 2), 'lineWidth', 3, 'color',rand(1,3),'DisplayName', names(3), 'linestyle', mark{3});
legend('location','southeast');
legend('location','southeast');
xlabel('$$\mu$$', 'Interpreter', 'LaTex');
ylabel('Power');
ylim([0, 1]);
xlim([0, max(mu)]);
grid on;
filename = sprintf('vec_pow_p_%d_K_%d_num_mu_%d_nrep_%d.png', p, K, num_mu, num_rep);
saveas(gcf, filename);

%%
%% t
p = 100;
num_mu = 20;
mu = linspace(0, 1, num_mu);
n1 = 200;
n2 = 200;

num_rep = 1000;
K_arr=[19,99];
names = ["Randomization K=19", "Randomization K=99"];
rej = zeros(length(names), num_mu, num_rep, 2);
tdof = 5;
%tdof = 3;

alpha= 0.05;
rng(2);
%% simulation
for k_ind = 1:2
    K = K_arr(k_ind);
    index = max(min(ceil((1-alpha)*(K+1)), K),1);
    for i=1:num_mu
        disp(i);
        S1 = zeros(p,1);
        S1(1)= mu(i);
        
        for j = 1:num_rep
            N1 = trnd(tdof,n1,p); 
            X1 = ones(n1,1)*(S1')+ N1;
            N2 = trnd(tdof,n2,p);
            X2 =  N2;
            Del = (mean(X2,1)-mean(X1,1))/sqrt(1/n1+1/n2);
            T = max(abs(Del));
            
            %Randomized
            gT = zeros(K,1);
            for k=1:K
                M1 = 2*binornd(1,ones(n1,1)/2)-1;
                M2 = 2*binornd(1,ones(n2,1)/2)-1;
                gX1 = diag(M1)*X1;
                gX2 = diag(M2)*X2;
                gX = (mean(gX2,1)-mean(gX1,1))/sqrt(1/n1+1/n2);
                gT(k) = max(abs(gX));
            end
            x  = sort(gT);
            thresh = x(index);
            if T> thresh
                rej(2,i,j,k_ind)=1;
            end
        end
    end
end

%%
pow = mean(rej(:, :, :, :),3);

%% plot
figure, hold on;
mark = {':', '-.', '-'};
rng(2);

%plot(mu, pow(1, :, 1), 'lineWidth', 3, 'color',rand(1,3), 'DisplayName', names(1), 'linestyle', mark{1});
plot(mu, pow(2, :, 1), 'lineWidth', 3, 'color',rand(1,3),'DisplayName', names(1), 'linestyle', mark{2});
plot(mu, pow(2, :, 2), 'lineWidth', 3, 'color',rand(1,3),'DisplayName', names(2), 'linestyle', mark{3});
legend('location','southeast');
legend('location','southeast');
xlabel('$$\mu$$', 'Interpreter', 'LaTex');
ylabel('Power');
ylim([0, 1]);
xlim([0, max(mu)]);
set(gca,'fontsize',18)
title('$$t$$-distributed noise', 'Interpreter', 'LaTex');
grid on;
filename = sprintf('vec_t_pow_p_%d_K_%d_num_mu_%d_nrep_%d_t_dof_%d.png', p, K, num_mu, num_rep, tdof);
saveas(gcf, filename);

