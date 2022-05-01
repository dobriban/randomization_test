% evaluatee rand test
cd('C:\Dropbox\Projects\Randomization test\exp\vector_rotation')
%p = 10;
p = 100;
num_mu = 20;
mu = linspace(0, 4, num_mu)*sqrt(log(p));

%num_rep = 100;
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
    index = ceil((1-alpha)*(K+1));
    for i=1:num_mu
        disp(i);
        S = zeros(p,1);
        S(1)= mu(i);
        
        for j = 1:num_rep
            N = randn(p,1);
            X = S + N;
            T = max(abs(X));
            
            %Deterministic
            rej(1,i,j,k_ind)=(T>t);
            
            %Randomized
            gT = zeros(K,1);
            for k=1:K
                M = randn(p);
                [O,~,~] = svd(M);
                gX = O*X;
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

plot(mu, pow(1, :, 1), 'lineWidth', 3, 'color',rand(1,3), 'DisplayName', names(1), 'linestyle', mark{1});
plot(mu, pow(2, :, 1), 'lineWidth', 3, 'color',rand(1,3),'DisplayName', names(2), 'linestyle', mark{2});
plot(mu, pow(2, :, 2), 'lineWidth', 3, 'color',rand(1,3),'DisplayName', names(3), 'linestyle', mark{3});
legend('location','southeast');
xlabel('$$\mu$$', 'Interpreter', 'LaTex');
ylabel('Power');
ylim([0, 1]);
xlim([0, max(mu)]);
set(gca,'fontsize',18)
% title(sprintf('$$n=%d, k=%d, \\gamma=%.2f, \\xi=%0.2f$$', n, k, gamma, xi), 'Interpreter', 'LaTex')
grid on;
filename = sprintf('pow_p_%d_K_%d_num_mu_%d_nrep_%d.png', p, K, num_mu, num_rep);
saveas(gcf, filename);


