% evaluatee rand test (two-sample test)
cd('C:\Dropbox\Projects\Randomization test\exp\t-test')

%choose the test type: "mean": difference in means; "t": t-test
%test_type = "mean";
test_type = "t";

%num_mu = 1;
num_mu = 20;

mu = linspace(0, 3, num_mu);
%mu = 0;
n1 = 15;
n2 = n1;

%num_rep = 1000;
num_rep = 100000;
K=99;
names = ["Deterministic", "Permutation"];
rej = zeros(length(names), num_mu, num_rep);

alpha= 0.05;
%t = tinv(1-alpha/2, n1+ n2-2);
index = ceil((1-alpha)*(K+1));
rng(2);

%% simulation
for i=1:num_mu
    disp(i);
    S1= mu(i);
    
    for j = 1:num_rep
        N1 = randn(n1,1);
        X1 = ones(n1,1)*S1+ N1;
        N2 = randn(n2,1);
        X2 = N2;
        
        %Deterministic
        if (test_type == "t")
            stat = t2(X1,X2,n1,n2);
            rej(1,i,j)=ttest2(X1,X2);
        end
        if (test_type == "mean")
            stat = (mean(X1)-mean(X2))/sqrt(1/n1+1/n2);
            pv = 1- normcdf(stat);
            rej(1,i,j)= (pv<alpha);
        end
        
        
        
        %Randomized
        perm_stat = zeros(K,1);
        for k=1:K
            X = [X1;X2];
            gX = X(randperm(length(X)));
            gX1 = gX(1:length(X1));
            gX2 = gX(length(X1)+1:length(X));
            
            if (test_type == "t")
                perm_stat(k) =  t2(gX1,gX2,n1,n2);
            end
            if (test_type == "mean")
                perm_stat(k) =  (mean(gX1)-mean(gX2))/sqrt(1/n1+1/n2);
            end
            
        end
        x  = sort(perm_stat);
        thresh = x(index);
        if stat> thresh
            rej(2,i,j)=1;
        end
    end
end

%%
pow = mean(rej(:, :, :),3);

%% plot
figure, hold on;
mark = {':', '-'};
rng(2);

plot(mu, pow(1, :), 'lineWidth', 2, 'color',rand(1,3), 'DisplayName', names(1), 'linestyle', mark{1});
plot(mu, pow(2, :), 'lineWidth', 2, 'color',rand(1,3), 'DisplayName', names(2), 'linestyle', mark{2});
legend('location','southeast');
xlabel('$$\mu$$', 'Interpreter', 'LaTex');
ylabel('Power');
ylim([0, 1]);
xlim([0, max(mu)]);
set(gca,'fontsize',18)
grid on;
if (test_type == "t")
    filename = sprintf('t-test_pow_n_%d_K_%d_num_mu_%d_nrep_%d.png', n1, K, num_mu, num_rep);
end
if (test_type == "mean")
    filename = sprintf('mean_pow_n_%d_K_%d_num_mu_%d_nrep_%d.png', n1, K, num_mu, num_rep);
end

saveas(gcf, filename);

%% optionally, save data to file
%save('t-test-exp.mat')

