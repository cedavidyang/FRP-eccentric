clear
load ../resistsmp.mat
Nusmp = resistmtr(1,:)'./1e3;
Musmp = resistmtr(2,:)'./1e6;

ftsize = 9;
data = {Nusmp, Musmp};
label = {'N_{u} sample', 'M_{u} sample'};

for i=1:2
    Rsmp = data{i};
    %% prepare data
    [ecdfRsmp, Rsmpunique] = ecdf(Rsmp);
    %lognormal
    mu = mean(Rsmp); sigma = std(Rsmp); cov = sigma/mu;
    logmu = log(mu/sqrt(1+cov^2)); logstd = sqrt(log(1+cov^2));
    startlogn = [logmu, logstd];
    [lognEsts,lognCIs] = mle(Rsmp, 'distribution','lognormal', 'start',startlogn);
    [H, p, ksstat] = kstest(Rsmp,'CDF',[Rsmp,logncdf(Rsmp,lognEsts(1), lognEsts(2))]);
    fcdfRsmp = logncdf(Rsmpunique,lognEsts(1), lognEsts(2));
    iae = sum(abs(ecdfRsmp - fcdfRsmp))/sum(ecdfRsmp);
    disp('lognormal')
    fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);
    
    %% postprocessing

    fig1 = figure;
    stairs(Rsmpunique, ecdfRsmp, 'color', 'k', 'linestyle', '--'); hold on;
    h = plot(Rsmpunique, logncdf(Rsmpunique, lognEsts(1), lognEsts(2)));
    set(h, 'lineWidth', 1);
    xlabel(label{i}); ylabel('CDF')
    postfigs(fig1, 'asce', false, ftsize);
    
    fig2 = figure;
    [~,BinEdge] = histcounts(Rsmp, 'Normalization', 'pdf');
    hLine = histogram(Rsmp, 'Normalization', 'pdf');
    set(hLine,'FaceColor','cyan','EdgeColor','k',...
        'LineStyle','-', 'LineWidth',0.5);
    hold on;
    Rsmpaug = sort([Rsmpunique; BinEdge']);
    plot(Rsmpaug, lognpdf(Rsmpaug, lognEsts(1), lognEsts(2)))
    xlabel(label{i}); ylabel('PDF')
    postfigs(fig2, 'asce', false, ftsize);
end