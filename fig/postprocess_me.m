clear;
databasetype = input('database?', 's');
switch databasetype
    case {'cfrp', 'c'}
        load ../data/modelerror-cfrp
    case {'gfrp', 'g'}
        load ../data/modelerror-gfrp
end
ftsize = 8;

%% linear plot (N)
fig1 = figure;
axes1 = axes('Parent',fig1, 'box', 'on');
hold(axes1,'all');
% N capacity (section)
NpreCell = {NsectionArray, NmodelArray};
datamarker = {'o', '^'};
regstyle = {'--', '-.'};
datacolor = {[0, 0.45, 0.74], [0,0.5,0]};
axis([200 1200 200 1200])
for i=1:2
    NprekN = NpreCell{i}'/1e3;
    NtestkN = NtestArray'/1e3;
    hdata = plot(axes1, NprekN, NtestkN, 'marker', datamarker{i}, 'markersize', 4,...
        'Color', datacolor{i}, 'LineStyle', 'none');
    % linear regression
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [-Inf 0];
    opts.Upper = [Inf 0];
    [cf, stats] = fit(NprekN, NtestkN, 'poly1', opts);
    hreg = refline(cf.p1, cf.p2); set(hreg, 'Color', datacolor{i}, 'LineStyle', regstyle{i});
    hdataArray(i) = hdata;
    hregArray(i) = hreg;
end
axis([200 1200 200 1200])
href = refline(1, 0); set(href, 'Color','r', 'LineStyle', '-');
xlabel('Predicted capacity, N_{u,pre} (kN)')
ylabel('Test capacity, N_{u,exp} (kN)')
lgd = legend([hdataArray(1), hdataArray(2), hregArray(1), hregArray(2), href], ...
    {'Section analysis', 'Jiang & Teng (2013)',...
     'Regression (section)', 'Regression (model)',...
     'N_{u,exp} = N_{u,pre}'},...
    'Location', 'NorthWest');
postfigs(fig1, 6.5/2, false, ftsize)

%% postprocessin of model error data (section)
meN = NtestArray'./NsectionArray';
fprintf(strcat('N model error mean = %.5f\n'), mean(meN));
fprintf(strcat('N model error std = %.5f\n'), std(meN));
fprintf(strcat('N model error cov = %.5f\n'), std(meN)/mean(meN));

[ecdfmeN, meNunique] = ecdf(meN);
%normal
[normEsts(1), normEsts(2) ,normCIs] = normfit(meN, 0.05);
[H, p, ksstat] = kstest(meN,'CDF',[meN,normcdf(meN,normEsts(1), normEsts(2))]);
fcdfmeN = normcdf(meNunique,normEsts(1), normEsts(2));
iae = sum(abs(ecdfmeN - fcdfmeN))/sum(ecdfmeN);
fprintf('normal\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);
%lognormal
[lognEsts,lognCIs] = lognfit(meN, 0.05);
[H, p, ksstat] = kstest(meN,'CDF',[meN,logncdf(meN,lognEsts(1), lognEsts(2))]);
fcdfmeN = logncdf(meNunique,lognEsts(1), lognEsts(2));
iae = sum(abs(ecdfmeN - fcdfmeN))/sum(ecdfmeN);
fprintf('lognormal\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);
%weibull
% mu = mean(meN); sigma = std(meN); cov = sigma/mu;
% wblparam = fsolve(@(x) [x(1)*gamma(1+1./x(2)) - mu ; x(1).^2 * (gamma(1+2./x(2)) - ( gamma(1+1./x(2)).^2)) - sigma^2],[mu;1.2/cov], optimset('Display','off'));
% startwbl = wblparam';
% [wblEsts,wblCIs] = mle(meN, 'distribution', 'weibull', 'start',startwbl);
[wblEsts,wblCIs] = wblfit(meN, 0.05);
[H, p, ksstat] = kstest(meN,'CDF',[meN,wblcdf(meN,wblEsts(1), wblEsts(2))]);
fcdfmeN = wblcdf(meNunique,wblEsts(1), wblEsts(2));
iae = sum(abs(ecdfmeN - fcdfmeN))/sum(ecdfmeN);
fprintf('Weibull\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);
%Gumbel (Maximum)
[gblEsts,gblCIs] = evfit(-meN, 0.05);
[H, p, ksstat] = kstest(meN,'CDF',[meN,1-evcdf(-meN,gblEsts(1), gblEsts(2))]);
fcdfmeN = 1-evcdf(-meNunique,gblEsts(1), gblEsts(2));
iae = sum(abs(ecdfmeN - fcdfmeN))/sum(ecdfmeN);
fprintf('Gumbel\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);

%% Probabilistic plots (section N)
% CDF
fig2 = figure;
axes2 = axes('Parent',fig2, 'box', 'on'); grid on;
hold(axes2,'all');
stairs(meNunique, ecdfmeN, 'color', 'k', 'linestyle', '--'); hold on;
h = plot(meNunique, normcdf(meNunique, normEsts(1), normEsts(2)),...
    meNunique, logncdf(meNunique, lognEsts(1), lognEsts(2)),...
    meNunique, 1-evcdf(-meNunique, gblEsts(1), gblEsts(2)));
set(h, 'lineWidth', 1);
xlabel('Model error data, \delta_m'); ylabel('CDF');
postfigs(fig2, 6.5/2, false, ftsize);
% PDF
fig3 = figure;
axes3 = axes('Parent',fig3, 'box', 'on'); grid on;
hold(axes3,'all');
[~,BinEdge] = histcounts(meN, 'Normalization', 'pdf');
hLine = histogram(meN, 'Normalization', 'pdf');
set(hLine,'FaceColor','cyan','EdgeColor','k',...
    'LineStyle','-', 'LineWidth',0.5);
hold on;
meNaug = sort([meNunique; BinEdge']);
h = plot(meNaug, normpdf(meNaug, normEsts(1), normEsts(2)),...
    meNaug, lognpdf(meNaug, lognEsts(1), lognEsts(2)),...
    meNaug, evpdf(-meNaug, gblEsts(1), gblEsts(2)));
set(h, 'lineWidth', 1);
xlabel('Model error data, \delta_m'); ylabel('PDF');
postfigs(fig3, 6.5/2, false, ftsize);

%% linear plot (e1)
load ../database/frpdatabase.mat
e0col = 19+1;
fmcol = 25+1;
failmode = cfrpdatabase(:,fmcol);
e0array = cfrpdatabase(failmode==1, e0col);
fig4 = figure;
axes4 = axes('Parent',fig4, 'box', 'on');
hold(axes4,'all');
e1preCell = {e1sectionArray+e0array', e1modelArray+e0array'};
datamarker = {'o', '^'};
regstyle = {'--', '-.'};
datacolor = {[0, 0.45, 0.74], [0,0.5,0]};
axis([0, 50, 0, 50])
for i=1:2
    e1pre = e1preCell{i}';
    e1test = e1testArray'+e0array;
    hdata = plot(axes4, e1pre, e1test, 'marker', datamarker{i}, 'markersize', 4,...
        'Color', datacolor{i}, 'LineStyle', 'none');
    % linear regression
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [-Inf 0];
    opts.Upper = [Inf 0];
    [cf, stats] = fit(e1pre, e1test, 'poly1', opts);
    hreg = refline(cf.p1, cf.p2); set(hreg, 'Color', datacolor{i}, 'LineStyle', regstyle{i});
    hdataArray(i) = hdata;
    hregArray(i) = hreg;
end
axis([0, 50, 0, 50])
href = refline(1, 0); set(href, 'Color','r', 'LineStyle', '-');
xlabel('Predicted extra eccentricity, e_{1,pre} (mm)')
ylabel('Extra eccentricity from tests, e_{1,pre} (mm)')
lgd = legend([hdataArray(1), hdataArray(2), hregArray(1), hregArray(2), href], ...
    {'Section analysis', 'Jiang & Teng (2013)',...
     'Regression (section)', 'Regression (model)',...
     'e_{1,exp} = e_{1,pre}'},...
    'Location', 'Southeast');
postfigs(fig4, 6.5/2, false, ftsize);

%% e1 comparison using different curvature models
hcol = 1+1;
Darray = cfrpdatabase(failmode==1, hcol);
fig41 = figure;
axes41 = axes('Parent',fig41, 'box', 'on');
hold(axes41,'all');
plot(e0array./Darray, e1testArray./e1sectionArray, '^',...
     e0array./Darray, e1testArray./e1modelArray, 'v',...
     e0array./Darray, e1testArray./e1modelxi1Array, 'o');
xlabel('Normalized initial eccentricity, e_0/D');
ylabel('Experimental / predicted deflection');
postfigs(fig41, 6.5/2, false, ftsize);

%% postprocessin of model error data of e1 (section)
mee1 = (e0col+e1testArray)'./(e0col+e1sectionArray)';
% mee1 = e1testArray' ./ e1sectionArray';
fprintf(strcat('e1 model error mean = %.5f\n'), mean(mee1));
fprintf(strcat('e1 model error std = %.5f\n'), std(mee1));
fprintf(strcat('e1 model error cov = %.5f\n'), std(mee1)/mean(mee1));

[ecdfmee1, mee1unique] = ecdf(mee1);
%normal
[normEsts(1), normEsts(2) ,normCIs] = normfit(mee1, 0.05);
[H, p, ksstat] = kstest(mee1,'CDF',[mee1,normcdf(mee1,normEsts(1), normEsts(2))]);
fcdfmee1 = normcdf(mee1unique,normEsts(1), normEsts(2));
iae = sum(abs(ecdfmee1 - fcdfmee1))/sum(ecdfmee1);
fprintf('normal\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);
%lognormal
[lognEsts,lognCIs] = lognfit(mee1, 0.05);
[H, p, ksstat] = kstest(mee1,'CDF',[mee1,logncdf(mee1,lognEsts(1), lognEsts(2))]);
fcdfmee1 = logncdf(mee1unique,lognEsts(1), lognEsts(2));
iae = sum(abs(ecdfmee1 - fcdfmee1))/sum(ecdfmee1);
fprintf('lognormal\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);
%weibull
[wblEsts,wblCIs] = wblfit(mee1, 0.05);
[H, p, ksstat] = kstest(mee1,'CDF',[mee1,wblcdf(mee1,wblEsts(1), wblEsts(2))]);
fcdfmee1 = wblcdf(mee1unique,wblEsts(1), wblEsts(2));
iae = sum(abs(ecdfmee1 - fcdfmee1))/sum(ecdfmee1);
fprintf('Weibull\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);
%Gumbel (Maximum)
[gblEsts,gblCIs] = evfit(-mee1, 0.05);
[H, p, ksstat] = kstest(mee1,'CDF',[mee1,1-evcdf(-mee1,gblEsts(1), gblEsts(2))]);
fcdfmee1 = 1-evcdf(-mee1unique,gblEsts(1), gblEsts(2));
iae = sum(abs(ecdfmee1 - fcdfmee1))/sum(ecdfmee1);
fprintf('Gumbel\n');
fprintf('h=%f, p=%f, kstat=%f, iae=%f\n', H, p, ksstat, iae);

%% Probabilistic plots (section e1)
% CDF
fig5 = figure;
axes5 = axes('Parent',fig5, 'box', 'on'); grid on;
hold(axes5,'all');
stairs(mee1unique, ecdfmee1, 'color', 'k', 'linestyle', '--'); hold on;
h = plot(mee1unique, normcdf(mee1unique, normEsts(1), normEsts(2)),...
    mee1unique, logncdf(mee1unique, lognEsts(1), lognEsts(2)),...
    mee1unique, 1-evcdf(-mee1unique, gblEsts(1), gblEsts(2)));
set(h, 'lineWidth', 1);
xlabel('Model error data, \delta_m'); ylabel('CDF');
postfigs(fig5, 6.5/2, false, ftsize);
% PDF
fig6 = figure;
axes6 = axes('Parent',fig6, 'box', 'on'); grid on;
hold(axes6,'all');
[~,BinEdge] = histcounts(mee1, 'Normalization', 'pdf');
hLine = histogram(mee1, 'Normalization', 'pdf');
set(hLine,'FaceColor','cyan','EdgeColor','k',...
    'LineStyle','-', 'LineWidth',0.5);
hold on;
mee1aug = sort([mee1unique; BinEdge']);
h = plot(mee1aug, normpdf(mee1aug, normEsts(1), normEsts(2)),...
    mee1aug, lognpdf(mee1aug, lognEsts(1), lognEsts(2)),...
    mee1aug, evpdf(-mee1aug, gblEsts(1), gblEsts(2)));
set(h, 'lineWidth', 1);
xlabel('Model error data, \delta_m'); ylabel('PDF');
postfigs(fig6, 6.5/2, false, ftsize);

%% model error correlation
fig7 = figure;
axes7 = axes('Parent',fig7, 'box', 'on'); grid on;
hold(axes7,'all');
plot(meN, mee1, 'o');
xlabel('Section model error data, \delta_s'); 
ylabel('Eccentricity model error data, \delta_e');
postfigs(fig7, 6.5/2, false, ftsize);