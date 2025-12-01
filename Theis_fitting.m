clear all;
clc;
close all;

datas=["testP30.m","testP90.m"];

for ii=1:2
P30 = load(datas(ii));
if ii==1
    r=30;
else 
    r=90;
end
Q = 788;

t_full = P30(2:end, 1);    
s_obs_full = P30(2:end, 2);

NS_max_globale = -inf;  
Opt_migliore = [];      
Start_Index_ottimale = 1; 
s_final_globale = [];  
t_finale = [];          

max_start_index = 50; 
if length(t_full) < max_start_index
    max_start_index = length(t_full) - 5; 
end

for x_start = 1:max_start_index
    
    t = t_full(x_start:end);
    s_obs = s_obs_full(x_start:end);
    
    iter = 10000; 
    NS_g_locale = -inf; 
    Opt_locale = [];  

    for j = 1:iter
        S = 1E-7 + (1E-3+1E-3)*rand; 
        T = randi([1,2000],1,1);   
        
        W = expint(r^2.*S./(4.*T.*t));
        s_sim = Q./(4.*pi())./T.*W;
        
        NS = 1 - sum((s_sim-s_obs).^2) / sum((s_obs-mean(s_obs)).^2);
        
        if NS > NS_g_locale
            s_final_locale = s_sim;
            Opt_locale = [S, T];
            NS_g_locale = NS;
        end
    end
    
    if NS_g_locale > NS_max_globale
        NS_max_globale = NS_g_locale;
        Opt_migliore = Opt_locale;
        Start_Index_ottimale = x_start;
        s_final_globale = s_final_locale;
        t_finale = t;
    end

end

disp(' ');
disp('==================================================');
disp('GLOBAL OPTIMIZATION RESULTS for' );
disp(datas(ii));
disp('==================================================');
fprintf('Best T: %.2f m2/d\n', Opt_migliore(2));
fprintf('Best S: %e [-]\n', Opt_migliore(1));
fprintf('Best starting index: %d\n', Start_Index_ottimale);
fprintf('Correspinding to initial t: %.4f days\n', t_full(Start_Index_ottimale));
fprintf('Maximum NS: %.4f\n', NS_max_globale);
disp('==================================================');

% --- 6. VISUALIZZAZIONE DEI RISULTATI ---
figure
sgtitle(['Optimal fit (NS: ' num2str(NS_max_globale, '%.4f') ') with t\_start = ' num2str(t_full(Start_Index_ottimale), '%.4f') ' days']);

% Grafico 1: Linear-Linear
subplot (1,3,1)
plot (t_full,s_obs_full,'ko','DisplayName','Observed (Complete)')
hold on
plot (t_finale,s_final_globale,'r-','LineWidth',2,'DisplayName','Simulated (Theis)')
plot (t_full(Start_Index_ottimale), s_obs_full(Start_Index_ottimale), 'b*', 'MarkerSize', 10, 'DisplayName', 'Optimal starting point');
xlabel ('t [days]')
ylabel ('s [meters]')
title('Linear Scale')
legend('Location','best')
grid on

% Grafico 2: Semilogaritmico (Scala Cooper-Jacob)
subplot (1,3,2)
semilogx (t_full,s_obs_full,'ko','DisplayName','Observed (Completo)')
hold on
semilogx (t_finale,s_final_globale,'r-','LineWidth',2,'DisplayName','Simulated (Theis)')
plot (t_full(Start_Index_ottimale), s_obs_full(Start_Index_ottimale), 'b*', 'MarkerSize', 10);
xlabel ('t [days] (Log)')
ylabel ('s [meters]')
title('Semilog scale')
grid on

% Grafico 3: Logaritmico-Logaritmico (Scala Theis)
subplot (1,3,3)
loglog (t_full,s_obs_full,'ko','DisplayName','Observed (Complete)')
hold on
loglog (t_finale,s_final_globale,'r-','LineWidth',2,'DisplayName','Simulated (Theis)')
plot (t_full(Start_Index_ottimale), s_obs_full(Start_Index_ottimale), 'b*', 'MarkerSize', 10);
xlabel ('t [days] (Log)')
ylabel ('s [meters] (Log)')
title('Logarithmic scale')
grid on


%%--- 6. VISUALIZZAZIONE DEI Dati ---
figure
sgtitle(['Qualitative analysis']);

sc=Q/(2*pi*Opt_migliore(2));
tc=r^2*Opt_migliore(1)/Opt_migliore(2);
sd=s_obs_full/sc;
td=t_full/tc;

if ii==1
   s_obs_1=s_obs_full;
   t_obs_1=t_full;
else
   s_obs_2=s_obs_full;
   t_obs_2=t_full;
end


for i = 2:length(t_full)-1
    dsdti(i) =(sd(i+1) - sd(i-1)) / (td(i+1) - td(i-1));
    num1 = log(td(i)/td(i-1)) * sd(i+1);
    den1 = log(td(i+1)/td(i)) * log(td(i+1)/td(i-1));
    num2 = log((td(i+1)*td(i-1))/td(i)^2) * sd(i);
    den2 = log(td(i+1)/td(i)) * log(td(i)/td(i-1));
    num3 = log(td(i+1)/td(i)) * sd(i-1);
    den3 = log(td(i)/td(i-1)) * log(td(i+1)/td(i-1));

    dsdti_log(i) = (num1/den1 + num2/den2 - num3/den3);
end

% --- Interpolazione della derivata ---
td_smooth = linspace(min(td(2:end-1)), max(td(2:end-1)), 200);
dsdti_smooth = interp1(td(2:end-1), dsdti(2:end), td_smooth, 'cubic');
dsdti_log_smooth = interp1(td(2:end-1), dsdti_log(2:end), td_smooth, 'cubic');

% --- PLOT ---
subplot(1,3,1)
plot(td, sd, 'ko', 'DisplayName', 'Observed (Complete)')
hold on
plot(td(2:end), dsdti, 'ro', 'DisplayName', 'Derivative (points)')
plot(td_smooth, dsdti_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Derivative (interp.)')
xlabel('t_d [days]')
ylabel('s [scaled]')
title('Linear Scale')
legend('Location','best')
grid on

subplot(1,3,2)
semilogx(td, sd, 'ko', 'DisplayName', 'Observed (Complete)')
hold on
semilogx(td(2:end), dsdti_log, 'ro', 'DisplayName', 'Derivative (points)')
semilogx(td_smooth, dsdti_log_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Derivative (interp.)')
xlabel('t_d [days] (Log)')
ylabel('s [scaled]')
title('Semilog Scale')
legend('Location','best')
grid on

subplot(1,3,3)
loglog(td, sd, 'ko', 'DisplayName', 'Observed (Complete)')
hold on
loglog(td(2:end), dsdti_log, 'ro', 'DisplayName', 'Derivative (points)')
loglog(td_smooth, dsdti_log_smooth, 'r-', 'LineWidth', 2, 'DisplayName', 'Derivative (interp.)')
xlabel('t_d [days] (Log)')
ylabel('s [scaled]')
title('Log-Log Scale')
legend('Location','best')
grid on

end

figure()
loglog(t_obs_1/30^2,s_obs_1)
hold on
loglog(t_obs_2/90^2,s_obs_2)
hold on

%%
% --- Assumo: t_obs_1, s_obs_1, t_obs_2, s_obs_2, r esistono nel workspace
% converti in vettori colonna
t1 = t_obs_1(:);
s1 = s_obs_1(:);
t2 = t_obs_2(:);
s2 = s_obs_2(:);

% r (scaling) già definito; crea le x normalizzate (come nel tuo plot)
x1 = t1./30.^2;
x2 = t2./90.^2;

% Parametri: lunghezza finestra (numero di punti), tolleranza e lunghezza minima run
window_pts = 9;      % finestra di regressione (prova 7-15)
tol = 0.009;          % tolleranza sulla differenza di slope (valore assoluto)
min_run = 3;         % numero minimo di finestre consecutive che devono soddisfare tol

% calcola slope locali
[tC1, slope1] = local_loglog_slope(x1, s1, window_pts);
[tC2, slope2] = local_loglog_slope(x2, s2, window_pts);

% Interpola slope2 sui centri di slope1 (per confronto) — usa extrap NaN per sicurezza
slope2_on_1 = interp1(tC2, slope2, tC1, 'linear', NaN);

% calcola differenza assoluta e trova run continui dove diff < tol
absdiff = abs(slope1 - slope2_on_1);

% Ignora posizioni con NaN (ad es. extrapolazione)
valid_mask = ~isnan(absdiff);

% trova indici dove condizione soddisfatta
good = find(valid_mask & (absdiff <= tol));

% funzione per trovare la prima sequenza continua di lunghezza >= min_run
start_idx = [];
if ~isempty(good)
    groups = diff(good) == 1;              % true quando consecutivi
    % troviamo le sequenze
    seq_starts = [1; find(~groups) + 1];
    seq_ends = [find(~groups); length(good)];
    found = false;
    for gg = 1:length(seq_starts)
        len = seq_ends(gg) - seq_starts(gg) + 1;
        if len >= min_run
            % primo indice globale nella seq
            local_start = good(seq_starts(gg));
            local_end = good(seq_ends(gg));
            start_idx = local_start;
            end_idx = local_end;
            found = true;
            break;
        end
    end
end

if isempty(start_idx)
    fprintf('Nessun intervallo continuo di almeno %d finestre con differenza <= %.3f trovato.\n', min_run, tol);
else
    % tempo di inizio (in unità t/r^2) e pendenza media comune
    t_start = tC1(start_idx);
    t_end = tC1(end_idx);
    slope_common = mean( (slope1(start_idx:end_idx) + slope2_on_1(start_idx:end_idx)) / 2 );
    
    fprintf('Intervallo trovato: da t = %.6g a t = %.6g (unità t/r^2).\n', t_start, t_end);
    fprintf('Primo punto di uguaglianza (start): t = %.6g (t/r^2).\n', t_start);
    fprintf('Pendenza comune (media sul tratto): %.4f  (d log s / d log t)\n', slope_common);
    
    % Plottiamo i risultati
    figure('Name','LogLog original data with start marker');
    loglog(x1, s1, 'b.-'); hold on
    loglog(x2, s2, 'r.-');
    % linea verticale al tempo di inizio
    xline(t_start, '--k', 'LineWidth', 1.5);
    legend('serie 1','serie 2','start equal slope','Location','best');
    xlabel('t / r^2'); ylabel('s');
    title(sprintf('Start equality at t/r^2 = %.6g, slope = %.4f, T= %.4f', t_start, slope_common,Q/(2*pi*slope_common)));
    grid on

    % Plot delle pendenze locali (vs centro finestra) per diagnosticare
    figure('Name','Local slopes (log-log)');
    semilogx(tC1, slope1, 'b.-', 'DisplayName','slope serie 1'); hold on
    semilogx(tC1, slope2_on_1, 'r.-', 'DisplayName','slope serie 2 (interp)');
    % evidenzia il tratto considerato
    semilogx(tC1(start_idx:end_idx), slope1(start_idx:end_idx), 'ko-', 'LineWidth',1.6, 'DisplayName','common region (1)');
    semilogx(tC1(start_idx:end_idx), slope2_on_1(start_idx:end_idx), 'ks--', 'LineWidth',1.6, 'DisplayName','common region (2)');
    xline(t_start, '--k', 'LineWidth', 1.2);
    xlabel('t / r^2'); ylabel('slope = d log s / d log t');
    legend('Location','best'); grid on
end
%%

%% Now that we have a T we can Run montecarlo just to find S


iter = 10000;

NS_g = -1000;


% Computation using Cooper - Jacob Method

for i=1:2
    if i==1
    s_obs=s_obs_1(end-7:end-1);
    t=t_obs_1(end-7:end-1);
    r=30;
    else
    s_obs=s_obs_2(end-9:end);
    t=t_obs_2(end-9:end);
    r=90;
    end
NS_g = -1000;
for j = 1:1:iter

S = 1E-7 + (1E-3+1E-3)*rand;

T = Q/(2*pi*slope_common);

% W = log(2.25.*T.*t./r.^2./S);
W = expint(r^2.*S./(4.*T.*t));

s_sim = Q./4./pi()./T.*W;

NS_1 = s_sim-s_obs;

NS_2 = s_obs-mean(s_obs);

NS = 1 - sum(NS_1.^2)./sum(NS_2.^2);

if NS >= NS_g

s_final = s_sim;

Opt = [S,T];

NS_g = 1 - sum((s_final-s_obs).^2)./sum((s_obs-mean(s_obs)).^2);

end

end



figure


subplot (1,3,1)

plot (t,s_obs,'o')

hold on

plot (t,s_final,'r*')

xlabel ('t [días]')

ylabel ('s [metros]')


subplot (1,3,2)

semilogx (t,s_obs,'o')

hold on

semilogy (t,s_final,'r*')

xlabel ('t [días]')

ylabel ('s [metros]')



subplot (1,3,3)

loglog (t,s_obs,'o')

hold on

loglog (t,s_final,'r*')

xlabel ('t [días]')

ylabel ('s [metros]')

disp(S)

end

%%
% Funzione che calcola slope locale su finestre mobili (log-log)
function [t_centers, slopes] = local_loglog_slope(x, s, window_pts)
    N = length(x);
    if N < window_pts
        error('Numero di punti inferiore alla finestra richiesta.');
    end
    slopes = nan(N - window_pts + 1, 1);
    t_centers = nan(N - window_pts + 1, 1);
    for k = 1:(N - window_pts + 1)
        idx = k:(k + window_pts - 1);
        lx = log(x(idx));
        ls = log(s(idx));
        p = polyfit(lx, ls, 1);        % p(1) = slope in log-log (d log s / d log t)
        slopes(k) = p(1);
        % centro della finestra come valore geometrico (utile per logscale)
        t_centers(k) = exp(mean(lx));
    end
end

