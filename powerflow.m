clear; clc;

% References
% Nelson Kagan
% William Kersting

% Choose the case system
ieee30;

tol = 1e-6;          % Tol
Sbase = 100;         % MVA
load_factor = 1;     
qg_lim = 1;          % 0 disabled, 1 enabled
ifr = line(:,1);    % from bus
ito = line(:,2);    % to bus
nl  = size(line,1); % number of lines
nb  = size(bus,1); % number of buses

% Identify PV, PQ and the slack bus
tipo = bus(:,2);   
pq  = find(tipo == 1);
pv  = find(tipo == 2);
ref = find(tipo == 3);

% Sort the vector
bp  = sort([pv; pq]);
bq  = pq;
bv = find(tipo == 2);

% Admitance matrix
zl = line(:,3) + 1i * line(:,4);
yl = ones(size(zl, 1), 1) ./ zl;
b = 1i * line(:,5) ./ 2;
a = line(:,6);             % transformer
shift = line(:,7);         % phase transformer
tap = a .* exp(1i .* pi/180 .* shift);
Gs = bus(:,9) ./ Sbase;
Bs = bus(:,10) ./ Sbase;

ysh = sparse(Gs + 1j * Bs); % shunt admitance

Y = zeros(nb,nb);

for idx = 1 : nl
    i = ifr(idx);
    k = ito(idx);
    
    Y(i, i) = Y(i, i) + yl(idx) / (tap(idx) .* conj(tap(idx))) + b(idx);
    Y(k, k) = Y(k, k) + yl(idx) + b(idx);
    Y(i, k) = Y(i, k) - yl(idx) / tap(idx);
    Y(k, i) = Y(k, i) - yl(idx) / conj(tap(idx));
end

Y = Y + diag(ysh);

% Newton-raphson
% Start the variables
delta   = zeros(nb,1);
vo      = ones(nb, 1);

% Controled voltage
vo(pv)  = bus(pv, 3);
vo(ref) = bus(ref,3);

% Power variables
Pg = bus(:,5);
Qg = bus(:,6);
Pl = bus(:,7);
Ql = bus(:,8);
Qgmax = bus(:,12) ./ Sbase;
Qgmin = bus(:,11) ./ Sbase;

Pg = Pg ./ Sbase;
Qg = Qg ./ Sbase;
Pl = Pl .* load_factor ./ Sbase;
Ql = Ql .* load_factor ./ Sbase;

% Control variables
flag_error = 0; % error flag
error = 0;      % power load error
idx = 0;
iter_wait = 0;  % iterator to convert a bus
iter_conv = 0;  % iterator to check convergence
 
% Start Newton-Rapshon
while true
    
    % Check if the reactive power went over the limit
    % Convert the buses where this violation occurred. If neccessary, change the slack bus
    if idx > 2 && qg_lim
        
        f = 0;
        if ~isempty(bv)
            for i = 1:length(bv)
                if f == 0
                    if Qg(bv(i)) > Qgmax(bv(i))
                        Qg(bv(i)) = Qgmax(bv(i));
                        f = bv(i);
                        tipo(f) = 1;
                        disp('================================================================================');
                        fprintf('  Bus %d changed to PQ \n',f);
                        disp('================================================================================');
                    end
                    if Qg(bv(i)) < Qgmin(bv(i))
                        Qg(bv(i)) = Qgmin(bv(i));
                        f = bv(i);
                        tipo(f) = 1;
                        disp('================================================================================');
                        fprintf('  Bus %d changed to PQ \n',f);
                        disp('================================================================================');
                    end
                end
            end
            
            if f > 0
                bv = find(tipo == 2);
                bq = find(tipo == 1);
            end
        end
        
        % Check and change the slack bus, if neccessary
        if f == 0
                % Upper limit
                if Qg(ref) > Qgmax(ref) && iter_wait > 4
                    
                    if isempty(bv)
                        disp('================================================================================');
                        fprintf(' Load factor: %4.3f \n', load_factor);
                        disp('================================================================================');
                        disp('================================================================================');
                        fprintf('  The power system does not converge \n');
                        fprintf('  There is not enough reactive power \n');
                        disp('================================================================================');
                        flag_error = 1;
                        break;
                    else
                        
                        Qg(ref) = Qgmax(ref);
                        f = ref;
                        tipo(f) = 1;
                        disp('================================================================================');
                        fprintf('  Bus %d changed to PQ \n',f);
                        disp('================================================================================');
                        
                        ref = bv(1);
                        tipo(bv(1)) = 3;
                        bv = find(tipo == 2);
                        bq = find(tipo == 1);
                        bp = find(tipo ~= 3);
                        fprintf('  Bus %d is the new slack bus \n',ref);
                        disp('================================================================================');
                        
                        iter_wait = 0;
                    end
                    
                end
                
                % Lower limit
                if Qg(ref) < Qgmin(ref) && iter_wait > 4
                    
                    if isempty(bv)
                        disp('================================================================================');
                        fprintf(' Load factor: %4.3f \n', load_factor);
                        disp('================================================================================');
                        disp('================================================================================');
                        fprintf('  The power system does not converge \n');
                        fprintf('  There is not enough reactive power \n');
                        disp('================================================================================');
                        flag_error = 1;
                        break;
                    else
                        
                        Qg(ref) = Qgmin(ref);
                        tipo(ref) = 1;
                        disp('================================================================================');
                        fprintf('  Bus %d changed to PQ \n',ref);
                        disp('================================================================================');
                        
                        ref = bv(1);
                        tipo(bv(1)) = 3;
                        bv = find(tipo == 2);
                        bq = find(tipo == 1);
                        bp = find(tipo ~= 3);
                        fprintf('  Bus %d is the new slack bus \n',ref);
                        disp('================================================================================');
                        
                        iter_wait = 0;
                    end
                
                end
        end
                
        
    end
    
    % Load balance
    v = vo .* exp(1i * delta);
    s = v .* conj(Y * v);
    m = [(Pg(bp) - Pl(bp) - real(s(bp)));
         (Qg(bq) - Ql(bq) - imag(s(bq)))];
     
    Qg(bv) = Ql(bv) + imag(s(bv));
    Qg(ref) = Ql(ref) + imag(s(ref));
    Pg(ref) = real(s(ref)) + Pl(ref);
    
    derror = max(abs(m)) - error;
    error = max(abs(m));
    
    % Check error
    if error < tol
        break;
        
    % Check convergence   
    elseif derror > 0
        iter_conv = iter_conv + 1;
        
        if iter_conv > 10
            flag_error = 1;
            disp('================================================================================');
            fprintf(' Load factor: %4.3f \n', load_factor);
            disp('================================================================================');
            disp('================================================================================');
            fprintf('  The power system did not converged in 10 iterations \n');
            disp('================================================================================');
            break;
        end
        
    elseif derror < 0
        if error > 1/tol
            iter_conv = iter_conv + 1;
        else
            iter_conv = 0;
        end
        
    end
    
    % Jacobian
    D_delta = diag(v) * conj(Y * 1i * diag(v)) + diag(conj(Y * v)) * 1i * diag(v);
    D_v     = diag(v) * conj(Y * diag(v ./ vo) + diag( conj(Y * v) ) * diag(v ./ vo));
    D_p_delta = real(D_delta(bp, bp));
    D_q_delta = imag(D_delta(bq, bp));
    D_p_v = real(D_v(bp, bq));
    D_q_v = imag(D_v(bq, bq));
    J = [D_p_delta, D_p_v;
         D_q_delta, D_q_v];
    dx = J \ m;
    
    % Variables delta increase
    delta(bp) = delta(bp) + dx(1:length(bp));
    vo(bq) = vo(bq) + dx(length(bp)+1:length(dx));
    
    idx = idx + 1;
    iter_wait = iter_wait + 1;
    
end

% pu to MW and MVAr
Pg = Pg .* Sbase;
Qg = Qg .* Sbase;
Pl = Pl .* Sbase;
Ql = Ql .* Sbase;

% Restore the bus type vector
tipo = bus(:,2);

% Fix the angle, if the slack bus has changed
delta = delta - delta(1);

%% Solution
if flag_error == 1
    disp('================================================================================');
    fprintf(' END \n');
    disp('================================================================================');
else
    disp('================================================================================');
    fprintf(' Load factor: %4.3f \n', load_factor);
    disp('================================================================================');
    disp('   Bus  ');
    disp('================================================================================');
    fprintf(' Bus \t Voltage \t\t\t Generation \t\t Load \n');
    fprintf(' num \t Mag(pu)  Ang(degrees) \t\t P(MW)  Q(MVAr) \t P(MW) \t Q(MVAR) \n');
    disp('--------------------------------------------------------------------------------');
    for m=1:max(bus(:,1))
        fprintf(' %3g   ' ,m);
        fprintf(' %6.3f \t' ,vo(m));
        fprintf(' %6.3f \t' ,delta(m) * 180 / pi);
        
        if tipo(m) == 2 || tipo(m) == 3
            fprintf('\t %4.2f\t' ,Pg(m));
            fprintf('%4.2f \t\t' ,Qg(m));
        else
            fprintf('\t   -  \t - \t\t');
        end
        
        fprintf('%4.2f \t' ,Pl(m));
        fprintf('%4.2f \t' ,Ql(m));
        fprintf('\n');
    end
    disp('--------------------------------------------------------------------------------');
    fprintf('  \t\t\t\t\t Total \t Total \t\t Total \t Total \n');
    disp('--------------------------------------------------------------------------------');
    fprintf('  \t\t\t\t\t %4.2f %4.2f \t\t %4.2f %4.2f \n', sum(Pg), sum(Qg), sum(Pl), sum(Ql));
    disp('--------------------------------------------------------------------------------');
    fprintf( 'Iterations: %3g \n',idx)
    
    % Short-circuit
    %bus = 2; % Faulty bus
    
    % Load admitance
    %Sl = (Pl + 1i*Ql) ./ Sbase;
    %yl = Sl ./ (abs(vo).^2);
    %Y = Y + diag(yl);
    
    % Generator admitance. Considered 0.01pu
    %yg = zeros(size(yl));
    %yg(tipo == 2) = 1i*0.01; % Bus PV admitance
    %yg(tipo == 3) = 1i*0.01; % slack bus admitance
    %Y = Y + diag(yg);
    
    %Z = inv(sparse(Y));
    %Icc = vo(bus) / Z(bus,bus);
end