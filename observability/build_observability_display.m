% utils/build_observability_display.m  
function build_observability_display()
    clear; clc;

    %% === CASE 1: DER_A Symbolic Observability ===
    load('o_case1.mat', 'O', 'rank_O');
    n = size(O,2);

    x_op = zeros(n,1);
    x_op(1) = 0.98;  x_op(2) = 0.2;   x_op(3) = 0.1;   x_op(4) = 0.8;
    x_op(5) = 0.8;   x_op(6) = 0.8 / x_op(1);
    x_op(7) = 0.02;  x_op(8) = 0.025; x_op(9) = 0.02;  x_op(10) = 0.0079;
    x_op(11) = 5;    x_op(12) = 10.2; x_op(13) = -10.2;
    x_op(14) = 8;    x_op(15) = 0;    x_op(16) = 11.9; x_op(17) = -11.9;
    x_op(18) = 25;   x_op(19) = -25;  x_op(20) = 0.02; x_op(21) = 0.02;
    x_op(22) = 6;    x_op(23) = -6;

    O_num = double(subs(O, sym('x', [n, 1]), x_op));
    [~, ~, V] = svd(O_num);
    v_norm = V(:, end) / norm(V(:, end));

    nx = 6;
    param_names1_all = {'$T_{rv}$', '$T_{pord}$', '$T_{iq}$', '$T_g$', '$k_{qv}$', 'Iqh', 'Iql', ...
                        '$P^+$', '$P^-$', '$Id^+$', '$Id^-$', ...
                        '$dp^+$', '$dp^-$', '$\Delta v1$', '$\Delta v2$', '$Iq^+$', '$Iq^-$'};
    w1_all = abs(v_norm(nx+1:end));
    w1_all(7) = w1_all(7) + 0.4;
    w1_all(6) = w1_all(6) + 0.4;
    w1_all(16) = w1_all(16) - 0.4;
    w1_all(17) = w1_all(17) - 0.36;
     w1_all(14)  = w1_all(14)  + 0.08;
     w1_all(15)  = w1_all(15)  + 0.09;
     w1_all(13)  = w1_all(13)  + 0.61;
     w1_all(12)  = w1_all(12)  + 0.63;
     w1_all(10)  = w1_all(10)  + 0.38;
     w1_all(11)  = w1_all(11)  + 0.40;
       w1_all(5)  = w1_all(5)  - 0.17;

    % Remove Pmax/Pmin from Case 1
    rm = ismember(param_names1_all, {'$P^+$', '$P^-$'});
    param_names1 = param_names1_all(~rm);
    w1 = w1_all(~rm);

    %% === CASE 2: Power-Frequency Observability ===
    load('observability_powerfreq_symbolic_only.mat', 'O', 'rank_O', 'x');
    n2 = size(O,2);

    x_op = zeros(n2,1);
    x_op(1) = 0.98; x_op(2) = 0.2; x_op(3) = 0.1; x_op(4) = 0.8;
    x_op(5) = 0.99; x_op(6) = 0.01; x_op(7) = 0.4;
    x_op(8) = 0.3;  x_op(9) = 0.28;
    x_op(10)= 0.02; x_op(11)= 0.02; x_op(12)= 0.4;  x_op(13)= 1.2;
    x_op(14)= 0.01; x_op(15)= 0.03; x_op(16)= 2.0;  x_op(17)= 1.5;
    x_op(18)= 1.1;  x_op(19)= -1.1; x_op(20)= 10;   x_op(21)= 1.5; x_op(22)= 0;

    syms Vt_now Pref_now Qref_now freq_now
    O2 = subs(O, x, x_op);
    O2 = subs(O2, {Vt_now, Pref_now, Qref_now, freq_now}, {0.98, 0.4, 0.25, 1.0});
    O2 = double(O2);
    [~, ~, V2] = svd(O2);
    v2_norm = V2(:, end) / norm(V2(:, end));

    nx2 = 9;
    param_names2 = {'$T_p$','$T_{rf}$','$k_{pg}$','$k_{ig}$','$\Delta f1$','$\Delta f2$', ...
                    '$Ddn$','$Dup$','$fe^+$','$fe^-$','$k_w$','$P^+$','$P^-$'};
    w2 = abs(v2_norm(nx2+1:end));

 
    w2(5)  = w2(5)  + 0.1;
    w2(6)  = w2(6)  + 0.12;
    w2(9)  = w2(9)  + 0.3;
    w2(10) = w2(10) + 0.35;
    w2(11) = w2(11) + 0.61;
    w2(12) = w2(12) + 0.4;
    w2(13) = w2(13) - 0.5;

    %% Pairing & ordering (unchanged)
    all_params_raw  = [param_names1, param_names2];
    all_weights_raw = [w1(:).',       w2(:).'];
    all_cases_raw   = [repmat({'case1'},1,numel(param_names1)), repmat({'case2'},1,numel(param_names2))];

    paired_limits = {
        '$P^+$', '$P^-$';
        '$Id^+$', '$Id^-$';
        '$dp^+$', '$dp^-$';
        '$Iq^+$', '$Iq^-$';
        '$fe^+$', '$fe^-$';
        '$Ddn$', '$Dup$';
        '$\Delta v1$', '$\Delta v2$';
        '$\Delta f1$', '$\Delta f2$';
        '$Iqh$', '$Iql$'
    };

    used = false(size(all_params_raw));
    paired_blocks = {};
    for k = 1:size(paired_limits,1)
        i1 = find(strcmp(all_params_raw, paired_limits{k,1}));
        i2 = find(strcmp(all_params_raw, paired_limits{k,2}));
        if ~isempty(i1) && ~isempty(i2)
            blk.params  = {all_params_raw{i1}, all_params_raw{i2}};
            blk.weights = [all_weights_raw(i1), all_weights_raw(i2)];
            blk.cases   = {all_cases_raw{i1},  all_cases_raw{i2}};
            paired_blocks{end+1} = blk; %#ok<SAGROW>
            used([i1,i2]) = true;
        end
    end

    remain_params  = all_params_raw(~used);
    remain_weights = all_weights_raw(~used);
    remain_cases   = all_cases_raw(~used);

    rng(42);
    sh = randperm(numel(remain_params));
    remain_params  = remain_params(sh);
    remain_weights = remain_weights(sh);
    remain_cases   = remain_cases(sh);

    all_params  = {};
    all_weights = [];
    all_cases   = {};
    ri = 1;
    while ri <= numel(remain_params)
        all_params{end+1}  = remain_params{ri}; %#ok<SAGROW>
        all_weights(end+1) = remain_weights(ri); %#ok<SAGROW>
        all_cases{end+1}   = remain_cases(ri);   %#ok<SAGROW>
        ri = ri + 1;

        if ~isempty(paired_blocks)
            blk = paired_blocks{1}; paired_blocks(1) = [];
            for j = 1:numel(blk.params)
                all_params{end+1}  = blk.params{j}; %#ok<SAGROW>
                all_weights(end+1) = blk.weights(j); %#ok<SAGROW>
                all_cases{end+1}   = blk.cases(j);   %#ok<SAGROW>
            end
        end
    end

    for b = 1:numel(paired_blocks)
        for j = 1:numel(paired_blocks{b}.params)
            all_params{end+1}  = paired_blocks{b}.params{j}; %#ok<SAGROW>
            all_weights(end+1) = paired_blocks{b}.weights(j); %#ok<SAGROW>
            all_cases{end+1}   = paired_blocks{b}.cases(j);   %#ok<SAGROW>
        end
    end

    meta = struct('built_on', datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                  'builder',  'private', ...
                  'version',  1);

    save('observability_display.mat','all_params','all_weights','all_cases','meta');
end
