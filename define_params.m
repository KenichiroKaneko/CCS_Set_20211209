function [PARAM, CONFIG] = define_params(type)
    %% PARAMとCONFIGを読み込む
    % 実験データexpか数値解solかで入力を変えられるようにしてある
    % 定数以外は自動的に決定したいので、実験データごとに入力ファイルを作らない方針

    CONFIG.DataType = type;
    CONFIG.Reload = 1;
    CONFIG.DevideFlux = 0;
    CONFIG.AutoCCSZPos = 1;
    CONFIG.AutoCCSRPos = 0;
    CONFIG.ShowFig = 1;
    CONFIG.CalibTF = 1;
    CONFIG.CalibPF3 = 1;
    CONFIG.RevSenPos = 1;

    % file direcory
    PARAMNUM = 1;
    % PARAM.utst_file_directory = '/Volumes';
    PARAM.utst_file_directory = './input';
    PARAM.input_file_directory = './input2';
    PARAM.temporary_file_directory = './temp';
    PARAM.output_file_directory = './output';
    PARAM.num_sol_name = '/UTST_numel_2033';
    % PARAM.num_sol_name = '/UTST_numel_1000Ip50';
    % PARAM.num_sol_name = '/UTST_numel_1000Ip55';
    % PARAM.num_sol_name = '/UTST_numel_1200Ip50';
    % PARAM.num_sol_name = '/UTST_numel_0840';
    PARAM.num_sol_name = '/UTST_numel_0800Ip50';
    PARAM.MP_pos = "/CCS_sensor_position_opt.txt";
    PARAM.FL_pos = "/CCS_sensor_position_opt.txt";

    % shot param
    PARAM.shotnum = '010';
    PARAM.shotnum_TF = '002';
    PARAM.date = '211217';
    % PARAM.shotnum = '006';
    % PARAM.shotnum_TF = '001';
    % PARAM.date = '210906';
    PARAM.shotnum_PF3 = '003';
    PARAM.date_calib_PF3 = '211217';
    PARAM.time_CCS = '9600';

    % 使わないセンサー
    PARAM.dead_BZ = [19, 20,27, 32, 39, 40];
    PARAM.dead_FL = [1,6,22,25,27,28,30,31,32];
    PARAM.dead_FL = [];
    PARAM.dead_BZ = [];

    % 打切り項数（０でLcurve法）
    PARAM.KUP = 44;
    PARAM.KUP = 85;
    PARAM.KUP = 70;

    % CCS parameters
    PARAM.IUTST = 5;
    PARAM.GETA_YN = 0;
    PARAM.NONC = 1;
    PARAM.SIGM = 0.00;
    PARAM.SEED = 3;
    PARAM.CCS = 1;
    PARAM.CCS = 2;

    if PARAM.CCS == 0
        PARAM.NE = 0;
    end

    z = 0.13;
    z = [-1, 1] * z;

    if PARAM.CCS == 1
        PARAM.IDECCS = 1; % D型のCCS
        i = 1;
        PARAM.R0(i) = 0.3;
        PARAM.Z0(i) = 0;
        PARAM.RR(i) = 0.06;
        PARAM.CAPPER(i) = 1.8;
        PARAM.NE(i) = 3;
        PARAM.MSEC(i, 1) = 1;
        PARAM.MSEC(i, 2) = 1;
        PARAM.MSEC(i, 3) = 1;
        PARAM.SOU(i) = 5;
    else
        PARAM.IDECCS = 1; % 丸いCCS
        for i = 1:PARAM.CCS
            PARAM.R0(i) = 0.22;
            PARAM.Z0(i) = z(i);
            PARAM.RR(i) = 0.06; % R方向の半径
            PARAM.CAPPER(i) = 3; % 楕円率
            % PARAM.R0(i) = 0.25;
            % PARAM.Z0(i) = z(i);
            % PARAM.RR(i) = 0.06; % R方向の半径
            % PARAM.CAPPER(i) = 3; % 楕円率
            PARAM.NE(i) = 3;
            PARAM.MSEC(i, 1) = 1;
            PARAM.MSEC(i, 2) = 1;
            PARAM.MSEC(i, 3) = 1;
            PARAM.SOU(i) = 5;
        end
    end

    PARAM.IPCONST = 0;
    PARAM.LCOND = 1;
    PARAM.ITSKP = 1;
    PARAM.ITRNC = 1;
    PARAM.IDCN = 1;
    PARAM.KUP0 = 0;
    PARAM.MTS = 0;
    PARAM.LXMAP = 1;
    PARAM.PLOTW = 1;
    PARAM.SERX = 1;
    PARAM.LOCAT = 1;
    PARAM.AAA1 = 0.184884;
    PARAM.BBB1 = 0.387214;
    PARAM.AAA2 = 0.164132;
    PARAM.BBB2 = -0.413749;
    PARAM.statR = 0;
    PARAM.endR = 0.1;
    PARAM.I0 = 1;
    PARAM.TRUTOT = 0.1;

    if type == 'exp'

    elseif type == 'sol'

    end

end