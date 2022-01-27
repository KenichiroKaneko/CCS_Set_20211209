function [PARAM, CONFIG] = define_params(type)
    %% PARAMとCONFIGを読み込む
    % 実験データexpか数値解solかで入力を変えられるようにしてある
    % 定数以外は自動的に決定したいので、実験データごとに入力ファイルを作らない方針

    CONFIG.DataType = type;
    CONFIG.Reload = 1;
    CONFIG.DevideFlux = 1;
    CONFIG.AutoCCSZPos = 0;
    CONFIG.AutoCCSRPos = 0;
    CONFIG.ShowFig = 1;
    CONFIG.ShowFig2 = 0;
    CONFIG.CalibTF = 1;
    CONFIG.CalibPF3 = 0;
    CONFIG.RevSenPos = 1;
    CONFIG.Nmrz = 'norm';

    % file direcory
    PARAMNUM = 1;
    % PARAM.utst_file_directory = '/Volumes';
    PARAM.utst_file_directory = './input';
    PARAM.input_file_directory = './input3';
    PARAM.temporary_file_directory = './temp';
    PARAM.output_file_directory = './output';
    PARAM.num_sol_name = '/UTST_numel_2033'; PARAM.CCS = 1; z = 0;
    % PARAM.num_sol_name = '/UTST_numel_0600Ip50'; PARAM.CCS = 2; % dame
    % PARAM.num_sol_name = '/UTST_numel_0640Ip50'; PARAM.CCS = 2; % dame
    % PARAM.num_sol_name = '/UTST_numel_0720Ip50'; PARAM.CCS = 2; z = 0.48; p = [0.2, 0.045, 3];
    % PARAM.num_sol_name = '/UTST_numel_0760Ip50'; PARAM.CCS = 2; z = 0.44; p = [0.22, 0.045, 3];
    % PARAM.num_sol_name = '/UTST_numel_0800Ip50'; PARAM.CCS = 2; z = 0.42; p = [0.22, 0.06, 3];
    % PARAM.num_sol_name = '/UTST_numel_0840Ip50'; PARAM.CCS = 2; z = 0.39; p = [0.22, 0.058, 3];
    % PARAM.num_sol_name = '/UTST_numel_0880Ip50'; PARAM.CCS = 2; z = 0.355; p = [0.22, 0.055, 3];
    % PARAM.num_sol_name = '/UTST_numel_0920Ip50'; PARAM.CCS = 2; z = 0.32; p = [0.22, 0.057, 3];
    % PARAM.num_sol_name = '/UTST_numel_1000Ip50'; PARAM.CCS = 2; z = 0.253; p = [0.25, 0.06, 3];
    % PARAM.num_sol_name = '/UTST_numel_1200Ip50'; % dame
    PARAM.MP_pos = "/CCS_sensor_position_opt.txt";
    PARAM.FL_pos = "/CCS_sensor_position_opt.txt";
    % PARAM.MP_pos = "/CCS_MP_sensor_position_i.txt";
    % PARAM.FL_pos = "/CCS_FLXLP_sensor_position_i.txt";

    % shot param
    % PARAM.shotnum = '010';
    % PARAM.shotnum_TF = '002';
    % PARAM.date = '211217';

    PARAM.shotnum = '004';
    PARAM.shotnum_TF = '002';
    PARAM.date = '211217';

    PARAM.shotnum_PF3 = '003';
    PARAM.date_calib_PF3 = '211217';

    PARAM.time_CCS = '9605';

    % 使わないセンサー
    PARAM.dead_BZ = [19, 20,27, 32, 39, 40];
    PARAM.dead_FL = [1,6,22,25,27,28,30,31,32];
    PARAM.dead_FL = [];
    PARAM.dead_BZ = [];

    % 打切り項数（０でLcurve法）
    PARAM.KUP = 44;
    PARAM.KUP = 85;
    PARAM.KUP = 75;

    % CCS parameters
    PARAM.IUTST = 5;
    PARAM.GETA_YN = 0;
    PARAM.NONC = 1;
    PARAM.SIGM = 0.00;
    PARAM.SEED = 3;
    % PARAM.CCS = 1;
    % PARAM.CCS = 2;

    if PARAM.CCS == 0
        PARAM.NE = 0;
    end

    % z1 z2         rr rr capper
    % z = 0.48; p = [0.2, 0.045, 3];
    z = [-1, 1] * z;

    if PARAM.CCS == 1
        PARAM.IDECCS = 1; % D型のCCS
        i = 1;
        PARAM.R0(i) = 0.34;
        PARAM.Z0(i) = 0;
        PARAM.RR(i) = 0.05;
        PARAM.CAPPER(i) = 1.6;
        PARAM.NE(i) = 3;
        PARAM.MSEC(i, 1) = 1;
        PARAM.MSEC(i, 2) = 1;
        PARAM.MSEC(i, 3) = 1;
        PARAM.SOU(i) = 5;
    else
        PARAM.IDECCS = 1; % 丸いCCS
        for i = 1:PARAM.CCS
            PARAM.Z0(i) = z(i);
            PARAM.R0(i) = p(1);
            PARAM.RR(i) = p(2); % R方向の半径
            PARAM.CAPPER(i) = p(3); % 楕円率
            % PARAM.Z0(i) = z(i);
            % PARAM.R0(i) = 0.22;
            % PARAM.RR(i) = 0.045; % R方向の半径
            % PARAM.CAPPER(i) = 3; % 楕円率
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
