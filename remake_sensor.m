function [FF, AA, SENSOR_TPRB, SENSOR_FLXLP] = remake_sensor(PARAM, FF, AA, SENSOR_FLXLP, SENSOR_TPRB)
    whos
    % 指定されたセンサーをFFから消去する
    dead_BZ = PARAM.dead_BZ;
    dead_FL = PARAM.dead_FL + SENSOR_TPRB.NUM;

    % FF行列から削除
    FF(dead_BZ) = [];
    FF(dead_FL) = [];

    % SENSORから削除
    i = PARAM.dead_BZ;
    SENSOR_TPRB.R(i) = [];
    SENSOR_TPRB.Z(i) = [];
    SENSOR_TPRB.TET(i) = [];
    SENSOR_TPRB.TPRB(i) = [];
    SENSOR_TPRB.NUM = length(SENSOR_TPRB.TPRB);
    
    i = PARAM.dead_FL;
    SENSOR_FLXLP.R(i) = [];
    SENSOR_FLXLP.Z(i) = [];
    SENSOR_FLXLP.TET(i) = [];
    SENSOR_FLXLP.FLXLP(i) = [];
    SENSOR_FLXLP.NUM = length(SENSOR_FLXLP.FLXLP);

    % AA行列から削除
    AA(dead_BZ, :) = [];
    AA(dead_FL, :) = [];
end