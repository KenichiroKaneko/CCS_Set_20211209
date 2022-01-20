close all;
load('vars_exp')
CONFIG.ShowFig = 1;
FFout = tsvd(PARAM, CONFIG, AA, FF, FC, ...
        SENSOR_TPRB, SENSOR_NPRB, SENSOR_FLXLP, CCSDAT, WALL);