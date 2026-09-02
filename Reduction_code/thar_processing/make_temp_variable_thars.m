function [temp1, temp2, temp3] = make_temp_variable_thars(info,s)
eval(['temp1=info.' info.list{s} '.HERCEXPT;'])
eval(['temp2=info.' info.list{s} '.OBJECT;'])
eval(['temp3=info.' info.list{s} '.EXPTIME;'])