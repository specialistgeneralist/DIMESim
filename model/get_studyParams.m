function [studyParams] = get_studyParams(INFILE)

%% Read in study params, and make available to the model.

eval(fileread(INFILE));     % -- create key study vars in current workspace

% -- make available to model via studyParams struct
studyParams.DIME_grad_means  = DIME_grad_means;
studyParams.DIME_grad_ses    = DIME_grad_ses;

studyParams.D_init_MEAN = D_init_MEAN;
studyParams.I_init_MEAN = I_init_MEAN;
studyParams.M_init_MEAN = M_init_MEAN;
studyParams.E_init_MEAN = E_init_MEAN;

studyParams.D_init_SD = D_init_SD;
studyParams.I_init_SD = I_init_SD;
studyParams.M_init_SD = M_init_SD;
studyParams.E_init_SD = E_init_SD;
