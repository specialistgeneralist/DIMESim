function [eta,fit] = predict_time(n,runT,T)

% PREDICT_TIME Intelligent estimation of experiment completion.
%     [ETA,FIT] = PREDICT_TIME(N,RUNT,T) for a given experiment parameter 
%     table (each row contains a unique experiment, as produced by 
%     MAKE_EXPS_TABLE), predict_time will fit an ensemble logistic regression 
%     trees model to the N completed rows, given a [N,1] vector of completed 
%     run-times. The model is then used to estimate, given remaining 
%     experimental rows (conditions) the likely row-wise run-time, to give a 
%     final ETA (Matlab datenum) and an estimate of the FIT of the model used 
%     (good fit --> 1).
%
%     If the statistics_toolbox is not installed, a simple average of the 
%     completed run-times is used instead to make the prediction.
%
% See also FITRENSEMBLE PREDICT RESUBLOSS

eta = [];
fit = [];

if license('checkout', 'statistics_toolbox')
	% .. train
	Y = runT;
	XTbl = T(1:n,:);
    try
        Mdl = fitrensemble(XTbl, Y, 'NumLearningCycles', 20);
        fit = sprintf('%.2f', 1-resubLoss(Mdl));
	
        % .. predict
        XTbl_left = T(n+1:end,:);
        pred_runT = predict(Mdl, XTbl_left);
    catch
        pred_runT = repmat(mean(runT), height(T)-n, 1);
        fit = 'fit NA, check docs';
    end
else
    pred_runT = repmat(mean(runT), height(T)-n, 1);
    fit = 'fit NA, check docs';
end

% .. output
togo = sum(pred_runT);
eta = datenum(now) + togo;
