function [R,t,rmse] = RANSAC_LeastQuantileOfSquareErrors(pl1,pl2,params)
% this function computes the Rotation (R) and translation (T) matrix from
% points in pl1 to points in pl2 using RANSAC methodology
%
% ------
% INPUTS
% ------
% pl1 and pl2: structure containing
% pl1.Location= 3xN matrix with coordinates of the points from object 1
% pl2.Location= 3xN matrix with coordinates of corresponding poinst from 
%               object 2   
% params: parameter structure containing
%       QuantilePercentage (0.5 = median) = percentage to consider for
%                                           calculating the rmse
%       NbTrialPoints = number of points used in every trial
%       NbTrials = number of trials
% 
% -------
% OUTPUTS
% -------
% R = 3x3 rotation matrix
% T = 3x1 translation matrix
% rmse = smallest rmse of all trials
% 
% Author: Jeroen Hermans
  
    quant = params.QuantilePercentage;
    nbTrialPoints = params.NbTrialPoints;
    nbTrials = params.NbTrials;

    pts1 = pl1.Location;
    pts2 = pl2.Location;

    nbPts = size(pts1,2);

    err = 10000000;
    R = zeros(3);
    t = zeros(3,1);
    i = 0;
    ii = 0;
  
    while i < 2*nbTrials
        i = i+1;
        while ii < nbTrials
            indsSel = randsample(nbPts,nbTrialPoints);
            indsLeft = setdiff(1:nbPts,indsSel);

            pl1Sel.Location = pl1.Location(:,indsSel);
            pl2Sel.Location = pl2.Location(:,indsSel);

            try
                [Rprobe,tprobe] = myprocrustes(pl1Sel,pl2Sel);  
            catch
                break
            end

            pts1Left = pts1(:,indsLeft);  
            pts2Left = pts2(:,indsLeft);    

            pts1LeftT = Rprobe*pts1Left + tprobe*ones(1,length(indsLeft));
            errs = sqrt(sum((pts1LeftT-pts2Left).^2));
            [serrs] = sort(errs,'ascend');
            selInd = round(quant*length(indsLeft));

            if(err>serrs(selInd))
                err = serrs(selInd);
                R = Rprobe;
                t = tprobe;
            end
            ii = ii+1;
        end
    end
    if err == 10000000
        error('Ransac failed to find a good solution')
    end
    pts1T = R*pts1 + t*ones(1,nbPts);
    rmse = sqrt(mean(sum((pts2-pts1T).^2)));