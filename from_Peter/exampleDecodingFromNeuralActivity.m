% Example decoding from neural activity
%If 2AFC: only use ZL
%If 2AUC: Use ZL and ZR

% 1) Get behavioural data
D = struct;
D.stimulus = [0 1; 1 0; 1 1; 0 0];
D.response = [2;1;1;3];
D.repeatNum = [1;1;1;1];

% 2) Fit behavioural model
g = GLM(D).setModel('C50-subset').fit;
behavParameterFit = g.parameterFits;

% 3) Do a cross-validated fit for the behavioural model, to provide a
% measure of the baseline behavioural model log likelihood
g = GLM(D).setModel('C50-subset').fitCV(10);
pL = g.p_hat(:,1);    pR = g.p_hat(:,2);    pNG = g.p_hat(:,3);  
likelihood = pL.*(g.data.response==1) + pR.*(g.data.response==2) + pNG.*(g.data.response==3);
loglik_bpt = mean(log2(likelihood));

% 4) Now do a cross-validated fit for the neural activity. Use the behavioural model
% parameter (non-crossval fit) to provide an 'offset' for each trial, which
% reflects the contribution of the behavioural model
D.offset_ZL = g.ZL(behavParameterFit, g.Zinput(g.data));
D.offset_ZR = g.ZR(behavParameterFit, g.Zinput(g.data));
D.neur = [20;-30;23;12];
g = GLM(D).setModel('neur').fitCV(10);
pL = g.p_hat(:,1);    pR = g.p_hat(:,2);    pNG = g.p_hat(:,3);  
likelihood = pL.*(g.data.response==1) + pR.*(g.data.response==2) + pNG.*(g.data.response==3);
loglik_bpt_neur = mean(log2(likelihood));

       