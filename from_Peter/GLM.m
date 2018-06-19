classdef GLM
    properties (Access=public)
        expRef;
        modelString;
        parameterLabels;
        parameterFits;
        parameterBounds;
        parameterStart;
        Zinput;
        ZL;
        ZR;
        regularise;
        data;
        p_hat;
        ContrastDimensions;
    end
    
    properties (Access=public)
        guess_bpt;
        lapseFlag=0;
    end
    
    methods
        function obj = GLM(inputData)
            if isa(inputData,'struct')
                %If input is a struct containing AT MINIMUM the fields:
                %                 -stimulus
                %                 -response
                %                 -repeatNum
                if isfield(inputData,'contrast_cond')
                    inputData.stimulus=inputData.contrast_cond(:,1:2);
                end
                
                obj.data = inputData;
                obj.expRef = 'none';
                
            elseif isa(inputData,'char')
                %if expRef, then load using the dat package
                obj.expRef = inputData;
                
                D = loadData(inputData);
                
                obj.data = D;
            else
                error('GLM:constructorFail', 'Must pass either an expRef or data struct to GLM constructor');
            end
            
            if ~isempty(obj.data.response)
                
                if any(min(obj.data.stimulus,[],2)>0)
                    obj.ContrastDimensions = 2;
                else
                    obj.ContrastDimensions = 1;
                end
                
                tab = tabulate(obj.data.response);
                tab = ( tab(:,3) + 0.5)/100; %Add 0.5 choice to each to regularise
                obj.guess_bpt=sum(tab.*log2(tab));
            else
                obj.ContrastDimensions = NaN;
                obj.guess_bpt = NaN;
            end
        end
        
        function obj = setModel(obj,modelString)
            %             Function which contains the model definitions. Each
            %             definition requires a set of parameter Labels, estimation
            %             bounds, and some anonymous functions. The functions ZL and
            %             ZR define the two linear models in 3-class multinomial
            %             regression for 3-choice behaviour. For 2-choice behaviour,
            %             only ZL needs to be defined for logistic regression. The
            %             Zinput function provides input data to ZL and ZR. Often
            %             this will be simply the trial contrast but it can also be
            %             any other feature of interest.
            %
            %             ZL and ZR always take to inputs P and IN. P is a vector of
            %             parameter values. IN is a matrix of model inputs as derived
            %             from the Zinput function
            %
            %             Zinput always takes in the data struct and outputs a matrix of
            %             variables to be used in the model
            
            obj.modelString = modelString;
            obj.parameterFits = [];
            obj.parameterStart = [];
            
            switch(modelString)
                case 'Offset' %Model guesses based on the proportion of responses in the data
                    %used as a baseline to compare other models
                    obj.parameterLabels = {'Offset_L','Offset_R'};
                    obj.parameterBounds = [-inf -inf; +inf +inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1)*ones(length(in(:,1)),1));
                    obj.ZR = @(P,in)(P(2)*ones(length(in(:,1)),1));
                case 'C-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf;
                        +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1)   );
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2)  );
                case 'C-subset-biasAsContrast'
                    obj.parameterLabels = {'Bias_L','ScaleL_L','Bias_R','ScaleR_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf;
                        +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(2).*(in(:,1) + P(1))   );
                    obj.ZR = @(P,in)( P(4).*(in(:,2) + P(3))  );
                case 'C'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf;
                        +inf +inf +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1) + P(3).*in(:,2));
                    obj.ZR = @(P,in)(P(4) + P(5).*in(:,1) + P(6).*in(:,2));
                case 'C^N'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R','N'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf 0;
                        +inf +inf +inf +inf +inf +inf +inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(7) + P(3).*in(:,2).^P(7));
                    obj.ZR = @(P,in)(P(4) + P(5).*in(:,1).^P(7) + P(6).*in(:,2).^P(7));
                case 'C^N-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0;
                        +inf +inf +inf +inf 3];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(5));
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2).^P(5));
                case 'C^N-subset-biasAsContrast'
                    obj.parameterLabels = {'Bias_L','ScaleL_L','Bias_R','ScaleR_R','N'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0;
                        +inf +inf +inf +inf 3];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(2).*(in(:,1).^P(5) + P(1)) );
                    obj.ZR = @(P,in)( P(4).*(in(:,2).^P(5) + P(3)) );
                case 'C^NL^NR-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N_L','N_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0;
                        +inf +inf +inf +inf 3 3];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*in(:,1).^P(5));
                    obj.ZR = @(P,in)(P(3) + P(4).*in(:,2).^P(6));
                case 'C50'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','ScaleR_L','Offset_R','ScaleL_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf +inf +inf 3 0.8];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(7))./(in(:,1).^P(7) + P(8)^P(7)) + P(3).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));
                    obj.ZR = @(P,in)(P(4) + P(5).*(in(:,1).^P(7))./(in(:,1).^P(7) + P(8)^P(7)) + P(6).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));
                case 'C50-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf 3 0.8];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)) );
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5) + P(6)^P(5)) );
                    
                case 'C50-subset-sharedS'
                    obj.parameterLabels = {'Offset_L','Offset_R','Scale','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf 0 0.001;
                        +inf +inf +inf 3 0.8];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(3).*(in(:,1).^P(4))./(in(:,1).^P(4) + P(5)^P(4)) );
                    obj.ZR = @(P,in)(P(2) + P(3).*(in(:,2).^P(4))./(in(:,2).^P(4) + P(5)^P(4)) );
                case 'C50-subset-neur'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50','neurL','neurR'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001 -inf -inf;
                        +inf +inf +inf +inf 3 0.8 inf inf];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2) D.neur(:,1)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)) + P(7).*in(:,3) );
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5) + P(6)^P(5)) + P(8).*in(:,3) );
                case 'neur'
                    obj.parameterLabels = {'Bias_L','Bias_R','neurL','neurR'};
                    obj.parameterBounds = [-inf -inf -inf -inf;
                        +inf +inf +inf +inf];
                    
                    obj.Zinput = @(D)([D.offset_ZL D.offset_ZR D.neur]);
                    obj.ZL = @(P,in)( in(:,1) + P(1) + P(3).*in(:,3) );
                    obj.ZR = @(P,in)( in(:,2) + P(2) + P(4).*in(:,3) );
                    
                    
                case 'C50-subset-biasAsContrast'
                    obj.parameterLabels = {'Bias_L','ScaleL_L','Bias_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf 3 0.8];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(2).*((in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)) + P(1) ) );
                    obj.ZR = @(P,in)( P(4).*((in(:,2).^P(5))./(in(:,2).^P(5) + P(6)^P(5)) + P(3) ) );
                case 'C50-subset-v1'
                    obj.parameterLabels = {'Offset_L','Alpha_L','ScaleL_L','Offset_R','Alpha_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf +inf +inf 3 0.8];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(1) + P(3).*((in(:,1).^P(7))./(in(:,1).^P(7) + P(8)^P(7)) + P(2) ) );
                    obj.ZR = @(P,in)( P(4) + P(6).*((in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)) + P(5) ) );
                case 'C50-subset-separate'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf;
                        +inf +inf +inf +inf inf inf];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2) D.response]);
                    obj.ZL = @(P,in)( P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)) ).*(in(:,3)~=2);
                    obj.ZR = @(P,in)( P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5) + P(6)^P(5)) ).*(in(:,3)~=1);
                    
                case 'C50-subset-contrastGain'
                    obj.parameterLabels = {'B_L','B_R','Alpha_L','Alpha_R','Sigma_L','Sigma_R','N'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0 1;
                        +inf +inf +inf +inf 1 1 1];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(1) + P(3).*(in(:,1).^P(7))./(in(:,1).^P(7) + P(5)) );
                    obj.ZR = @(P,in)( P(2) + P(4).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(6)) );
                case 'C50LR-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N_L','C50_L', 'N_R', 'C50_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001 0 0.001;
                        +inf +inf +inf +inf 3 0.8 3 0.8];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)));
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(7))./(in(:,2).^P(7) + P(8)^P(7)));
                case 'Clogistic-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','k','c0_L','c0_R'};
                    obj.parameterBounds = [-inf -inf -inf -inf -inf -inf -inf;
                        +inf +inf +inf +inf +inf +inf +inf];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2)./(1 + exp(-P(5).*(in(:,1) - P(6)))) );
                    obj.ZR = @(P,in)(P(3) + P(4)./(1 + exp(-P(5).*(in(:,2) - P(7)))) );
                case 'Supersaturation-subset'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50','Q'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0 0;
                        +inf +inf +inf +inf +inf +inf 2];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5)*P(7) + P(6)^P(5)*P(7)));
                    obj.ZR = @(P,in)(P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5)*P(7) + P(6)^P(5)*P(7)));
                case 'AFC'
                    obj.parameterLabels = {'Offset','ScaleL','ScaleR'};
                    obj.parameterBounds = [-inf -inf -inf; +inf +inf +inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2)*in(:,1) + P(3)*in(:,2));
                    obj.ZR = [];
                case 'AFC_diff'
                    obj.parameterLabels = {'Offset','Scale'};
                    obj.parameterBounds = [-inf -inf;inf inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(1) + P(2)*(in(:,1)-in(:,2)) );
                    obj.ZR = [];
                case 'AFC_diff_alternate'
                    obj.parameterLabels = {'Offset','Scale'};
                    obj.parameterBounds = [-inf -inf;inf inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)(P(2)*P(1) + P(2)*(in(:,1)-in(:,2)) );
                    obj.ZR = [];
                    %                 case 'AFC_multimodalDiff'
                    %                     obj.parameterLabels = {'OffsetV','OffsetA','ScaleV','ScaleA'};
                    %                     obj.parameterBounds = [-inf -inf -inf -inf ;inf inf inf inf];
                    %                     obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2) D.auditory_cond(:,1) D.auditory_cond(:,2)]);
                    %                     obj.ZL = @(P,in)(P(1) + P(2)*(in(:,1)-in(:,2)) );
                    %                     obj.ZR = [];
                case 'C^N-subset-hist-response'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','HistL_L','HistR_L','HistNG_L','HistL_R','HistR_R','HistNG_R','N'};
                    obj.parameterBounds = [-inf(1,10) 0;
                        +inf(1,10) inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2) D.hist(:,1) D.hist(:,2) D.hist(:,3)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(11) + P(5).*in(:,3) + P(6).*in(:,4) + P(7).*in(:,5) );
                    obj.ZR = @(P,in)( P(3) + P(4).*in(:,2).^P(11) + P(8).*in(:,3) + P(9).*in(:,4) + P(10).*in(:,5) );
                case 'C^N-subset-hist-contrast'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','HistL_L','HistR_L','HistNG_L','HistL_R','HistR_R','HistNG_R','N'};
                    obj.parameterBounds = [-inf(1,10) 0;
                        +inf(1,10) inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2) D.hist(:,1) D.hist(:,2) D.hist(:,3)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(11) + P(5).*in(:,3) + P(6).*in(:,4) + P(7).*in(:,5) );
                    obj.ZR = @(P,in)( P(3) + P(4).*in(:,2).^P(11) + P(8).*in(:,3) + P(9).*in(:,4) + P(10).*in(:,5) );
                case 'C^N-subset-hist-success'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','Hist_L','Hist_R','N'};
                    obj.parameterBounds = [-inf(1,6) 0;
                        +inf(1,6) inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2) D.hist]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(7) + P(5).*in(:,3) );
                    obj.ZR = @(P,in)( P(3) + P(4).*in(:,2).^P(7) + P(6).*in(:,3) );
                case 'AFC-C^N-subset-hist-success'
                    obj.parameterLabels = {'Offset','ScaleL_L','ScaleR_R','Hist','N'};
                    obj.parameterBounds = [-inf(1,4) 0;
                        +inf(1,4) inf];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2) D.hist]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(5) + P(3).*in(:,2).^P(5) + P(4).*in(:,3) );
                case 'C^N-subset-Qlearning'
                    obj.parameterLabels = {'ScaleL_L','ScaleR_R','Q_L','Q_R','N'};
                    obj.parameterBounds = [-inf(1,2) 1 1 0.1;
                        inf(1,2) 1 1 1];
                    obj.Zinput = @(D,B)(behav.GLM_QlearningInputs(D,B));
                    obj.ZL = @(P,in)(P(1).*in(:,1).^P(5) + P(3).*in(:,3) ); %+ P(5).*in(:,5));
                    obj.ZR = @(P,in)(P(2).*in(:,2).^P(5) + P(4).*in(:,4) ); %+ P(5).*in(:,5));
                    obj.parameterStart = @()([-1+2*rand(1,4) 0.5]);
                    
                case 'C^N-subset-2AFC'
                    obj.parameterLabels = {'Offset','ScaleL','ScaleR','N'};
                    obj.parameterBounds = [-inf -inf -inf 0.1;
                        +inf +inf +inf 1];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1).^P(4) + P(3).*in(:,2).^P(4)  );
                    obj.parameterStart = [0 0 0 0.1];
                    
                case 'C50-subset-2AFC'
                    obj.parameterLabels = {'Offset','ScaleL','ScaleR','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf 0 0.01;
                        +inf +inf +inf 3 8];
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*(in(:,1).^P(4))./(in(:,1).^P(4) + P(5)^P(4)) + P(3).*(in(:,2).^P(4))./(in(:,2).^P(4) + P(5)^P(4))  );
                    obj.parameterStart = [0 0 0 0.1 0.1];
                    
                case 'C^N-subset-Qlearning-2AFC'
                    obj.parameterLabels = {'OffsetL-R','ScaleL-R','QL-R'};
                    obj.parameterBounds = [-inf -inf -inf;
                        +inf +inf +inf];
                    obj.Zinput = @(D)(behav.GLM_QlearningInputs(D));
                    N = 0.4;
                    obj.ZL = @(P,in)( P(1) + P(2).*in(:,1) + P(3).*in(:,2) );
                    obj.parameterStart = @()([10*rand(1,3)]);
                    
                case 'C^N-subset-Qlearning-noBias-2AFC'
                    obj.parameterLabels = {'ScaleL-R','QL-R'};
                    obj.parameterBounds = [-inf -inf;
                        +inf +inf];
                    obj.Zinput = @(D)(behav.GLM_QlearningInputs(D));
                    N = 0.4;
                    obj.ZL = @(P,in)( P(1).*in(:,1) + P(2).*in(:,2) );
                    obj.parameterStart = @()([10*rand(1,2)]);
                    
                case 'C50-subset-centred'
                    obj.parameterLabels = {'Offset_L','ScaleL_L','Offset_R','ScaleR_R','N','C50'};
                    obj.parameterBounds = [-inf -inf -inf -inf 0 0.001;
                        +inf +inf +inf +inf 3 0.8];
                    
                    obj.Zinput = @(D)([D.stimulus(:,1) D.stimulus(:,2)]);
                    obj.ZL = @(P,in)( P(1) + P(2).*(in(:,1).^P(5))./(in(:,1).^P(5) + P(6)^P(5)) - P(2)*mean(in(:,1).^P(5)./(in(:,1).^P(5) + P(6)^P(5))) );
                    obj.ZR = @(P,in)( P(3) + P(4).*(in(:,2).^P(5))./(in(:,2).^P(5) + P(6)^P(5)) - P(4)*mean(in(:,2).^P(5)./(in(:,2).^P(5) + P(6)^P(5))) );
                otherwise
                    error('Model does not exist');
                    
            end
            
            if isempty(obj.parameterStart)
                obj.parameterStart = zeros(1,length(obj.parameterLabels));
            end
            
        end
        
        function obj = addLapse(obj)
            if isempty(obj.ZL)
                error('Set model first');
            end
            
            obj.parameterLabels = [obj.parameterLabels 'LapseRate'];
            obj.parameterBounds = [obj.parameterBounds, [0;1]];
            obj.parameterStart = [obj.parameterStart 0];
            
            obj.lapseFlag=obj.ContrastDimensions;
        end
        
        function obj = removeRepeats(obj,varargin)
            if nargin > 1
                threshold = varargin{1};
            else
                threshold = 1;
            end
            
            obj.data = obj.getrow(obj.data,obj.data.repeatNum<=threshold);
        end
        
        function obj = fit(obj)
            %Non crossvalidated fitting
            
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            %Trim first 5 trials
            %             obj.data = obj.getrow(obj.data,6:length(obj.data.response));
            
            %Remove trials with repeats
            obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            %             options = optimoptions('fmincon','algorithm','interior-point','MaxFunEvals',100000,'MaxIter',10000);
            options = optimset('algorithm','interior-point','MaxFunEvals',100000,'MaxIter',10000);
            
            responses = obj.data.response;
            
            if isempty(obj.regularise)
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), responses));
            else
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(obj.data), responses) + obj.regularise(b));
            end
            
            %             if exist('opti','file')==2 %use opti toolbox if loaded
            %                 options = optiset('display','final','solver','NOMAD');
            %                 Opt = opti('fun',objective,'bounds',obj.parameterBounds(1,:),obj.parameterBounds(2,:),'x0',obj.parameterStart(),'options',options);
            %                 [p,~,exitflag,~] = Opt.solve;
            %                 obj.parameterFits = p';
            %                 if exitflag < 0
            %                     obj.parameterFits = nan(1,length(obj.parameterLabels));
            %                 end
            %             else
            [obj.parameterFits,~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
            if ~any(exitflag == [1,2])
                obj.parameterFits = nan(1,length(obj.parameterLabels));
            end
            %             end
            
        end
        
        function r2(obj)
            if isempty(obj.parameterFits)
                error('Fit first');
            end
            
            %non-cv log likelihood for fitted model
            ll_model = obj.calculateLogLik(obj.parameterFits,obj.data.stimulus,obj.data.response)/length(obj.data.response);
            
            %another thing to do is just take the mode of the pred
            %probabilities and see if they match the data
            phats = obj.calculatePhat(obj.parameterFits,obj.data.stimulus);
            classify = nan(length(obj.data.response),1);
            for t = 1:length(obj.data.response)
                [~,rhat] = max(phats(t,:));
                
                if rhat == obj.data.response(t)
                    classify(t)=1;
                else
                    classify(t)=0;
                end
            end
            
            %             keyboard;
            
            %             %non-cv log likelihood for naive intercept model (no contrast)
            %             obj = obj.setModel('Offset');
            %             obj=obj.fit;
            %             ll_naive = obj.calculateLogLik(obj.parameterFits,obj.data.stimulus,obj.data.response);
            
            r2 = 1 + ll_model/obj.guess_bpt;
            disp(['McFadden Pseudo-R^2: ' num2str(r2)]);
            disp(['Proportion correctly classified: ' num2str(mean(classify))]);
            
            
        end
        
        function [obj,varargout] = fitCV(obj,varargin)
            %Crossvalidated fitting
            
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            %Trim first 5 trials
            %             obj.data = obj.getrow(obj.data,6:length(obj.data.response));
            
            %Remove trials with repeats
            %             obj.data = obj.getrow(obj.data,obj.data.repeatNum==1);
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',2000,'Display','off');
            
            if isempty(varargin)
                C = cvpartition(length(obj.data.response),'LeaveOut');
            else
                C = cvpartition(obj.data.response,'KFold',varargin{1});
            end
            
            obj.parameterFits = nan(C.NumTestSets,length(obj.parameterLabels));
            obj.p_hat = nan(length(obj.data.response),3);
            for f=1:C.NumTestSets
                %                 disp(['Model: ' obj.modelString '. Fold: ' num2str(f) '/' num2str(C.NumTestSets)]);
                trainIdx = find(C.training(f)==1);
                testIdx = find(C.test(f)==1);
                
                inputs = obj.Zinput(obj.data);
                trainInputs = inputs(trainIdx,:);
                testInputs = inputs(testIdx,:);
                
                trainResponses = obj.data.response(trainIdx);
                testResponse = obj.data.response(testIdx);
                
                if isempty(obj.regularise)
                    objective = @(b) ( obj.calculateLogLik(b, trainInputs, trainResponses) );
                else
                    objective = @(b) ( obj.calculateLogLik(b, trainInputs, trainResponses) + obj.regularise(b));
                end
                
                %                 if exist('opti','file')==2 %use opti toolbox if loaded
                %                     options = optiset('solver','NLOPT');
                %                     Opt = opti('fun',objective,'bounds',obj.parameterBounds(1,:),obj.parameterBounds(2,:),'x0',obj.parameterStart(),'options',options);
                %                     [p,~,exitflag,~] = Opt.solve;
                %                     obj.parameterFits(f,:) = p';
                %                     if exitflag < 0
                %                         obj.parameterFits = nan(1,length(obj.parameterLabels));
                %                     end
                %                 else
                [obj.parameterFits(f,:),~,exitflag] = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
                if ~any(exitflag == [1,2])
                    obj.parameterFits(f,:) = nan(1,length(obj.parameterLabels));
                end
                %                 end
                
                
                phat = obj.calculatePhat(obj.parameterFits(f,:), testInputs);
                obj.p_hat(testIdx,:) = phat;
                %
                %                 for i = 1:length(testResponse)
                %                     obj.p_hat(testIdx(i),1) = phat(i,testResponse(i));
                %                 end
                
                %                 LL(f)=mean(-log2(obj.p_hat(testIdx)));
            end
            
            %             varargout = {LL};
        end
        
        function h = plotData(obj)
            
            switch(obj.ContrastDimensions)
                case 1
                    contrast1D = obj.data.stimulus(:,2) - obj.data.stimulus(:,1);
                    uniqueC1D = unique(contrast1D);
                    prop=[];
                    prop_ci=[];
                    for c = 1:length(uniqueC1D)
                        D = obj.getrow(obj.data,contrast1D == uniqueC1D(c));
                        p = sum([D.response==1 D.response==2 D.response==3],1)/length(D.response);
                        
                        [p,pci]=binofit(sum([D.response==1 D.response==2 D.response==3],1),length(D.response),0.05);
                        
                        prop_ci(:,:,c) = pci;
                        
                        prop = [prop;p];
                        %                         bse = sqrt((p.*(1-p)/length(D.response)));
                        %                         binom_se = [binom_se;bse];
                    end
                    
                    
                    
                    %                     plot(uniqueC1D,prop,'.','MarkerSize',20);
                    
                    err = [prop - squeeze(prop_ci(:,1,:))', squeeze(prop_ci(:,2,:))' - prop];
                    
                    if max(obj.data.response) == 2 %for 2AFC tasks
                        r = 2;
                    else
                        r = 3;
                    end
                    
                    if length(uniqueC1D)<50
                        hx=errorbar(repmat(uniqueC1D,1,r),prop(:,1:r),err(:,[1:r]),err(:,[4:(3+r)]),'.','MarkerSize',20);
                    else
                        plot(repmat(uniqueC1D,1,r),prop(:,1:r),'o');
                    end
                    xlabel('contrast');
                    ylabel('P( choice | contrast)');
                    
                    set(gca,'box','off');
                    h=gca;
                    
                case 2
                    uniqueCL = unique(obj.data.stimulus(:,1));
                    uniqueCR = unique(obj.data.stimulus(:,2));
                    prop=nan(length(uniqueCL),length(uniqueCR),3);
                    
                    for cl = 1:length(uniqueCL)
                        for cr = 1:length(uniqueCR)
                            E = obj.getrow(obj.data,obj.data.stimulus(:,1) == uniqueCL(cl) & obj.data.stimulus(:,2) == uniqueCR(cr));
                            for i=1:3
                                prop(cl,cr,i) = sum(E.response==i)/length(E.response);
                            end
                        end
                    end
                    
                    titles = {'p( Left | c)','p( Right | c)','p( NoGo | c)'};
                    for i=1:3
                        h(i)=subplot(2,3,i);
                        %                         p_plot = prop(:,:,i);
                        %                         p_plot = [p_plot; nan(1,length(uniqueCR))];
                        %                         p_plot = [p_plot, nan(length(uniqueCL)+1,1)];
                        %                         pcolor([uniqueCR; 1],[uniqueCL; 1],p_plot); caxis([0 1]); shading('flat');
                        %                         set(gca,'box','off');
                        %                         set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1);
                        imagesc(prop(:,:,i)); caxis([0 1]);
                        set(gca,'xtick',1:length(uniqueCR),'ytick',1:length(uniqueCL),'xticklabels',uniqueCR,'yticklabels',uniqueCL);
                        set(gca,'YDir','normal','box','off');
                        %                         xlabel('Contrast right');
                        %                         ylabel('Contrast left');
                        title(titles{i});
                        axis square;
                        %                         set(gca,'XTick','','YTick',0:0.1:0.5);
                        if i > 1
                            set(gca,'XTick',[],'ytick',[]);
                        end
                        
                        if i == 1
                            %                                 xlabel('Contrast right');
                            ylabel('Contrast left');
                        end
                    end
            end
            
            set(gcf,'color','w');
            
        end
        
        function plotRT(obj,varargin)
            if cell2mat(strfind(varargin,'log')) == 1
                logFlag=1;
                obj.data.RT = log(obj.data.RT);
                %                 fcn = @(rt)(log(rt));
            else
                logFlag=0;
            end
            
            figure('color','w');
            resplabel = {'Left choices','Right choices'};
            
            cVal = unique(obj.data.stimulus(:));
            i=1;
            meanRT = nan(length(cVal),length(cVal),2);
            medianRT = nan(length(cVal),length(cVal),2);
            for cl = 1:length(cVal)
                for cr = 1:length(cVal)
                    subplot(length(cVal),length(cVal),i);
                    
                    try
                        idx = obj.data.stimulus(:,1) == cVal(cl) & obj.data.stimulus(:,2) == cVal(cr) & obj.data.response==1;
                        distributionPlot(obj.data.RT(idx),'histOpt',0,'widthDiv',[2 1],'histOri','left','color','b','showMM',0,'globalNorm',2)
                        meanRT(cl,cr,1) = mean(obj.data.RT(idx));
                        medianRT(cl,cr,1) = median(obj.data.RT(idx));
                    catch
                    end
                    
                    try
                        idx = obj.data.stimulus(:,1) == cVal(cl) & obj.data.stimulus(:,2) == cVal(cr) & obj.data.response==2;
                        distributionPlot(obj.data.RT(idx),'histOpt',0,'widthDiv',[2 2],'histOri','right','color','r','showMM',0,'globalNorm',2)
                        set(gca,'xtick','','xcolor','w','box','off','ytick','','ycolor','w');
                        meanRT(cl,cr,2) = mean(obj.data.RT(idx));
                        medianRT(cl,cr,2) = median(obj.data.RT(idx));
                        %                     ylim([0 1.5]);
                    catch
                    end
                    
                    if cl == length(cVal)
                        set(gca,'xcolor','k');
                        xlabel(['CR=' num2str(cVal(cr))]);
                    end
                    
                    if cr == 1
                        set(gca,'ycolor','k','ytick',0:0.4:1.5);
                        ylabel(['CL=' num2str(cVal(cl))]);
                    end
                    
                    %                     if cr == length(cVal)
                    %                         set(gca,'ycolor','k');
                    %                         xlabel(['CR=' num2str(cVal(cr))]);
                    %                     end
                    %
                    i=i+1;
                end
            end
            
            figure('color','w');
            subplot(2,2,1);
            imagesc(meanRT(:,:,1)); xlabel('CR'); ylabel('CL'); axis square; title('Left choices');
            set(gca,'xtick',1:length(cVal),'xticklabels',cVal,'ytick',1:length(cVal),'yticklabels',cVal);
            set(gca,'ydir','normal'); caxis([min(meanRT(:)) max(meanRT(:))]);
            subplot(2,2,2);
            imagesc(meanRT(:,:,2)); xlabel('CR'); ylabel('CL'); axis square; title('Right choices');
            set(gca,'xtick',1:length(cVal),'xticklabels',cVal,'ytick',1:length(cVal),'yticklabels',cVal);
            set(gca,'ydir','normal'); caxis([min(meanRT(:)) max(meanRT(:))]);
            subplot(2,2,3);
            imagesc(medianRT(:,:,1)); xlabel('CR'); ylabel('CL'); axis square; title('Left choices');
            set(gca,'xtick',1:length(cVal),'xticklabels',cVal,'ytick',1:length(cVal),'yticklabels',cVal);
            set(gca,'ydir','normal'); caxis([min(meanRT(:)) max(meanRT(:))]);
            subplot(2,2,4);
            imagesc(medianRT(:,:,2)); xlabel('CR'); ylabel('CL'); axis square; title('Right choices');
            set(gca,'xtick',1:length(cVal),'xticklabels',cVal,'ytick',1:length(cVal),'yticklabels',cVal);
            set(gca,'ydir','normal'); caxis([min(meanRT(:)) max(meanRT(:))]);
            colormap(flipud(hot));
            
            
            %Plot histograms at the 4 extreme contrast settings
            figure('color','w');
            cVal = cVal([1,end]);
            labels = {'C=0','High CR','High CL','High CL&CR'};
            choiceLab = {'Chose L','Chose R'};
            for r = 1:2
                subplot(1,2,r); hold on;
                for cl = 1:length(cVal)
                    for cr = 1:length(cVal)
                        idx = obj.data.stimulus(:,1) == cVal(cl) & obj.data.stimulus(:,2) == cVal(cr) & obj.data.response==r;
                        histogram(obj.data.RT(idx),'DisplayStyle','stairs');
                    end
                end
                xlabel(['RT ' choiceLab{r}]);
                legend(labels);
            end
            %             for r = 1:2
            %                 subplot(1,2,r);
            %                 hist(obj.data.RT(obj.data.response==r));
            %                 xlim([0 1.5]);
            %                 title(resplabel{r});
            %                 set(gca,'box','off');
            %             end
            
        end
        
        function fig = plotFit(obj)
            if size(obj.parameterFits,1)>0
                
                fig = figure('color','w');
                ha = tight_subplot(4,4,0.05,[0.05 0.05],[0.05 0.01]);
                set(ha([1 5 9 13]),'dataaspectratio',[1 1 1],'xlim',[0 1],'ylim',[0 1]);
                set(ha,'xcolor','none','ycolor','none');
%                 set(ha([3 4 7 8 11 12 15 16]),'ytick','','ycolor','none');
%                 set(ha(14:16),'xticklabelmode','auto');
%                 set(ha([2 6 10 14]),'yticklabelmode','auto');
                
                cont = obj.data.stimulus;
                resp = obj.data.response;
%                 cVals = unique(cont(:));
                cVals = unique(min(cont,[],2));
                
                if length(cVals)>4
                    %Take top 4 instead
                    tab=tabulate(min(cont,[],2));
                    tab=sortrows(tab,3,'descend');
                    cVals = tab(1:4,1);
                    cVals = sort(cVals);
%                     keyboard;
                end
                
                numPedestals = length(cVals);
                cols = [0 0.4470 0.7410;
                    0.8500 0.3250 0.0980;
                    0.4940    0.1840    0.5560];
                
                set(ha,'xlim',[-1 1]*(max(cont(:))+0.1),'ylim',[0 1]);
               
                for ped = 1:numPedestals

                    %Add cartoon on the left of the contrast conditions
                    %plotted
                    ha_cartoon = ha( 4*(ped-1) + 1);
                    ped_idx = min(cont,[],2)==cVals(ped);
                    ped_c = unique(cont(ped_idx,:),'rows');
                    
                    plot(ha_cartoon, ped_c(:,2), ped_c(:,1), 'k.','markersize',15);
                    set(ha_cartoon,'dataaspectratio',[1 1 1],'xlim',[0 max(cont(:))+0.1],'ylim',[0 max(cont(:))+0.1],'box','off');
                    set(ha_cartoon,'ytick',unique(ped_c(:)),'xtick',unique(ped_c(:)));
                    
                    ha_ped = [ ha(4*(ped-1) + 2) ha(4*(ped-1) + 4) ha(4*(ped-1) + 3)];
                    set(ha_ped(1),'yticklabelmode','auto','ycolor','k');
                    if ped == numPedestals
                        set(ha_ped,'xticklabelmode','auto','xcolor','k');
                    end
                    
                    hold(ha_ped(1),'on');
                    hold(ha_ped(2),'on');
                    hold(ha_ped(3),'on');
                    
                    set(ha_ped,'colororder',cols);
                    
                    %Plot actual datapoints
                    ped_idx = min(cont,[],2)==cVals(ped);
                    ped_c_diff = diff(cont(ped_idx,:),[],2);
                    ped_r = resp(ped_idx);
                    uC = unique(ped_c_diff);
                    ph=[];
                    for c = 1:length(uC)
                        r = ped_r(ped_c_diff==uC(c));
                        [ph(c,:),pci] = binofit(sum([r==1 r==2 r==3],1),length(r));
                        for ch=1:3
                            l(1)=line(ha_ped(ch),[1 1]*uC(c),pci(ch,:));
                            %                                         l(2)=line(ha(ped),[uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,1));
                            %                                         l(3)=line(ha(ped),[uC(c)-0.03 uC(c)+0.03],[1 1]*pci(ch,2));
                            set(l,'Color',cols(ch,:),'Linewidth',0.5);
                            %                             l.Color = cols{ch};
                            %                             l.LineWidth=1;
                        end
                    end
                    set(ha(ped),'ColorOrderIndex',1);
                    
                    for ch=1:3
                        plot(ha_ped(ch),uC,ph(:,ch),'.','markersize',15,'color',cols(ch,:));
                    end
                    
                    %Plot predictions
                    testCont = [linspace(max(abs(uC))+0.1,0,100)' zeros(100,1); zeros(100,1) linspace(0,max(abs(uC))+0.1,100)'] + cVals(ped);
                    p_hat = obj.calculatePhat(obj.parameterFits,testCont);
                    set(gca,'ColorOrderIndex',1);
                    for ch=1:3
                        plot(ha_ped(ch),diff(testCont,[],2),p_hat(:,ch),'linewidth',1,'color',cols(ch,:));
                    end
                    
                    
                    %                                 if ped == 1
                    %                                     legend(ha_ped(3),{'L','R','NG'}); legend boxoff;
                    %                                 end
                    
                    
                    %                                 ylabel(['ped: ' num2str(cVals(ped))]);
                    
                    if ped==1
                        title(ha_ped(1),'pLeft');
                        title(ha_ped(2),'pRight');
                        title(ha_ped(3),'pNoGo');
                    end
                    
                    
                    
                end
                
                set(get(fig,'children'),'fontsize',6);
                
%                 keyboard;
            else
                error('Model not fitted (non-crossvalidated) yet');
            end
        end
        
        function plotPedestal(obj)
            if isempty(obj.parameterFits)
                error('Need to fit model first');
            end
            
            cVals = unique(obj.data.stimulus(:));
            prop=nan(length(cVals),length(cVals),3);
            
            for cl = 1:length(cVals)
                for cr = 1:length(cVals)
                    E = obj.getrow(obj.data,obj.data.stimulus(:,1) == cVals(cl) & obj.data.stimulus(:,2) == cVals(cr));
                    for i=1:3
                        prop(cl,cr,i) = sum(E.response==i)/length(E.response);
                    end
                    pd.propChooseLeft(cl,cr) = prop(cl,cr,1);
                    pd.propChooseRight(cl,cr) = prop(cl,cr,2);
                    pd.propChooseNogo(cl,cr) = prop(cl,cr,3);
                end
            end
            
            cTestVals = min(cVals):0.01:2;
            selectContrast{1} = [reshape(reshape(repmat(cVals, 1, length(cTestVals)), length(cVals), length(cTestVals))', 1, length(cVals)*length(cTestVals)); ...
                repmat(cTestVals, 1, length(cVals))];
            selectContrast{2} = selectContrast{1}([2 1], :);
            
            predictionsSelect{1} = obj.calculatePhat(obj.parameterFits, selectContrast{1}')';
            predictionsSelect{2} = obj.calculatePhat(obj.parameterFits, selectContrast{2}')';
            f3 = figure; %set(f3, 'Position', [  -1896         507        1058         405]);
            
            colors = [        0    0.4470    0.7410;
                0.8500    0.3250    0.0980;
                0.9290    0.6940    0.1250];
            
            for ped = 1:length(cVals)-1
                
                subplot(1, length(cVals)-1, ped)
                
                for r = 1:3
                    plot(-(cVals(ped:end)-cVals(ped)), prop(ped:end,ped,r), 'o','Color',colors(r,:));
                    
                    hold on;
                    plot((cVals(ped:end)-cVals(ped)), prop(ped,ped:end,r), 'o','Color',colors(r,:));
                    
                    
                    
                    plot(-(cTestVals(cTestVals>=cVals(ped))-cVals(ped)), predictionsSelect{2}(r,selectContrast{1}(1,:)==cVals(ped)&selectContrast{1}(2,:)>=cVals(ped)),'Color',colors(r,:));
                    
                    plot((cTestVals(cTestVals>=cVals(ped))-cVals(ped)), predictionsSelect{1}(r,selectContrast{2}(2,:)==cVals(ped)&selectContrast{2}(1,:)>=cVals(ped)),'Color',colors(r,:));
                    
                end
                
                title(['pedestal = ' num2str(cVals(ped))]);
                xlabel('delta C');
                ylim([0 1]);
                xlim([-1 1]);
                set(gca,'box','off');
                %                 makepretty
            end
            set(gcf,'color','w');
        end
        
        function plotPredVsActual(obj, varargin)
            
            plotParams.Color = 'r';
            plotParams.ax = gca;
            if ~isempty(varargin)
                plotParams = mergeStructs(varargin{1},plotParams);
            end
            
            switch(obj.ContrastDimensions)
                case 1
                    contrast1D = diff(obj.data.stimulus, [], 2);
                    uniqueC1D = unique(contrast1D);
                    nC = length(uniqueC1D);
                    prop=zeros(nC,3);
                    prop_ci=zeros(nC,3,2);
                    for c = 1:length(uniqueC1D)
                        D = obj.getrow(obj.data,contrast1D == uniqueC1D(c));
                        respSum = sum([D.response==1 D.response==2 D.response==3],1);
                        p = respSum/length(D.response);
                        
                        [p,pci]=binofit(respSum,length(D.response),0.05);
                        
                        prop_ci(c,:,:) = pci;
                        
                        prop(c,:) = p;
                    end
                    
                    if max(obj.data.response) == 2 %for 2AFC tasks
                        rMax = 2;
                    else
                        rMax = 3;
                    end
                    
                    evalCon = unique(obj.data.stimulus,'rows');
                    evalC1d = evalCon(:,2) - evalCon(:,1);
                    [~,sortIdx]=sort(evalC1d);
                    evalCon = evalCon(sortIdx,:);
                    phat = obj.calculatePhat(obj.parameterFits,evalCon);
                    
                    rSymbols = {'o', '.', 'x'};
                    
                    for c = 1:length(uniqueC1D)
                        for r = 1:rMax
                            plot(plotParams.ax, phat(c,r), prop(c,r), rSymbols{r}, 'Color', plotParams.Color)
                            hold on;
                        end
                        for r = 1:rMax
                            plot(plotParams.ax, phat(c,r)*ones(1,2), squeeze(prop_ci(c,r,:)), 'Color', plotParams.Color)
                        end
                    end
                    plot(plotParams.ax, [0 1], [0 1], 'k--');
                    
                    ylabel('actual probability');
                    xlabel('predicted probability');
                    legend({'left resp', 'right resp', 'nogo'}, 'Location', 'Best');
                    
                    axis square
                    box off
                    
                case 2
                    fprintf(1, 'plotPredVsActual not yet implemented for 2D task\n')
            end
        end
        
        function h = plotParams(obj)
            if size(obj.parameterFits,1)==1
                bar(obj.parameterFits);
                set(gca,'XTickLabel',obj.parameterLabels,'XTick',1:numel(obj.parameterLabels));
                title(obj.modelString);
                h=gca;
            end
        end
        
        function phat = calculatePhat(obj,testParams,inputs)
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            switch(obj.lapseFlag)
                case 0
                    lapse = 0;
                case 1
                    lapse = testParams(end);
                case 2
                    lapse = testParams(end);
            end
            
            if isempty(obj.ZR) %if a AFC task then no ZR is defined, only pL vs pR
                zl = obj.ZL(testParams,inputs);
                pL = lapse + (1-2*lapse)*exp(zl)./(1+exp(zl));
                pR = 1 - pL;
                N = length(pL);
                phat = [pL pR zeros(N,1)];
            else
                zl = obj.ZL(testParams,inputs);
                zr = obj.ZR(testParams,inputs);
                pL = (1-lapse)*exp(zl)./(1+exp(zl)+exp(zr));
                pR = (1-lapse)*exp(zr)./(1+exp(zl)+exp(zr));
                pNG = 1 - (pL + pR);
                
                phat = [pL pR pNG];
            end
            
            %             if any(phat<0)
            %                 keyboard;
            %             end
        end
        
        function obj = addData(obj,type)
            %Adds extra data depending on what is required
            
            if strcmp(obj.expRef,'none')
                error('not coded for custom struct inputs');
            end
            
            switch(type)
                case 'lick'
                    block = dat.loadBlock(obj.expRef);
                    trials = block.trial;
                    
                    disp('Loading lick data...');
                    try
                        L=load(dat.expFilePath(obj.expRef, 'Timeline', 'm'));
                        tseries = L.Timeline.rawDAQTimestamps;
                        lickseries = L.Timeline.rawDAQData(:,7);
                        
                        for t=1:block.numCompletedTrials
                            start = trials(t).trialStartedTime;
                            finish = trials(t).trialEndedTime;
                            idx = (start < tseries) & (tseries < finish);
                            
                            obj.data.lickenergy(t,1) = sum(lickseries(idx).^2);
                        end
                    catch
                        warning('No lick data found.. Setting to NaNs');
                        obj.data.lickenergy = nan(block.numCompletedTrials,1);
                    end
                    
                case 'pupil'
                    block = dat.loadBlock(obj.expRef);
                    trials = block.trial;
                    
                    disp('Loading eye movie...');
                    
                    try
                        vidFilename = ['\\zserver\Data\EyeCamera\' obj.expRef(14:end) '\' obj.expRef(1:10) '\' obj.expRef(12) '\eye.mj2'];
                        v = VideoReader(vidFilename);
                        
                        for t=1:block.numCompletedTrials
                            %                             disp(t);
                            start = trials(t).trialStartedTime;
                            finish = trials(t).trialEndedTime;
                            
                            v.CurrentTime = start; a=0; pixCounts=0;
                            while v.CurrentTime < finish
                                frame = 255 - (v.readFrame*5.5);
                                frame(frame>0)=1;
                                frame = frame(50:110,80:180);
                                pixCounts = pixCounts + sqrt(sum(frame(:)));
                                a=a+1;
                            end
                            
                            obj.data.pupilenergy(t,1) = pixCounts/a;
                        end
                    catch
                        warning('No eye data found.. Setting to NaNs');
                        obj.data.pupilenergy = nan(block.numCompletedTrials,1);
                        
                    end
            end
            
            disp('Done!');
            
        end
        
        function cov=paramCov(obj,varargin)
            if isempty(obj.parameterFits)
                error('Fit model first!');
            end
            
            c = obj.data.stimulus;
            H = obj.hessian(obj.parameterFits,c);
            
            %Calculate fisher info
            F = -sum(H,3);
            cov = inv(F);
            
            figure;
            ax1=subplot(1,2,1);
            imagesc(cov); title('Covariance'); axis square;
            ax2=subplot(1,2,2);
            imagesc(corrcov(cov)); caxis([-1 1]); title('Correlation'); axis square;
            set(get(gcf,'children'),'XTickLabel',{'B_L','B_R','S_L','S_R'},'xtick',1:4,'yTickLabel',{'B_L','B_R','S_L','S_R'},'ytick',1:4);
            
            cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
                linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
            colormap(ax2,cmap); colorbar;
        end
        
        function paramCovSimulate(obj)
            %Plots covariance/variance of parameters as a function of
            %varying stimulus and behavioural configurations. This is to
            %explore other ways of setting up the task such that the
            %parameter estimation can be less variable/less covariable.
            if isempty(obj.parameterFits)
                error('Fit first');
            end
            
            %1) Same parameters but vary the distribution of contrasts
            if obj.ContrastDimensions == 1
                n = obj.parameterFits(end-1);
                c50 = obj.parameterFits(end);
                cfn = @(c)(c.^n)./(c.^n + c50.^n);
                numTrials = length(obj.data.response);
                
                widths = linspace(0.01,0.5,50);
                
                meanVarS = []; cv = [];
                for w = 1:length(widths)
                    varS = [];
                    for iter = 1:50
                        crnd = widths(w)*randn(numTrials,1);
                        c = zeros(length(crnd),2);
                        c(crnd<0,1) = -crnd(crnd<0);
                        c(crnd>0,2) = crnd(crnd>0);
                        
                        H = obj.hessian(obj.parameterFits,c);
                        
                        F = -sum(H,3);
                        cov = inv(F);
                        varS(iter,:) = [cov(1,1) cov(2,2) cov(3,3) cov(4,4) cov(3,1) cov(2,1) cov(4,3)];
                        
                    end
                    meanVarS(w,:) = mean(varS);
                    
                    corrSLBL(w,1) = meanVarS(w,5)/sqrt(meanVarS(w,1)*meanVarS(w,3));
                    
                    cv(w,1) = sqrt(meanVarS(w,1))/obj.parameterFits(1);
                    cv(w,2) = sqrt(meanVarS(w,2))/obj.parameterFits(3);
                    cv(w,3) = sqrt(meanVarS(w,3))/obj.parameterFits(2);
                    cv(w,4) = sqrt(meanVarS(w,4))/obj.parameterFits(4);
                end
                figure;
                subplot(3,1,1);
                plot(widths,[meanVarS corrSLBL],'LineWidth',2); xlabel('Standard Deviation of Contrasts');
                legend({'Var[B_L]','Var[B_R]','Var[S_L]','Var[S_R]','Cov[B_L,S_L]','Cov[B_L,B_R]','Cov[S_L,S_R]','Corr[B_L,S_L]'});
                
                hold on;
                this = std(diff(obj.data.stimulus,[],2));
                line([this,this],[-1 1]);
            end
            
            %2) Same contrasts but vary GovsNoGo bias
            fitted_p = obj.parameterFits;
            Biases = linspace(-5,5,50);
            varS = []; CVs = [];
            for bias = 1:length(Biases)
                new_p = fitted_p;
                new_p([1 3]) = Biases(bias);
                H = obj.hessian(new_p,obj.data.stimulus);
                F = -sum(H,3);
                cov = inv(F);
                varS(bias,:) = [cov(1,1) cov(2,2) cov(3,3) cov(4,4) cov(3,1) cov(2,1) cov(4,3)];
                CVs(bias,1) = sqrt(cov(1,1))/new_p(1);
                CVs(bias,2) = sqrt(cov(2,2))/new_p(3);
                CVs(bias,3) = sqrt(cov(3,3))/new_p(2);
                CVs(bias,4) = sqrt(cov(4,4))/new_p(4);
            end
            corrS = varS(:,5)./sqrt(varS(:,1).*varS(:,3));
            subplot(3,1,2);
            pGO = 1- 1./(1+2*exp(Biases)); xlim([0 1]);
            plot(pGO,[varS,corrS],'LineWidth',2); xlabel('pGO');
            %             legend({'cv[B_L]','cv[B_R]','cv[S_L]','cv[S_R]'});
            
            this = mean(fitted_p([1 3]));
            this = 1- 1./(1+2*exp(this));
            line([this,this],[-1 1]);
            
            %3) Same contrast but vary LvNG bias
            fitted_p = obj.parameterFits;
            Biases = linspace(-5,5,50);
            varS = [];
            for bias = 1:length(Biases)
                new_p = fitted_p;
                new_p(1) = Biases(bias);
                H = obj.hessian(new_p,obj.data.stimulus);
                F = -sum(H,3);
                cov = inv(F);
                varS(bias,:) = [cov(1,1) cov(2,2) cov(3,3) cov(4,4) cov(3,1) cov(2,1) cov(4,3)];
            end
            corrS = varS(:,5)./sqrt(varS(:,1).*varS(:,3));
            subplot(3,1,3);
            plot(Biases,[varS corrS],'LineWidth',2); xlabel(['bL . bR=' num2str(fitted_p(3))]);
            this = fitted_p(1);
            %             this = 1- 1./(1+2*exp(this));
            line([this,this],[-1 1]);
        end
        
        function plotPspace_mnrvsnested(obj)
            %use mnrfit to fit rudimentary models MNR and NESTED and
            %compare directly
            cont = obj.data.stimulus;
            resp = obj.data.response;
            resp_nes = resp; resp_nes(resp_nes==3)=0; resp_nes=resp_nes+1;
            
            B_mnr = mnrfit(cont,resp);
            B_nes = mnrfit(cont,resp_nes,'model','hierarchical');
            
            cVals = linspace(0,0.54,1000);
            [cr,cl]=meshgrid(cVals);
            cont = [cl(:) cr(:)];
            
            P_mnr = mnrval(B_mnr,cont);
            P_nes = mnrval(B_nes,cont,'model','hierarchical');
            P_nes = [P_nes(:,2:3) P_nes(:,1)];
            
            %Plot
            figure('color','w');
            labels = {'pL','pR','pNG'};
            %             cVals=log(cVals);
            for r = 1:3
                subplot(1,3,r);
                [~,ax]=contour(cVals,cVals,reshape(P_mnr(:,r),size(cl)));
                ax.LineWidth=1; hold on;
                [~,ax]=contour(cVals,cVals,reshape(P_nes(:,r),size(cl)));
                ax.LineWidth=2;
                ax.LineStyle=':';
                
                set(gca,'ydir','normal','box','off');
                xlabel('CR'); ylabel('CL'); title(labels{r}); axis square; caxis([0 1]);
                %
                % %                 subplot(2,3,r+3);
                %                 imagesc(cVals,cVals,reshape(P_nes(:,r),size(cl)));
                %                 set(gca,'ydir','normal');
                %                 xlabel('CR'); ylabel('CL'); title(labels{r});axis square;  caxis([0 1]);
            end
            %             keyboard;
        end
        
        function bootstrapFit(obj)
            if isempty(obj.ZL)
                error('Please set a model first using method setModel(...)');
            end
            
            numIter = 300;
            
            options = optimoptions('fmincon','UseParallel',0,'MaxFunEvals',100000,'MaxIter',10000);
            bootstrap_params = nan(numIter,length(obj.parameterLabels));
            
            T = struct2table(obj.data);
            for iter = 1:numIter
                datTemp = datasample(T,height(T));
                objective = @(b) (obj.calculateLogLik(b, obj.Zinput(datTemp), datTemp.response));
                bootstrap_params(iter,:) = fmincon(objective, obj.parameterStart(), [], [], [], [], obj.parameterBounds(1,:), obj.parameterBounds(2,:), [], options);
            end
            
            gplotmatrix(bootstrap_params,[],[],[],[],[],[],[],obj.parameterLabels);
        end
        
    end
    
    
    
    methods (Access= {?GLM})
        function logLik = calculateLogLik(obj,testParams, inputs, responses)
            phat = obj.calculatePhat(testParams, inputs);
            logLik = -sum(log2( phat(sub2ind(size(phat), [1:length(responses)]', responses)) ));
        end
        function row = getrow(~,D,numrow)
            % Version 1.0 9/18/03
            % by Joern Diedrichsen
            % http://www.icn.ucl.ac.uk/motorcontrol/toolboxes/toolbox_util.htm
            
            if (~isstruct(D))
                error('D must be a struct');
            end;
            
            field = fieldnames(D);
            row=[];
            for f=1:length(field)
                F = getfield(D,field{f});
                if iscell(F)
                    row = setfield(row,field{f},F(numrow,:));
                else
                    row = setfield(row,field{f},F(numrow,:));
                end
            end
            
        end
        
        function H = hessian(obj,params,contrasts)
            
            %Use untransformed contrasts in ZL and ZR eqns
            zl=obj.ZL(params,contrasts);
            zr=obj.ZR(params,contrasts);
            pL = exp(zl)./(1+exp(zl)+exp(zr));
            pR = exp(zr)./(1+exp(zl)+exp(zr));
            
            n = obj.parameterFits(end-1);
            c50 = obj.parameterFits(end);
            cfn = @(c)(c.^n)./(c.^n + c50.^n);
            CL = cfn(contrasts(:,1));
            CR = cfn(contrasts(:,2));
            
            H = zeros(4,4,length(contrasts));
            H(1,1,:) = pL.*(pL-1);
            H(2,2,:) = pR.*(pR-1);
            H(3,3,:) = (CL.^2).*pL.*(pL-1);
            H(4,4,:) = (CR.^2).*pR.*(pR-1);
            H(2,1,:) = pL.*pR;
            H(1,2,:) = pL.*pR;
            H(3,1,:) = CL.*pL.*(pL-1);
            H(1,3,:) = CL.*pL.*(pL-1);
            H(4,1,:) = CR.*pL.*pR;
            H(1,4,:) = CR.*pL.*pR;
            H(3,2,:) = CL.*pL.*pR;
            H(2,3,:) = CL.*pL.*pR;
            H(4,2,:) = CR.*pR.*(pR-1);
            H(2,4,:) = CR.*pR.*(pR-1);
            H(4,3,:) = CL.*CR.*pL.*pR;
            H(3,4,:) = CL.*CR.*pL.*pR;
        end
    end
end