
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% This code is for the Development of Generalized Potenial%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for Vinyl Bromide%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of epochs
epochs= 20;

% Goal to reach
goal=0;

% Gradient 'cutoff'
min_grad = 1e-10;


mu = 0.001;
mu_dec = 0.1;
mu_inc = 10;
mu_max = 1e10;


% Defining NNs for fc(r), f(r) and f(theta)
S2=15;
net_fc = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');

% All the 2-body terms
net_fr_HH = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');
% net_fr_HH.iw{1,1}= [  40.0000;  -40.0000;   40.0000;  -40.0135;  -39.7665;];
% net_fr_HH.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_fr_HH.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_fr_HH.b{2,1}= [-0.1620]; %#ok<NBRAK>
net_fr_HBr = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');
% net_fr_HBr.iw{1,1}= [  40.0000;  -40.0000;   40.0000;  -40.0135;  -39.7665;];
% net_fr_HBr.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_fr_HBr.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_fr_HBr.b{2,1}= [-0.1620]; %#ok<NBRAK>
net_fr_HC = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');
% net_fr_HC.iw{1,1}= [  40.0000;  -40.0000;   40.0000;  -40.0135;  -39.7665;];
% net_fr_HC.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_fr_HC.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_fr_HC.b{2,1}= [-0.1620]; %#ok<NBRAK>
net_fr_CC = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');
% net_fr_CC.iw{1,1}= [  40.0000;  -40.0000;   40.0000;  -40.0135;  -39.7665;];
% net_fr_CC.lw{2,1}=  [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_fr_CC.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_fr_CC.b{2,1}= [-0.1620]; %#ok<NBRAK>
net_fr_CBr = newff([1.2 4.7],[S2 1],{'logsig' 'purelin'},'trainlm');
% net_fr_CBr.iw{1,1}= [  40.0000;  -40.0000;   40.0000;  -40.0135;  -39.7665;];
% net_fr_CBr.lw{2,1}=  [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_fr_CBr.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_fr_CBr.b{2,1}= [-0.1620]; %#ok<NBRAK>

% All the 3-body terms
S3=25;
net_ftheta_HHH = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7 ],[S3 1],{'logsig' 'purelin'},'trainlm');
% net_ftheta_HHH.iw{1,1}= [40.0000  40.0000 40.0000;  -40.0000 -40.0000 -40.0000;   40.0000 40.0000 40.0000;  -40.0135 -40.0135 -40.0135;  -39.7665 -39.7665 -39.7665;];
% net_ftheta_HHH.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_HHH.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_ftheta_HHH.b{2,1}= [-0.1620];

net_ftheta_HHC = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7 ],[S3 1],{'logsig' 'purelin'},'trainlm');
% net_ftheta_HHC.iw{1,1}= [40.0000  40.0000 40.0000;  -40.0000 -40.0000 -40.0000;   40.0000 40.0000 40.0000;  -40.0135 -40.0135 -40.0135;  -39.7665 -39.7665 -39.7665;];
% net_ftheta_HHC.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_HHC.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_ftheta_HHC.b{2,1}= [-0.1620];

net_ftheta_HHBr = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7 ],[S3 1],{'logsig' 'purelin'},'trainlm');
% net_ftheta_HHBr.iw{1,1}= [40.0000  40.0000 40.0000;  -40.0000 -40.0000 -40.0000;   40.0000 40.0000 40.0000;  -40.0135 -40.0135 -40.0135;  -39.7665 -39.7665 -39.7665;];
% net_ftheta_HHBr.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_HHBr.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_ftheta_HHBr.b{2,1}= [-0.1620];

net_ftheta_HCBr = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7 ],[S3 1],{'logsig' 'purelin'},'trainlm');
% net_ftheta_HCBr.iw{1,1}= [40.0000  40.0000 40.0000;  -40.0000 -40.0000 -40.0000;   40.0000 40.0000 40.0000;  -40.0135 -40.0135 -40.0135;  -39.7665 -39.7665 -39.7665;];
% net_ftheta_HCBr.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_HCBr.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_ftheta_HCBr.b{2,1}= [-0.1620];

net_ftheta_HCC = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7 ],[S3 1],{'logsig' 'purelin'},'trainlm');
% net_ftheta_HCC.iw{1,1}= [40.0000  40.0000 40.0000;  -40.0000 -40.0000 -40.0000;   40.0000 40.0000 40.0000;  -40.0135 -40.0135 -40.0135;  -39.7665 -39.7665 -39.7665;];
% net_ftheta_HCC.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_HCC.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_ftheta_HCC.b{2,1}= [-0.1620];

net_ftheta_CCBr = newff([1.2 4.7; 1.2 4.7 ;1.2 4.7 ],[S3 1],{'logsig' 'purelin'},'trainlm');
% net_ftheta_CCBr.iw{1,1}= [40.0000  40.0000 40.0000;  -40.0000 -40.0000 -40.0000;   40.0000 40.0000 40.0000;  -40.0135 -40.0135 -40.0135;  -39.7665 -39.7665 -39.7665;];
% net_ftheta_CCBr.lw{2,1}= [-0.1146   -0.7523    0.9238   -0.9597    0.0854];
% net_ftheta_CCBr.b{1,1}= [ -188.0000;  182.1667; -176.3333;  170.4967;  164.7197;];
% net_ftheta_CCBr.b{2,1}= [-0.1620];

% All the 4-body terms
S4=40;
net_fphi_HHHC = newff([1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7],[S4 1],{'logsig' 'purelin'},'trainlm');
% net_fphi_HHHC.iw{1,1} = [40.0000  40.0000 40.0000 40.0000  40.0000 40.0000;
%                          -40.0000 -40.0000 -40.0000 -40.0000 -40.0000 -40.0000;
%                          40.0000 40.0000 40.0000 40.0000 40.0000 40.0000;  
%                          -40.0135 -40.0135 -40.0135 -40.0135 -40.0135 -40.0135;
%                          -39.7665 -39.7665 -39.7665 -39.7665 -39.7665 -39.7665;];
net_fphi_HHHBr = newff([1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7],[S4 1],{'logsig' 'purelin'},'trainlm');
% net_fphi_HHHBr.iw{1,1} = [40.0000  40.0000 40.0000 40.0000  40.0000 40.0000;
%                          -40.0000 -40.0000 -40.0000 -40.0000 -40.0000 -40.0000;
%                          40.0000 40.0000 40.0000 40.0000 40.0000 40.0000;  
%                          -40.0135 -40.0135 -40.0135 -40.0135 -40.0135 -40.0135;
%                          -39.7665 -39.7665 -39.7665 -39.7665 -39.7665 -39.7665;];
net_fphi_HHCC = newff([1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7],[S4 1],{'logsig' 'purelin'},'trainlm');
% net_fphi_HHCC.iw{1,1} = [40.0000  40.0000 40.0000 40.0000  40.0000 40.0000;
%                          -40.0000 -40.0000 -40.0000 -40.0000 -40.0000 -40.0000;
%                          40.0000 40.0000 40.0000 40.0000 40.0000 40.0000;  
%                          -40.0135 -40.0135 -40.0135 -40.0135 -40.0135 -40.0135;
%                          -39.7665 -39.7665 -39.7665 -39.7665 -39.7665 -39.7665;];
net_fphi_HHCBr = newff([1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7],[S4 1],{'logsig' 'purelin'},'trainlm');
% net_fphi_HHCBr.iw{1,1} = [40.0000  40.0000 40.0000 40.0000  40.0000 40.0000;
%                          -40.0000 -40.0000 -40.0000 -40.0000 -40.0000 -40.0000;
%                          40.0000 40.0000 40.0000 40.0000 40.0000 40.0000;  
%                          -40.0135 -40.0135 -40.0135 -40.0135 -40.0135 -40.0135;
%                          -39.7665 -39.7665 -39.7665 -39.7665 -39.7665 -39.7665;];
net_fphi_HCCBr = newff([1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7; 1.2 4.7],[S4 1],{'logsig' 'purelin'},'trainlm');
% net_fphi_HCCBr.iw{1,1} = [40.0000  40.0000 40.0000 40.0000  40.0000 40.0000;
%                          -40.0000 -40.0000 -40.0000 -40.0000 -40.0000 -40.0000;
%                          40.0000 40.0000 40.0000 40.0000 40.0000 40.0000;  
%                          -40.0135 -40.0135 -40.0135 -40.0135 -40.0135 -40.0135;
%                          -39.7665 -39.7665 -39.7665 -39.7665 -39.7665 -39.7665;];

data = xlsread('vinylBromide_6000points.xlsx');

%no of data pts
Q=length(data);
% Q=1;

%in future clust_size and type will be read from the data file
total_columns=length(data(1,:));

% column 1 in the data file -> cluster size
max_clust_size=max(data(1:Q,1)); %5


% For now, all the atoms are of same type. In future it will be for different types
type(1:max_clust_size)=[6,6,1,1,1,35];


for iQ=1:1:Q
        
    clust_size(iQ)=data(iQ,1);     %#ok<AGROW>
    colcount=2;
    for it=1:1:clust_size(iQ)
        for jt=it:1:clust_size(iQ)
            if jt>it
                
                rs(it,jt,iQ)=data(iQ,colcount); %#ok<AGROW>
                rs(jt,it,iQ)=rs(it,jt,iQ); %#ok<AGROW>
                colcount=colcount+1;
            end
            
        end
    end
    V(iQ)=data(iQ,total_columns); %#ok<AGROW>
end


X_fc = getx(net_fc);

X_fr_HH = getx(net_fr_HH); X_fr_HBr = getx(net_fr_HBr); X_fr_HC = getx(net_fr_HC); X_fr_CC = getx(net_fr_CC); X_fr_CBr = getx(net_fr_CBr);

X_ftheta_HHH = getx(net_ftheta_HHH); X_ftheta_HHC = getx(net_ftheta_HHC); X_ftheta_HHBr = getx(net_ftheta_HHBr); X_ftheta_HCBr = getx(net_ftheta_HCBr); X_ftheta_HCC = getx(net_ftheta_HCC); X_ftheta_CCBr=getx(net_ftheta_CCBr);

X_fphi_HHHC = getx(net_fphi_HHHC); X_fphi_HHHBr = getx(net_fphi_HHHBr); X_fphi_HHCC = getx(net_fphi_HHCC); X_fphi_HHCBr = getx(net_fphi_HHCBr); X_fphi_HCCBr = getx(net_fphi_HCCBr);




% X(1:3*S1+1)=X_fc;
% X(3*S1+1+1:3*S1+1+3*S1+1)=X_fr(1:3*S1+1);
% X(3*S1+1+3*S1+1+1:3*S1+1+3*S1+1+5*S1+1)=X_ftheta(1:5*S1+1);

% X(1:3*S1+1)=X_fc;
% X_fr=vertcat(X_fr_HH, X_fr_HBr, X_fr_HC, X_fr_CC, X_fr_CBr); 
% X(1:5*(3*S1+1),1)=X_fr(1:5*(3*S1+1));
% X_ftheta=vertcat(X_ftheta_HHH, X_ftheta_HHC, X_net_ftheta_HHBr, X_ftheta_HCBr, X_ftheta_HCC, X_ftheta_CCBr);
% X(5*(3*S1+1)+1:5*(3*S1+1)+6*(5*S1+1),1)=X_ftheta(1:6*(5*S1+1));
% X_fphi=vertcat(X_fphi_HHHC, X_fphi_HHHBr, X_fphi_HHCC, X_fphi_HHCBr, X_fphi_HCCBr);
% X(5*(3*S1+1)+6*(5*S1+1)+1:5*(3*S1+1)+6*(5*S1+1)+5*(8*S1+1),1)=X_fphi(1:5*(8*S1+1));

X_fr=vertcat(X_fr_HH, X_fr_HBr, X_fr_HC, X_fr_CC, X_fr_CBr); 
X(1:length(X_fr),1)=X_fr;
X_ftheta=vertcat(X_ftheta_HHH, X_ftheta_HHC, X_ftheta_HHBr, X_ftheta_HCBr, X_ftheta_HCC, X_ftheta_CCBr);
X(length(X_fr)+1:length(X_fr)+ length(X_ftheta),1)=X_ftheta;
X_fphi=vertcat(X_fphi_HHHC, X_fphi_HHHBr, X_fphi_HHCC, X_fphi_HHCBr, X_fphi_HCCBr);
% X(length(X_fr)+ length(X_ftheta)+1:length(X_fr)+
% length(X_ftheta)+length(X_fphi),1)=X_fphi;


numParameters=length(X_fr)+length(X_ftheta);
% numParameters=length(X_fr)+ length(X_ftheta) + length(X_fphi);

ii = sparse(1:numParameters,1:numParameters,ones(1,numParameters));%%%

[perf,Ex, Vhat] = calcperf_GPES(net_fc, net_fr_HH, net_fr_HBr, net_fr_HC, net_fr_CC, net_fr_CBr,...
                                net_ftheta_HHH, net_ftheta_HHC, net_ftheta_HHBr, net_ftheta_HCBr, net_ftheta_HCC, net_ftheta_CCBr,...
                                net_fphi_HHHC, net_fphi_HHHBr, net_fphi_HHCC,net_fphi_HHCBr, net_fphi_HCCBr, rs, V, Q, clust_size, type);


% [gXt,jjt,normgX]=calcjejj_GPES(net_fc, net_fr_HH, net_fr_HBr, net_fr_HC, net_fr_CC, net_fr_CBr,...
%                                 net_ftheta_HHH, net_ftheta_HHC, net_ftheta_HHBr, net_ftheta_HCBr, net_ftheta_HCC, net_ftheta_CCBr,...
%                                 net_fphi_HHHC, net_fphi_HHHBr, net_fphi_HHCC,net_fphi_HHCBr, net_fphi_HCCBr, rs, clust_size,Q,type,Ex);%%%%


for epoch=0:epochs
    
    [je,jj,normgX]=calcjejj_GPES(net_fc, net_fr_HH, net_fr_HBr, net_fr_HC, net_fr_CC, net_fr_CBr,...
                                net_ftheta_HHH, net_ftheta_HHC, net_ftheta_HHBr, net_ftheta_HCBr, net_ftheta_HCC, net_ftheta_CCBr,...
                                net_fphi_HHHC, net_fphi_HHHBr, net_fphi_HHCC,net_fphi_HHCBr, net_fphi_HCCBr, rs, clust_size,Q,type,Ex);%%%%
    
    
    normgX;
  
    if(isnan(normgX) == 1)
        normgX
    end

    % Training Record
    epochPlus1 = epoch+1;
    tr.perf(epochPlus1) = perf;
    tr.mu(epochPlus1) = mu;
    tr.gradient(epochPlus1) = normgX;

    % Stopping Criteria
    %   currentTime = etime(clock,startTime);
    if (perf <= goal)
        stop = 'Performance goal met.';
    elseif (epoch == epochs)
        stop = 'Maximum epoch reached, performance goal was not met.';
        %   elseif (currentTime > time)
        %     stop = 'Maximum time elapsed, performance goal was not met.';
    elseif (normgX < min_grad)
        stop = 'Minimum gradient reached, performance goal was not met.';
    elseif (mu > mu_max)
        stop = 'Maximum MU reached, performance goal was not met.';
        %   elseif (doValidation) & (VV.numFail > max_fail)
        %     stop = 'Validation stop.';
        %   elseif flag_stop
        %     stop = 'User stop.';
    end

    % Progress
    %   if isfinite(show) & (~rem(epoch,show) | length(stop))
    %     fprintf('%s%s%s',this,'-',gradientFcn);
    if isfinite(epochs) fprintf(', Epoch %g/%g',epoch, epochs); end
    %   if isfinite(time) fprintf(', Time %4.1f%%',currentTime/time*100); end
    
    if isfinite(goal) fprintf(', %s %g/%g',upper(net_fc.performFcn),perf,goal); end
    
      
    
    if isfinite(min_grad) fprintf(', Gradient %g/%g',normgX,min_grad); end
    fprintf('\n')
    %   flag_stop=plotperf(tr,goal,this,epoch);
    %     if length(stop) fprintf('%s, %s\n\n',this,stop); end
    %   end

    % Stop when criteria indicate its time
    %   if length(stop)
    %     if (doValidation)
    %     net = VV.net;
    %   end
    %     break
    %   end

    % Levenberg Marquardt
    while (mu <= mu_max)
        % CHECK FOR SINGULAR MATRIX
        [msgstr,msgid] = lastwarn;
        lastwarn('MATLAB:nothing','MATLAB:nothing')
        warnstate = warning('off','all');
        dX = -(jj+ii*mu) \ je;
        [msgstr1,msgid1] = lastwarn;
        flag_inv = isequal(msgid1,'MATLAB:nothing');
        if flag_inv, lastwarn(msgstr,msgid); end;
        warning(warnstate)
        X2 = X + dX;

%         net2_fc = setx(net_fc,X2(1:3*S1+1));
%         net2_fr = setx(net_fr,X2(3*S1+1+1:3*S1+1+3*S1+1));
%         net2_ftheta = setx(net_ftheta,X2(3*S1+1+3*S1+1+1:3*S1+1+3*S1+1+5*S1+1));

        net2_fc = net_fc;
        
        X2_fr = X2(1:length(X_fr),1);
        net2_fr_HH = setx(net_fr_HH,X2_fr(1:length(X_fr_HH)));
        net2_fr_HBr= setx(net_fr_HBr,X2_fr(length(X_fr_HH)+1:length(X_fr_HH)+length(X_fr_HBr)));
        net2_fr_HC = setx(net_fr_HC,X2_fr(length(X_fr_HH)+length(X_fr_HBr)+1:length(X_fr_HH)+length(X_fr_HBr)+length(X_fr_HC))); 
        net2_fr_CC = setx(net_fr_CC,X2_fr(length(X_fr_HH)+length(X_fr_HBr)+length(X_fr_HC)+1:length(X_fr_HH)+length(X_fr_HBr)+length(X_fr_HC)+length(X_fr_CC))); 
        net2_fr_CBr= setx(net_fr_CBr,X2_fr(length(X_fr_HH)+length(X_fr_HBr)+length(X_fr_HC)+length(X_fr_CC)+1:length(X_fr_HH)+length(X_fr_HBr)+length(X_fr_HC)+length(X_fr_CC)+length(X_fr_CBr))); 
        
        
        X2_ftheta = X2(length(X_fr)+1:length(X_fr)+ length(X_ftheta),1);
        net2_ftheta_HHH=setx(net_ftheta_HHH,X2_ftheta(1:length(X_ftheta_HHH)));
        net2_ftheta_HHC=setx(net_ftheta_HHC,X2_ftheta(length(X_ftheta_HHH)+1:length(X_ftheta_HHH)+length(X_ftheta_HHC)));
        net2_ftheta_HHBr=setx(net_ftheta_HHBr,X2_ftheta(length(X_ftheta_HHH)+length(X_ftheta_HHC)+1:length(X_ftheta_HHH)+length(X_ftheta_HHC)+length(X_ftheta_HHBr)));
        net2_ftheta_HCBr=setx(net_ftheta_HCBr,X2_ftheta(length(X_ftheta_HHH)+length(X_ftheta_HHC)+length(X_ftheta_HHBr)+1:length(X_ftheta_HHH)+length(X_ftheta_HHC)+length(X_ftheta_HHBr)+length(X_ftheta_HCBr)));
        net2_ftheta_HCC=setx(net_ftheta_HCC,X2_ftheta(length(X_ftheta_HHH)+length(X_ftheta_HHC)+length(X_ftheta_HHBr)+length(X_ftheta_HCBr)+1:length(X_ftheta_HHH)+length(X_ftheta_HHC)+length(X_ftheta_HHBr)+length(X_ftheta_HCBr)+length(X_ftheta_HCC)));
        net2_ftheta_CCBr=setx(net_ftheta_HCC,X2_ftheta(length(X_ftheta_HHH)+length(X_ftheta_HHC)+length(X_ftheta_HHBr)+length(X_ftheta_HCBr)+length(X_ftheta_HCC)+1:length(X_ftheta_HHH)+length(X_ftheta_HHC)+length(X_ftheta_HHBr)+length(X_ftheta_HCBr)+length(X_ftheta_HCC)+length(X_ftheta_CCBr)));
        
        
%         X2_fphi = X2(length(X_fr)+ length(X_ftheta)+1:length(X_fr)+ length(X_ftheta)+length(X_fphi),1);
%         net2_fphi_HHHC = setx(net_fphi_HHHC,X2_fphi(1:length(X_fphi_HHHC)));
%         net2_fphi_HHHBr= setx(net_fphi_HHHBr,X2_fphi(length(X_fphi_HHHC)+1:length(X_fphi_HHHC)+length(X_fphi_HHHBr)));
%         net2_fphi_HHCC = setx(net_fphi_HHCC,X2_fphi(length(X_fphi_HHHC)+length(X_fphi_HHHBr)+1:length(X_fphi_HHHC)+length(X_fphi_HHHBr)+length(X_fphi_HHCC)));
%         net2_fphi_HHCBr= setx(net_fphi_HHCBr,X2_fphi(length(X_fphi_HHHC)+length(X_fphi_HHHBr)+length(X_fphi_HHCC)+1:length(X_fphi_HHHC)+length(X_fphi_HHHBr)+length(X_fphi_HHCC)+length(X_fphi_HHCBr)));
%         net2_fphi_HCCBr= setx(net_fphi_HCCBr,X2_fphi(length(X_fphi_HHHC)+length(X_fphi_HHHBr)+length(X_fphi_HHCC)+length(X_fphi_HHCBr)+1:length(X_fphi_HHHC)+length(X_fphi_HHHBr)+length(X_fphi_HHCC)+length(X_fphi_HHCBr)+length(X_fphi_HCCBr)));

            net2_fphi_HHHC = net_fphi_HHHC;
            net2_fphi_HHHBr= net_fphi_HHHBr;
            net2_fphi_HHCC = net_fphi_HHCC;
            net2_fphi_HHCBr = net_fphi_HHCBr;
            net2_fphi_HCCBr = net_fphi_HCCBr;

        
        [perf2,Ex, Vhat] = calcperf_GPES(net2_fc, net2_fr_HH, net2_fr_HBr, net2_fr_HC, net2_fr_CC, net2_fr_CBr,...
                                net2_ftheta_HHH, net2_ftheta_HHC, net2_ftheta_HHBr, net2_ftheta_HCBr, net2_ftheta_HCC, net2_ftheta_CCBr,...
                                net2_fphi_HHHC, net2_fphi_HHHBr, net2_fphi_HHCC,net2_fphi_HHCBr, net2_fphi_HCCBr, rs, V, Q, clust_size, type);
        V-Vhat;
        if (perf2 < perf) && flag_inv
%             X = X2; net = net2; %Zb = Zb2; Zi = Zi2; Zl = Zl2;
            X = X2; 
            
            net_fc = net2_fc; 
            
            net_fr_HH = net2_fr_HH;
            net_fr_HBr= net2_fr_HBr;
            net_fr_HC = net2_fr_HC; 
            net_fr_CC = net2_fr_CC; 
            net_fr_CBr= net2_fr_CBr; 
            
            net_ftheta_HHH = net2_ftheta_HHH;
            net_ftheta_HHC = net2_ftheta_HHC;
            net_ftheta_HHBr = net2_ftheta_HHBr;
            net_ftheta_HCBr = net2_ftheta_HCBr;
            net_ftheta_HCC = net2_ftheta_HCC;
            net_ftheta_CCBr = net2_ftheta_CCBr; 
            
            net_fphi_HHHC = net2_fphi_HHHC;
            net_fphi_HHHBr= net2_fphi_HHHBr;
            net_fphi_HHCC = net2_fphi_HHCC;
            net_fphi_HHCBr = net2_fphi_HHCBr;
            net_fphi_HCCBr = net2_fphi_HCCBr;
        
            %N = N2; Ac = Ac2; El = El2;
            perf = perf2;
            mu = mu * mu_dec;
            if (mu < 1e-20)
                mu = 1e-20;
            end
            break   % Must be after the IF
        end
        mu = mu * mu_inc;
    end

end


