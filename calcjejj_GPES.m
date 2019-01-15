%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% This code is for the Development of Generalized Potenial%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for Vinyl Bromide%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [JE,JtJ,normJE]=calcjejj_GPES(net_fc, net_fr_HH, net_fr_HBr, net_fr_HC, net_fr_CC, net_fr_CBr,...
                                net_ftheta_HHH, net_ftheta_HHC, net_ftheta_HHBr, net_ftheta_HCBr, net_ftheta_HCC, net_ftheta_CCBr,...
                                net_fphi_HHHC, net_fphi_HHHBr, net_fphi_HHCC,net_fphi_HHCBr, net_fphi_HCCBr, rs, clust_size, Q, type, Ex)

% There will be 3 Jacobian matrices, one each for fc, fr, and ftheta

%%%% Jacobian matrix for fr
%J_fr will have dimensions of Q x ( 3N + 1 )
% First layer W (N x 1), b (N x 1), Second Layer W (1 x N ), b (1x1)
W1HH = net_fr_HH.IW{1,1};
W2HH = net_fr_HH.LW{2,1};
B1HH = net_fr_HH.b{1};
B2HH = net_fr_HH.b{2};

W1HBr = net_fr_HBr.IW{1,1};
W2HBr = net_fr_HBr.LW{2,1};
B1HBr = net_fr_HBr.b{1};
B2HBr = net_fr_HBr.b{2};

W1HC = net_fr_HC.IW{1,1};
W2HC = net_fr_HC.LW{2,1};
B1HC = net_fr_HC.b{1};
B2HC = net_fr_HC.b{2};

W1CC = net_fr_CC.IW{1,1};
W2CC = net_fr_CC.LW{2,1};
B1CC = net_fr_CC.b{1};
B2CC = net_fr_CC.b{2};

W1CBr = net_fr_CBr.IW{1,1};
W2CBr = net_fr_CBr.LW{2,1};
B1CBr = net_fr_CBr.b{1};
B2CBr = net_fr_CBr.b{2};




S1 = 15;
S2 = 1;

% %%%%%%%
for iQ = 1:Q

    ciQ = 0;
    ciQHH = 0; ciQHBr = 0; ciQHC = 0; ciQCC = 0; ciQCBr = 0;

    for i=1:1:clust_size(iQ)
        for j=1:1:clust_size(iQ)
            if j>i%j~= i %
                ciQ = ciQ+1;
                identify2=sprintf('%d %d', type(i),type(j));
                switch identify2
                   case {'1 1'}
                      ciQHH = ciQHH+1;
                      r_fr_HH(ciQHH) = rs(i,j,iQ);                      
                   case {'1 35', '35 1'}
                      ciQHBr = ciQHBr+1;
                      r_fr_HBr(ciQHBr) = rs(i,j,iQ);
                   case {'1 6', '6 1'}
                      ciQHC = ciQHC+1;
                      r_fr_HC(ciQHC) = rs(i,j,iQ);
                   case {'6 6'}
                      ciQCC = ciQCC+1;
                      r_fr_CC(ciQCC) = rs(i,j,iQ);
                    case {'6 35', '35 6'}
                      ciQCBr = ciQCBr+1;
                      r_fr_CBr(ciQCBr) = rs(i,j,iQ);
                end
                r_fr(ciQ) = rs(i,j,iQ);

            end
        end
    end




    A1HH = logsig(W1HH*r_fr_HH + B1HH*ones(1,(3)*1));
    A2HH = W2HH*A1HH + B2HH*ones(1,(3)*1);
    % FIND JACOBIAN
    A1HH = kron(A1HH,ones(1,S2));
    D2HH = nnmdlin(A2HH);

    D1HH = nnmdlog(A1HH,D2HH,W2HH);
    jac1 = nnlmarq(kron(r_fr_HH,ones(1,S2)),D1HH);
    jac2 = nnlmarq(A1HH,D2HH);

    jac_fr_HH = [jac1,D1HH',jac2,D2HH']; % Make jacbian for current iQ configuration
    
    for j= 1: 3*S1+1
        jac(1,j)= 0;
        for i = 1:(3)
            jac(1,j) = jac(1,j)+jac_fr_HH(i,j);  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    
    A1HBr = logsig(W1HBr*r_fr_HBr + B1HBr*ones(1,(3)*1));
    A2HBr = W2HBr*A1HBr + B2HBr*ones(1,(3)*1);
    % FIND JACOBIAN
    A1HBr = kron(A1HBr,ones(1,S2));
    D2HBr = nnmdlin(A2HBr);

    D1HBr = nnmdlog(A1HBr,D2HBr,W2HBr);
    jac1 = nnlmarq(kron(r_fr_HBr,ones(1,S2)),D1HBr);
    jac2 = nnlmarq(A1HBr,D2HBr);

    jac_fr_HBr = [jac1,D1HBr',jac2,D2HBr']; % Make jacbian for current iQ configuration

    for j= (3*S1+1)+1: 2*(3*S1+1)
        jac(1,j)= 0;
        for i = 1:(3)
            jac(1,j) = jac(1,j)+jac_fr_HBr(i,j-(3*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    
    
    A1HC = logsig(W1HC*r_fr_HC + B1HC*ones(1,(6)*1));
    A2HC = W2HC*A1HC + B2HC*ones(1,(6)*1);
    % FIND JACOBIAN
    A1HC = kron(A1HC,ones(1,S2));
    D2HC = nnmdlin(A2HC);

    D1HC = nnmdlog(A1HC,D2HC,W2HC);
    jac1 = nnlmarq(kron(r_fr_HC,ones(1,S2)),D1HC);
    jac2 = nnlmarq(A1HC,D2HC);

    jac_fr_HC = [jac1,D1HC',jac2,D2HC']; % Make jacbian for current iQ configuration

    for j= 2*(3*S1+1)+1: 3*(3*S1+1)
        jac(1,j)= 0;
        for i = 1:(6)
            jac(1,j) = jac(1,j)+jac_fr_HC(i,j-2*(3*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    
    
    A1CC = logsig(W1CC*r_fr_CC + B1CC*ones(1,(1)*1));
    A2CC = W2CC*A1CC + B2CC*ones(1,(1)*1);
    % FIND JACOBIAN
    A1CC = kron(A1CC,ones(1,S2));
    D2CC = nnmdlin(A2CC);

    D1CC = nnmdlog(A1CC,D2CC,W2CC);
    jac1 = nnlmarq(kron(r_fr_CC,ones(1,S2)),D1CC);
    jac2 = nnlmarq(A1CC,D2CC);

    jac_fr_CC = [jac1,D1CC',jac2,D2CC']; % Make jacbian for current iQ configuration

    for j= 3*(3*S1+1)+1: 4*(3*S1+1)
        jac(1,j)= 0;
        for i = 1:(1)
            jac(1,j) = jac(1,j)+jac_fr_CC(i,j-3*(3*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    
    
    A1CBr = logsig(W1CBr*r_fr_CBr + B1CBr*ones(1,(2)*1));
    A2CBr = W2CBr*A1CBr + B2CBr*ones(1,(2)*1);
    % FIND JACOBIAN
    A1CBr = kron(A1CBr,ones(1,S2));
    D2CBr = nnmdlin(A2CBr);

    D1CBr = nnmdlog(A1CBr,D2CBr,W2CBr);
    jac1 = nnlmarq(kron(r_fr_CBr,ones(1,S2)),D1CBr);
    jac2 = nnlmarq(A1CBr,D2CBr);

    jac_fr_CBr = [jac1,D1CBr',jac2,D2CBr']; % Make jacbian for current iQ configuration

    for j= 4*(3*S1+1)+1: 5*(3*S1+1)
        jac(1,j)= 0;
        for i = 1:(2)
            jac(1,j) = jac(1,j)+jac_fr_CBr(i,j-4*(3*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    
   

    for j= 1:5*(3*S1+1)
        J_fr(iQ,j)= jac(1,j);
    end

end


% %%%%%%%


% J_theta will have dimensions of Q x ( 5N + 1)
% First layer W (N x 3), b (N x 1), Second Layer W (1 x N ), b (1x1)

W1HHH = net_ftheta_HHH.IW{1,1};
W2HHH = net_ftheta_HHH.LW{2,1};
B1HHH = net_ftheta_HHH.b{1};
B2HHH = net_ftheta_HHH.b{2};

W1HHC = net_ftheta_HHC.IW{1,1};
W2HHC = net_ftheta_HHC.LW{2,1};
B1HHC = net_ftheta_HHC.b{1};
B2HHC = net_ftheta_HHC.b{2};

W1HHBr = net_ftheta_HHBr.IW{1,1};
W2HHBr = net_ftheta_HHBr.LW{2,1};
B1HHBr = net_ftheta_HHBr.b{1};
B2HHBr = net_ftheta_HHBr.b{2};

W1HCBr = net_ftheta_HCBr.IW{1,1};
W2HCBr = net_ftheta_HCBr.LW{2,1};
B1HCBr = net_ftheta_HCBr.b{1};
B2HCBr = net_ftheta_HCBr.b{2};

W1HCC = net_ftheta_HCC.IW{1,1};
W2HCC = net_ftheta_HCC.LW{2,1};
B1HCC = net_ftheta_HCC.b{1};
B2HCC = net_ftheta_HCC.b{2};

W1CCBr = net_ftheta_CCBr.IW{1,1};
W2CCBr = net_ftheta_CCBr.LW{2,1};
B1CCBr = net_ftheta_CCBr.b{1};
B2CCBr = net_ftheta_CCBr.b{2};



S1 = 25;
S2 = 1;




for iQ = 1:Q
    ciQ = 0;
    ciQHHH=0;ciQHHC=0;ciQHHBr=0;ciQHCBr=0;ciQHCC=0;ciQCCBr=0;
    for i=1:1:clust_size(iQ)
        for j=1:1:clust_size(iQ)
            for k=1:1:clust_size(iQ)
                if (j>i && k>i && k>j)
                    ciQ = ciQ+1;
                    identify3=sprintf('%d %d %d', type(i),type(j),type(k));
                    switch identify3
                       case {'1 1 1'}
                           ciQHHH=ciQHHH+1;
                           r_ftheta_HHH(1,ciQHHH)=rs(i,j,iQ);
                           r_ftheta_HHH(2,ciQHHH)=rs(i,k,iQ);
                           r_ftheta_HHH(3,ciQHHH)=rs(j,k,iQ);
                       case {'1 1 6', '1 6 1', '6 1 1'}
                           ciQHHC=ciQHHC+1;
                           r_ftheta_HHC(1,ciQHHC)=rs(i,j,iQ);
                           r_ftheta_HHC(2,ciQHHC)=rs(i,k,iQ);
                           r_ftheta_HHC(3,ciQHHC)=rs(j,k,iQ);
                          
                      case {'1 1 35','1 35 1','35 1 1'}
                           ciQHHBr=ciQHHBr+1;
                           r_ftheta_HHBr(1,ciQHHBr)=rs(i,j,iQ);
                           r_ftheta_HHBr(2,ciQHHBr)=rs(i,k,iQ);
                           r_ftheta_HHBr(3,ciQHHBr)=rs(j,k,iQ);
                          
                      case {'1 6 35','1 35 6', '6 1 35', '6 35 1', '35 1 6', '35 6 1'}
                           ciQHCBr=ciQHCBr+1;
                           r_ftheta_HCBr(1,ciQHCBr)=rs(i,j,iQ);
                           r_ftheta_HCBr(2,ciQHCBr)=rs(i,k,iQ);
                           r_ftheta_HCBr(3,ciQHCBr)=rs(j,k,iQ);
                      case {'1 6 6','6 6 1', '6 1 6'}
                           ciQHCC=ciQHCC+1;
                           r_ftheta_HCC(1,ciQHCC)=rs(i,j,iQ);
                           r_ftheta_HCC(2,ciQHCC)=rs(i,k,iQ);
                           r_ftheta_HCC(3,ciQHCC)=rs(j,k,iQ);
                      case {'6 6 35', '6 35 6', '35 6 6'}
                           ciQCCBr=ciQCCBr+1;
                           r_ftheta_CCBr(1,ciQCCBr)=rs(i,j,iQ);
                           r_ftheta_CCBr(2,ciQCCBr)=rs(i,k,iQ);
                           r_ftheta_CCBr(3,ciQCCBr)=rs(j,k,iQ); 
                    end
%                     r_ftheta(1,ciQ) = rs(i,j,iQ);
%                     r_ftheta(2,ciQ) = rs(i,k,iQ);
%                     r_ftheta(3,ciQ) = rs(j,k,iQ);

                end
            end
        end
    end

    A1HHH = logsig(W1HHH*r_ftheta_HHH + B1HHH*ones(1,(1)*1));
    A2HHH = W2HHH*A1HHH+B2HHH*ones(1,(1)*1);

    % FIND JACOBIAN
    A1HHH = kron(A1HHH,ones(1,S2));
    D2HHH = nnmdlin(A2HHH);

    D1HHH = nnmdlog(A1HHH,D2HHH,W2HHH);
    jac1 = nnlmarq(kron(r_ftheta_HHH,ones(1,S2)),D1HHH);
    jac2 = nnlmarq(A1HHH,D2HHH);

    jac_ftheta_HHH = [jac1,D1HHH',jac2,D2HHH'];
    %%%%
    for j= 1: 5*S1+1
        jac(1,j)= 0;
        for i = 1:(1)
            jac(1,j) = jac(1,j)+jac_ftheta_HHH(i,j);  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    
    A1HHC = logsig(W1HHC*r_ftheta_HHC + B1HHC*ones(1,(6)*1));
    A2HHC = W2HHC*A1HHC+B2HHC*ones(1,(6)*1);

    % FIND JACOBIAN
    A1HHC = kron(A1HHC,ones(1,S2));
    D2HHC = nnmdlin(A2HHC);

    D1HHC = nnmdlog(A1HHC,D2HHC,W2HHC);
    jac1 = nnlmarq(kron(r_ftheta_HHC,ones(1,S2)),D1HHC);
    jac2 = nnlmarq(A1HHC,D2HHC);

    jac_ftheta_HHC = [jac1,D1HHC',jac2,D2HHC'];
    %%%%
    for j= (5*S1+1)+1: 2*(5*S1+1)
        jac(1,j)= 0;
        for i = 1:(6)
            jac(1,j) = jac(1,j)+jac_ftheta_HHC(i,j-(5*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    A1HHBr = logsig(W1HHBr*r_ftheta_HHBr + B1HHBr*ones(1,(3)*1));
    A2HHBr = W2HHBr*A1HHBr+B2HHBr*ones(1,(3)*1);

    % FIND JACOBIAN
    A1HHBr = kron(A1HHBr,ones(1,S2));
    D2HHBr = nnmdlin(A2HHBr);

    D1HHBr = nnmdlog(A1HHBr,D2HHBr,W2HHBr);
    jac1 = nnlmarq(kron(r_ftheta_HHBr,ones(1,S2)),D1HHBr);
    jac2 = nnlmarq(A1HHBr,D2HHBr);

    jac_ftheta_HHBr = [jac1,D1HHBr',jac2,D2HHBr'];
    %%%%
    for j= 2*(5*S1+1)+1: 3*(5*S1+1)
        jac(1,j)= 0;
        for i = 1:(3)
            jac(1,j) = jac(1,j)+jac_ftheta_HHBr(i,j-2*(5*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    A1HCBr = logsig(W1HCBr*r_ftheta_HCBr + B1HCBr*ones(1,(6)*1));
    A2HCBr = W2HCBr*A1HCBr+B2HCBr*ones(1,(6)*1);

    % FIND JACOBIAN
    A1HCBr = kron(A1HCBr,ones(1,S2));
    D2HCBr = nnmdlin(A2HCBr);

    D1HCBr = nnmdlog(A1HCBr,D2HCBr,W2HCBr);
    jac1 = nnlmarq(kron(r_ftheta_HCBr,ones(1,S2)),D1HCBr);
    jac2 = nnlmarq(A1HCBr,D2HCBr);

    jac_ftheta_HCBr = [jac1,D1HCBr',jac2,D2HCBr'];
    %%%%
    for j= 3*(5*S1+1)+1: 4*(5*S1+1)
        jac(1,j)= 0;
        for i = 1:(6)
            jac(1,j) = jac(1,j)+jac_ftheta_HCBr(i,j-3*(5*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    A1HCC = logsig(W1HCC*r_ftheta_HCC + B1HCC*ones(1,(3)*1));
    A2HCC = W2HCC*A1HCC+B2HCC*ones(1,(3)*1);

    % FIND JACOBIAN
    A1HCC = kron(A1HCC,ones(1,S2));
    D2HCC = nnmdlin(A2HCC);

    D1HCC = nnmdlog(A1HCC,D2HCC,W2HCC);
    jac1 = nnlmarq(kron(r_ftheta_HCC,ones(1,S2)),D1HCC);
    jac2 = nnlmarq(A1HCC,D2HCC);

    jac_ftheta_HCC = [jac1,D1HCC',jac2,D2HCC'];
    %%%%
    for j= 4*(5*S1+1)+1: 5*(5*S1+1)
        jac(1,j)= 0;
        for i = 1:(3)
            jac(1,j) = jac(1,j)+jac_ftheta_HCC(i,j-4*(5*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
    A1CCBr = logsig(W1CCBr*r_ftheta_CCBr + B1CCBr*ones(1,(1)*1));
    A2CCBr = W2CCBr*A1CCBr+B2CCBr*ones(1,(1)*1);

    % FIND JACOBIAN
    A1CCBr = kron(A1CCBr,ones(1,S2));
    D2CCBr = nnmdlin(A2CCBr);

    D1CCBr = nnmdlog(A1CCBr,D2CCBr,W2CCBr);
    jac1 = nnlmarq(kron(r_ftheta_CCBr,ones(1,S2)),D1CCBr);
    jac2 = nnlmarq(A1CCBr,D2CCBr);

    jac_ftheta_CCBr = [jac1,D1CCBr',jac2,D2CCBr'];
    %%%%
    for j= 5*(5*S1+1)+1: 6*(5*S1+1)
        jac(1,j)= 0;
        for i = 1:(1)
            jac(1,j) = jac(1,j)+jac_ftheta_CCBr(i,j-5*(5*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row

        end
    end
    
     for j= 1: 6*(5*S1+1)
        J_ftheta(iQ,j)= jac(1,j);
    end
    
end




 
    
W1HHHC = net_fphi_HHHC.IW{1,1};
W2HHHC = net_fphi_HHHC.LW{2,1};
B1HHHC = net_fphi_HHHC.b{1};
B2HHHC = net_fphi_HHHC.b{2};

W1HHHBr = net_fphi_HHHBr.IW{1,1};
W2HHHBr = net_fphi_HHHBr.LW{2,1};
B1HHHBr = net_fphi_HHHBr.b{1};
B2HHHBr = net_fphi_HHHBr.b{2};

W1HHCC = net_fphi_HHCC.IW{1,1};
W2HHCC = net_fphi_HHCC.LW{2,1};
B1HHCC = net_fphi_HHCC.b{1};
B2HHCC = net_fphi_HHCC.b{2};

W1HHCBr = net_fphi_HHCBr.IW{1,1};
W2HHCBr = net_fphi_HHCBr.LW{2,1};
B1HHCBr = net_fphi_HHCBr.b{1};
B2HHCBr = net_fphi_HHCBr.b{2};

W1HCCBr = net_fphi_HCCBr.IW{1,1};
W2HCCBr = net_fphi_HCCBr.LW{2,1};
B1HCCBr = net_fphi_HCCBr.b{1};
B2HCCBr = net_fphi_HCCBr.b{2};



S1 = 40;
S2 = 1;



% 
% for iQ = 1:Q
%     if clust_size(iQ)<4
%         for j= 1: 8*S1+1
%             J_fphi(iQ,j)= 0.0;
%         end
%         continue;
%     end
%     ciQ = 0;
%     ciQHHHC=0;ciQHHHBr=0;ciQHHCC=0;ciQHHCBr=0;ciQHCCBr=0;
%     for i=1:1:clust_size(iQ)
%         for j=1:1:clust_size(iQ)
%             for k=1:1:clust_size(iQ)
%                 for l=1:1:clust_size(iQ)
%                     if (j>i && k>i && k>j && l>i && l>j && l>k)
%                         ciQ = ciQ+1;
% 
%                         identify4=sprintf('%d %d %d %d', type(i),type(j),type(k), type(l));
%                         switch identify4
% 
%                             case {'1 1 1 6', '1 1 6 1', '1 6 1 1', '6 1 1 1'}
%                                 ciQHHHC=ciQHHHC+1;
%                                 r_fphi_HHHC(1,ciQHHHC) = rs(i,j,iQ);
%                                 r_fphi_HHHC(2,ciQHHHC) = rs(i,k,iQ);
%                                 r_fphi_HHHC(3,ciQHHHC) = rs(i,l,iQ);
%                                 r_fphi_HHHC(4,ciQHHHC) = rs(j,k,iQ);
%                                 r_fphi_HHHC(5,ciQHHHC) = rs(j,l,iQ);
%                                 r_fphi_HHHC(6,ciQHHHC) = rs(k,l,iQ);
%                             case {'1 1 1 35', '1 1 35 1', '1 35 1 1', '35 1 1 1'}
%                                 ciQHHHBr=ciQHHHBr+1;
%                                 r_fphi_HHHBr(1,ciQHHHBr) = rs(i,j,iQ);
%                                 r_fphi_HHHBr(2,ciQHHHBr) = rs(i,k,iQ);
%                                 r_fphi_HHHBr(3,ciQHHHBr) = rs(i,l,iQ);
%                                 r_fphi_HHHBr(4,ciQHHHBr) = rs(j,k,iQ);
%                                 r_fphi_HHHBr(5,ciQHHHBr) = rs(j,l,iQ);
%                                 r_fphi_HHHBr(6,ciQHHHBr) = rs(k,l,iQ);
%                             case {'1 1 6 6', '1 6 6 1', '6 6 1 1', '1 6 1 6', '6 1 6 1', '6 1 1 6'}
%                                 ciQHHCC=ciQHHCC+1;
%                                 r_fphi_HHCC(1,ciQHHCC) = rs(i,j,iQ);
%                                 r_fphi_HHCC(2,ciQHHCC) = rs(i,k,iQ);
%                                 r_fphi_HHCC(3,ciQHHCC) = rs(i,l,iQ);
%                                 r_fphi_HHCC(4,ciQHHCC) = rs(j,k,iQ);
%                                 r_fphi_HHCC(5,ciQHHCC) = rs(j,l,iQ);
%                                 r_fphi_HHCC(6,ciQHHCC) = rs(k,l,iQ);
%                             case {'1 1 6 35', '1 1 35 6', '1 6 35 1','1 35 6 1', '6 1 1 35', '35 1 1 6', '6 35 1 1', '35 6 1 1'}
%                                 ciQHHCBr=ciQHHCBr+1;
%                                 r_fphi_HHCBr(1,ciQHHCBr) = rs(i,j,iQ);
%                                 r_fphi_HHCBr(2,ciQHHCBr) = rs(i,k,iQ);
%                                 r_fphi_HHCBr(3,ciQHHCBr) = rs(i,l,iQ);
%                                 r_fphi_HHCBr(4,ciQHHCBr) = rs(j,k,iQ);
%                                 r_fphi_HHCBr(5,ciQHHCBr) = rs(j,l,iQ);
%                                 r_fphi_HHCBr(6,ciQHHCBr) = rs(k,l,iQ);
%                             case {'1 6 6 35', '1 6 35 6', '1 35 6 6', '35 6 6 1', '35 6 1 6', '35 1 6 6', '6 6 1 35'}
%                                 ciQHCCBr=ciQHCCBr+1;
%                                 r_fphi_HCCBr(1,ciQHCCBr) = rs(i,j,iQ);
%                                 r_fphi_HCCBr(2,ciQHCCBr) = rs(i,k,iQ);
%                                 r_fphi_HCCBr(3,ciQHCCBr) = rs(i,l,iQ);
%                                 r_fphi_HCCBr(4,ciQHCCBr) = rs(j,k,iQ);
%                                 r_fphi_HCCBr(5,ciQHCCBr) = rs(j,l,iQ);
%                                 r_fphi_HCCBr(6,ciQHCCBr) = rs(k,l,iQ);
%                         end
%                         %                         r_fphi(1,ciQ) = rs(i,j,iQ);
%                         %                         r_fphi(2,ciQ) = rs(i,k,iQ);
%                         %                         r_fphi(3,ciQ) = rs(i,l,iQ);
%                         %                         r_fphi(4,ciQ) = rs(j,k,iQ);
%                         %                         r_fphi(5,ciQ) = rs(j,l,iQ);
%                         %                         r_fphi(6,ciQ) = rs(k,l,iQ);
% 
%                     end
%                 end
%             end
%         end
%     end
% 
%     A1HHHC = logsig(W1HHHC*r_fphi_HHHC + B1HHHC*ones(1,(2)*1));
%     A2HHHC = W2HHHC*A1HHHC+B2HHHC*ones(1,(2)*1);
% 
%     % FIND JACOBIAN
%     A1HHHC = kron(A1HHHC,ones(1,S2));
%     D2HHHC = nnmdlin(A2HHHC);
% 
% 
%     D1HHHC = nnmdlog(A1HHHC,D2HHHC,W2HHHC);
%     jac1 = nnlmarq(kron(r_fphi_HHHC,ones(1,S2)),D1HHHC);
%     jac2 = nnlmarq(A1HHHC,D2HHHC);
% 
%     jac_fphi_HHHC = [jac1,D1HHHC',jac2,D2HHHC'];
% 
%     for j= 1: 8*S1+1
%         jac(1,j)= 0;
%         for i = 1:(2)
%             jac(1,j) = jac(1,j)+jac_fphi_HHHC(i,j);  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
% 
% 
%     A1HHHBr = logsig(W1HHHBr*r_fphi_HHHBr + B1HHHBr*ones(1,(1)*1));
%     A2HHHBr = W2HHHBr*A1HHHBr+B2HHHBr*ones(1,(1)*1);
% 
%     % FIND JACOBIAN
%     A1HHHBr = kron(A1HHHBr,ones(1,S2));
%     D2HHHBr = nnmdlin(A2HHHBr);
% 
% 
%     D1HHHBr = nnmdlog(A1HHHBr,D2HHHBr,W2HHHBr);
%     jac1 = nnlmarq(kron(r_fphi_HHHBr,ones(1,S2)),D1HHHBr);
%     jac2 = nnlmarq(A1HHHBr,D2HHHBr);
% 
%     jac_fphi_HHHBr = [jac1,D1HHHBr',jac2,D2HHHBr'];
% 
%     for j= (8*S1+1)+1: 2*(8*S1+1)
%         jac(1,j)= 0;
%         for i = 1:(1)
%             jac(1,j) = jac(1,j)+jac_fphi_HHHBr(i,j-(8*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
% 
%     A1HHCC = logsig(W1HHCC*r_fphi_HHCC + B1HHCC*ones(1,(3)*1));
%     A2HHCC = W2HHCC*A1HHCC+B2HHCC*ones(1,(3)*1);
% 
%     % FIND JACOBIAN
%     A1HHCC = kron(A1HHCC,ones(1,S2));
%     D2HHCC = nnmdlin(A2HHCC);
% 
% 
%     D1HHCC = nnmdlog(A1HHCC,D2HHCC,W2HHCC);
%     jac1 = nnlmarq(kron(r_fphi_HHCC,ones(1,S2)),D1HHCC);
%     jac2 = nnlmarq(A1HHCC,D2HHCC);
% 
%     jac_fphi_HHCC = [jac1,D1HHCC',jac2,D2HHCC'];
% 
%     for j= 2*(8*S1+1)+1: 3*(8*S1+1)
%         jac(1,j)= 0;
%         for i = 1:(3)
%             jac(1,j) = jac(1,j)+jac_fphi_HHCC(i,j-2*(8*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
% 
%     A1HHCBr = logsig(W1HHCBr*r_fphi_HHCBr + B1HHCBr*ones(1,(6)*1));
%     A2HHCBr = W2HHCBr*A1HHCBr+B2HHCBr*ones(1,(6)*1);
% 
%     % FIND JACOBIAN
%     A1HHCBr = kron(A1HHCBr,ones(1,S2));
%     D2HHCBr = nnmdlin(A2HHCBr);
% 
% 
%     D1HHCBr = nnmdlog(A1HHCBr,D2HHCBr,W2HHCBr);
%     jac1 = nnlmarq(kron(r_fphi_HHCBr,ones(1,S2)),D1HHCBr);
%     jac2 = nnlmarq(A1HHCBr,D2HHCBr);
% 
%     jac_fphi_HHCBr = [jac1,D1HHCBr',jac2,D2HHCBr'];
% 
%     for j= 3*(8*S1+1)+1: 4*(8*S1+1)
%         jac(1,j)= 0;
%         for i = 1:(6)
%             jac(1,j) = jac(1,j)+jac_fphi_HHCBr(i,j-3*(8*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
% 
% 
%     A1HCCBr = logsig(W1HCCBr*r_fphi_HCCBr + B1HCCBr*ones(1,(3)*1));
%     A2HCCBr = W2HCCBr*A1HCCBr+B2HCCBr*ones(1,(3)*1);
% 
%     % FIND JACOBIAN
%     A1HCCBr = kron(A1HCCBr,ones(1,S2));
%     D2HCCBr = nnmdlin(A2HCCBr);
% 
% 
%     D1HCCBr = nnmdlog(A1HCCBr,D2HCCBr,W2HCCBr);
%     jac1 = nnlmarq(kron(r_fphi_HCCBr,ones(1,S2)),D1HCCBr);
%     jac2 = nnlmarq(A1HCCBr,D2HCCBr);
% 
%     jac_fphi_HCCBr = [jac1,D1HCCBr',jac2,D2HCCBr'];
% 
%     for j= 4*(8*S1+1)+1: 5*(8*S1+1)
%         jac(1,j)= 0;
%         for i = 1:(3)
%             jac(1,j) = jac(1,j)+jac_fphi_HCCBr(i,j-4*(8*S1+1));  % sum jacobian for current iQ configuration along all rows to get 1 row
% 
%         end
%     end
% 
% 
%     for j= 1: 5*(8*S1+1)
%         J_fphi(iQ,j)= jac(1,j);
%     end
% end




%%%Make cobined Jacobian matrix
% J = zeros(Q,5*(3*S1+1)+6*(5*S1+1)+5*(8*S1+1) );

% J(:,1:5*(3*S1+1))=J_fr;
% J(:,5*(3*S1+1)+1:5*(3*S1+1)+6*(5*S1+1))=J_ftheta;
% J(:,(3*S1+1)+(5*S1+1)+1:(3*S1+1)+(5*S1+1)+(8*S1+1))=J_fphi;
% J=horzcat(J_fr,J_ftheta,J_fphi);
J=horzcat(J_fr,J_ftheta);
JE=J'*Ex;
JtJ=J'*J;
normJE=sqrt(JE'*JE);



i=1; %dummy line for breakpoint




