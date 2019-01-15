%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% This code is for the Development of Generalized Potenial%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for Vinyl Bromide%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [perf,Ex,Vhat]=calcperf_GPES(net_fc, net_fr_HH, net_fr_HBr, net_fr_HC, net_fr_CC, net_fr_CBr,...
                                net_ftheta_HHH, net_ftheta_HHC, net_ftheta_HHBr, net_ftheta_HCBr, net_ftheta_HCC, net_ftheta_CCBr,...
                                net_fphi_HHHC, net_fphi_HHHBr, net_fphi_HHCC,net_fphi_HHCBr, net_fphi_HCCBr, rs, V, Q, clust_size, type)

perf=0.0;

for iQ=1:1:Q
    
    %converting rs to my format
    for i=1:1:clust_size(iQ)
         for j=1:1:clust_size(iQ)
             if i~=j
                 r(i,j)=rs(i,j,iQ); %#ok<AGROW>
             end
         end
    end
    
    
    %calculating potential energy for iQ configuration
     
    frHH=0.0; frHBr=0.0; frHC=0.0; frCC=0.0; frCBr=0.0;
    fthetaHHH=0.0; fthetaHHC=0.0; fthetaHHBr=0.0; fthetaHCBr=0.0; fthetaHCC=0.0; fthetaCCBr=0.0;
    fphiHHHC=0.0; fphiHHHBr=0.0; fphiHHCC=0.0; fphiHHCBr=0.0; fphiHCCBr=0.0;
    fc=1;%sim(net_fc,rij);
    for i=1:1:clust_size(iQ)
        for j=1:1:clust_size(iQ)
            if j>i
                rij=r(i,j);
%                 identify=horzcat(type(i),type(j));
                identify2=sprintf('%d %d', type(i),type(j));
                switch identify2
                   case {'1 1'}
                      frHH=frHH+sim(net_fr_HH,rij);
                   case {'1 35', '35 1'}
                      frHBr=frHBr+sim(net_fr_HBr,rij);
                   case {'1 6', '6 1'}
                      frHC=frHC+sim(net_fr_HC,rij);
                   case {'6 6'}
                      frCC=frCC+sim(net_fr_CC,rij);
                    case {'6 35', '35 6'}
                      frCBr=frCBr+sim(net_fr_CBr,rij);   
                    otherwise
                        disp(identify2)
                   
                end
                
            end
            for k=1:1:clust_size(iQ)
                if (j>i && k>i && k>j)
                    rik=r(i,k);
                    rjk=r(j,k);
                    identify3=sprintf('%d %d %d', type(i),type(j),type(k));
                    switch identify3
                       case {'1 1 1'}
                          fthetaHHH=fthetaHHH+sim(net_ftheta_HHH,[rij; rik ;rjk]);
                       case {'1 1 6', '1 6 1', '6 1 1'}
                          fthetaHHC=fthetaHHC+sim(net_ftheta_HHC,[rij; rik ;rjk]);
                      case {'1 1 35','1 35 1','35 1 1'}
                          fthetaHHBr=fthetaHHBr+sim(net_ftheta_HHBr,[rij; rik ;rjk]);
                      case {'1 6 35','1 35 6', '6 1 35', '6 35 1', '35 1 6', '35 6 1'}
                          fthetaHCBr=fthetaHCBr+sim(net_ftheta_HCBr,[rij; rik ;rjk]);
                      case {'1 6 6','6 6 1', '6 1 6'}
                          fthetaHCC=fthetaHCC+sim(net_ftheta_HCC,[rij; rik ;rjk]);
                      case {'6 6 35', '6 35 6', '35 6 6'}
                          fthetaCCBr=fthetaCCBr+sim(net_ftheta_CCBr,[rij; rik ;rjk]);  
                      otherwise
                        disp(identify3)
                    end
                                                      
                end
%                 for l=1:1:clust_size(iQ)
%                     if (j>i && k>i && k>j && l>i && l>j && l>k)
%                         ril=r(i,l);
%                         rjl=r(j,l);
%                         rkl=r(k,l);
%                         identify4=sprintf('%d %d %d %d', type(i),type(j),type(k), type(l));
%                         switch identify4
%                             
%                            case {'1 1 1 6', '1 1 6 1', '1 6 1 1', '6 1 1 1'}
%                               fphiHHHC=fphiHHHC+sim(net_fphi_HHHC,[rij; rik ; ril; rjk; rjl; rkl]);
%                            case {'1 1 1 35', '1 1 35 1', '1 35 1 1', '35 1 1 1'}
%                               fphiHHHBr=fphiHHHBr+sim(net_fphi_HHHBr,[rij; rik ; ril; rjk; rjl; rkl]);
%                            case {'1 1 6 6', '1 6 6 1', '6 6 1 1', '1 6 1 6', '6 1 6 1', '6 1 1 6'}
%                               fphiHHCC=fphiHHCC+sim(net_fphi_HHCC,[rij; rik ; ril; rjk; rjl; rkl]);
%                            case {'1 1 6 35', '1 1 35 6', '1 6 35 1','1 35 6 1', '6 1 1 35', '35 1 1 6', '6 35 1 1', '35 6 1 1'}
%                               fphiHHCBr=fphiHHCBr+sim(net_fphi_HHCBr,[rij; rik ; ril; rjk; rjl; rkl]);
%                            case {'1 6 6 35', '1 6 35 6', '1 35 6 6', '35 6 6 1', '35 6 1 6', '35 1 6 6', '6 6 1 35'}
%                               fphiHCCBr=fphiHCCBr+sim(net_fphi_HCCBr,[rij; rik ; ril; rjk; rjl; rkl]);  
%                             otherwise
%                                disp(identify4);
%                         end
%                        
%                     end
%                     
%                 end
            end
           
        end
    end
    
%     
%     for i=1:1:clust_size(iQ)
%         for j=1:1:clust_size(iQ)
%                                        
%            
%         end
%     end
    
%     Vhat_iQ=fc*(frHH + frHBr + frHC + frCC + frCBr + ...
%                 fthetaHHH + fthetaHHC + fthetaHHBr + fthetaHCBr + fthetaHCC + fthetaCCBr + ...
%                 fphiHHHC + fphiHHHBr + fphiHHCC + fphiHHCBr + fphiHCCBr);
Vhat_iQ=fc*(frHH + frHBr + frHC + frCC + frCBr + ...
                fthetaHHH + fthetaHHC + fthetaHHBr + fthetaHCBr + fthetaHCC + fthetaCCBr);
   
    Vhat(iQ)=Vhat_iQ; %#ok<AGROW>
    Ex(iQ,1)=V(iQ)-Vhat(iQ);  %#ok<AGROW>
    %SSE
    perf=perf+Ex(iQ,1)*Ex(iQ,1);
     
                 
    
end


i=1;%dummy line

