                                                                                                                                                                                                                                                               %Victor Adrian Sosa Hernandez
%ITESM-CEM
%vsosa@tec.mx

%Multiple stepsize for solving a given MOP using the gradient hypervolume
%(GH) and the approximation of the Hessian (HH)
%
%Notation: 
%    k: number objectives 
%    n: number variables

%INPUT:
%  opc: Select the test problem
% iter: Number of iterations of the algorithm
% khes: (1) to use hessian approximation and (2) to use real hessian
%  out: name of an archive to store results (this parameter can be empty)

%Example
%multiple_step_size_test_ultimate(1,6,1,2)


function []=Newton_methodEC(opc,iter,out)
    warning off;
    close all;
    %tic;
    global FH;  %Multiobjective optimization problem
    global DFH;
    global DDFH;
    global HH;
    global DHH;
    global DDHH;
    global R;   %Reference point
    global n;   %Number of variables
    global mu;  %Number of points in the archive
    global k;   %Number of objectives
    global nus; %Directions in objective space
   
    if opc==1
        xi=load('Test_5_Convex2.txt');
        %xi=[ -0.9845    0.1538;
   % 0.9542   -0.3047;
   %-0.2449   -0.9741;
   %-0.2643    0.9688;
   % 0.6805    0.7316]
        %idealxi=load('Ideal_Convex2_5.txt');
        
        R=[20 20];
        n=2;
        mu=5;
        k=2;
        [FH,DFH,DDFH,HH,DHH,DDHH]=Convex2GHH();
        ps=load('Convex2GHH_RPS.txt');
        pf=load('Convex2GHH_RPF.txt');  
        hvmax= 373.1665222413705000 ;  
        lambda=ones(1,mu)*1/mu;
        refPS=load('reference_Convex2GHH_PS.txt');
        refPF=load('reference_Convex2GHH_PF.txt');
    elseif opc==2
        xi=load('Test_5_Convex2_2.txt');
        %idealxi=load('Ideal_Convex2_5.txt');
        R=[20 20];
        n=2;
        mu=3;
        k=2;
        [FH,DFH,DDFH,HH,DHH,DDHH]=Convex3GHH();
        ps=load('Convex3GHH_RPS.txt');
        pf=load('Convex2GHH_RPF.txt');  
        hvmax= 373.1665222413705000 ;  
        lambda=ones(1,mu)*1/mu;
        refPS=load('reference_Convex2GHH_PS.txt');
    elseif opc==3
        xi=load('Test_2_ZDT1MH.txt');
        %idealxi=load('Ideal_Convex2_5.txt');
        R=[11 11];
        n=2;
        mu=2;
        k=2;
        [FH,DFH,DDFH,HH,DHH,DDHH]=ZDT1MGHH(n);
        ps=load('ZDT1M-0.5_0_-0.25_PS.txt');
        pf=load('ZDT1M-0.5_0_-0.25_PF.txt');  
        hvmax= 373.1665222413705000;  
        lambda=ones(1,mu)*1/mu;
        refPS=load('reference2_ZDT1_-0.5_0_-0.25.txt');
    elseif opc==4
        xi=load('Test_5_ZDT1MH.txt');
        %idealxi=load('Ideal_Convex2_5.txt');
        R=[11 11];
        n=2;
        mu=3;
        k=2;
        [FH,DFH,DDFH,HH,DHH,DDHH]=ZDT1MGHH(n);
        ps=load('ZDT1M-0.5_0_-0.25_PS.txt');
        pf=load('ZDT1M-0.5_0_-0.25_PF.txt'); 
        hvmax= 373.1665222413705000;  
        lambda=ones(1,mu)*1/mu;
        refPS=load('reference3_ZDT1_-0.5_0_-0.25.txt');
    elseif opc==5
         %xi=load('Test_5_Convex2.txt');
        xi=[ 0.9853   -0.1589;
   -0.9932   -0.0907;
   -0.9178    0.4074;
   -0.4787    0.8821;
   -0.7400   -0.6776;
   -0.3326   -0.9417;
    0.8790    0.4709;
    0.1667   -0.9821;
    0.0518    0.9946;
    0.5121    0.8594];
        %idealxi=load('Ideal_Convex2_5.txt');
        
        R=[20 20];
        n=2;
        mu=10;
        k=2;
        [FH,DFH,DDFH,HH,DHH,DDHH]=Convex2GHH();
        ps=load('Convex2GHH_RPS.txt');
        pf=load('Convex2GHH_RPF.txt');  
        hvmax= 373.1665222413705000 ;  
        lambda=ones(1,mu)*1/mu;
        refPS=load('reference_Convex2GHH_PS.txt');
        refPF=load('reference_Convex2GHH_PF.txt');
    end
    
    
       
    
    %Defining pictures configuration
    figure(1)
    hold on;
    axis equal;
    
    xlabel(['$x_1$'],'interpreter','latex', 'FontSize', 22)
    ylabel(['$x_2$'],'interpreter','latex', 'FontSize', 22)
   
    
    
    if(opc==6)
        zlabel(['$x_3$'],'interpreter','latex', 'FontSize', 22)
        plot3(ps(:,1),ps(:,2),ps(:,3),'.k','LineWidth',2)
    else
        plot(ps(:,1),ps(:,2),'-k','LineWidth',2)        
    end
   
    figure(2)
    hold on;
    axis equal;
    
    xlabel(['$f_1$'],'interpreter','latex', 'FontSize', 22)
    ylabel(['$f_2$'],'interpreter','latex', 'FontSize', 22)
    plot(pf(:,1),pf(:,2),'-k','LineWidth',2)
    
    %Establishing image per setting
    fig_PS=1;
    fig_OS=2;
    
    size_mark=9;
    
    %Resize archive according the multiset representation
    multisets_x=reshape(xi',1,n*mu);
    %Plotting multisets in parameter space
    
    if(opc==6)        
        plot_multiset3D(multisets_x,fig_PS,'r','o',size_mark+2);
    else        
        plot_multiset2D(multisets_x,fig_PS,'r','o',size_mark);
    end
    
    
    %Evaluating multisets
    multisets_fx=evaluate_multiset(multisets_x);
    
    
    %Ordering the first multiset
    xaux=vec2mat(multisets_x,length(multisets_x)/mu);
    faux=vec2mat(multisets_fx,length(multisets_fx)/mu);
    [faux,ind]=sortrows(faux);
    xaux=xaux(ind,:);
    multisets_x=reshape(xaux',1,n*mu);
    multisets_fx=reshape(faux',1,k*mu);
    
    %Plotting multisets in objective space
    plot_multiset2D(multisets_fx,fig_OS,'r','o',size_mark);
    
    if nargin <5
        fout=1;
    else
        fout = fopen(out,'wt');
    end
    
    hv=evaluate_hypervolume(multisets_fx);    
    
    [dh,LH] = gradient_hypervolume(multisets_x,multisets_fx,0); 
    xx=vec2mat(multisets_x,length(multisets_x)/mu);
    %xfx=vec2mat(multisets_fx,length(multisets_fx)/mu);
    
    dp=distAB(refPS,xx);
    res=computenormH(multisets_x);
    error(1,:)=[0 ,hv, hvmax, dp, res];
   
    fprintf(fout,'%d & \t %5.10f & \t %5.10f \\\\ \n',error(1,1),error(1,4),error(1,5));
    
    ddh=hessian_hypervolume(multisets_x,LH);
   
    xslambda=[multisets_x,lambda];
    for i=1:iter
        H=jacobian_equality(multisets_x);
        G=computeG(multisets_x,dh,H,lambda);
        DG=computeDG(multisets_x,ddh,H,lambda);
        nus=(DG\G)';
        
        auxmultisets_x=multisets_x;
        auxmultisets_fx=multisets_fx;
        
        xslambda=xslambda-nus;
        multisets_x=xslambda(1:mu*n);
        multisets_fx=evaluate_multiset(multisets_x);        
        lambda=xslambda(mu*n+1:end);
        hv=evaluate_hypervolume(multisets_fx);
        [dh,LH] = gradient_hypervolume(multisets_x,multisets_fx,0);    
        ddh=hessian_hypervolume(multisets_x,LH);
        xx=vec2mat(multisets_x,length(multisets_x)/mu);
        %xfx=vec2mat(multisets_fx,length(multisets_fx)/mu);
        dp=distAB(refPS,xx);
        res=computenormH(multisets_x);
        error(i+1,:)=[i ,hv, hvmax, dp,res];
        fprintf(fout,'%d & \t %5.10f & \t %5.10f \\\\ \n',error(i+1,1),error(i+1,4),error(i+1,5));
    
        if(opc==6)        
            plot_multiset3D(multisets_x,fig_PS,'y','^',size_mark);
            quiver_multiset3(auxmultisets_x,multisets_x,fig_PS,'r',3)
        else        
            plot_multiset2D(multisets_x,fig_PS,'y','^',size_mark);
            quiver_multiset(auxmultisets_x,multisets_x,fig_PS,'r',3)
        end
        plot_multiset2D(multisets_fx,fig_OS,'y','^',size_mark); 
        quiver_multiset(auxmultisets_fx,multisets_fx,fig_OS,'r',3)
    end
    last1=vec2mat(multisets_x,length(multisets_x)/mu);
    dlmwrite('reference1.txt', last1,'delimiter',' ', 'precision', '%1.10f');
    
    figure(1)
    set(gca, 'FontSize', 16)
    if(opc==6)
        view(-6,20)
    end
    figure(2)
    set(gca, 'FontSize', 16)
    figure(3)
    hold on;
    axis equal;
    
    xlabel(['$f_1$'],'interpreter','latex', 'FontSize', 22)
    ylabel(['$f_2$'],'interpreter','latex', 'FontSize', 22)
    plot(pf(:,1),pf(:,2),'-k','LineWidth',2)
    plot_multiset2D(multisets_fx,3,'y','^',size_mark+2); 
    set(gca, 'FontSize', 16)
    
    figure(4)
    hold on;
    axis equal;
    
    xlabel(['$x_1$'],'interpreter','latex', 'FontSize', 22)
    ylabel(['$x_2$'],'interpreter','latex', 'FontSize', 22)
    
    if(opc==6)
        zlabel(['$x_3$'],'interpreter','latex', 'FontSize', 22)
        plot3(ps(:,1),ps(:,2),ps(:,3),'.k','LineWidth',2)
        plot_multiset3D(multisets_x,4,'y','^',size_mark+2);
        view(-6,20)
    else
        plot(ps(:,1),ps(:,2),'-k','LineWidth',2)
        plot_multiset2D(multisets_x,4,'y','^',size_mark+2); 
    end
    set(gca, 'FontSize', 16)
    %time=toc;
end

function [res]=computenormH(multisets_x)
    global HH;
    global mu;    
    xaux=vec2mat(multisets_x,length(multisets_x)/mu);
    nxs=[];
    for i=1:mu
        nxs=[nxs, HH(xaux(i,:)')];     
    end
    res=norm(nxs);
end

function [G]=computeG(multisets_x,dh,H,lam)
    global HH;
    global mu;    
    xaux=vec2mat(multisets_x,length(multisets_x)/mu);
    nxs=dh;
    for i=1:mu
        nxs=nxs+H(i,:)*lam(i);                
    end
    for i=1:mu
        nxs=[nxs, HH(xaux(i,:)')];     
    end
    G=nxs';
end


function DG=computeDG(multisets_x,ddh,H,lam)
    
    global mu;    
    nxs=ddh;
    for i=1:mu
        nxs=nxs+hessian_equality(multisets_x,i) *lam(i);                
    end
    
    DG=[nxs,H'; H zeros(mu,mu)];
end

%Function to compute the  multiple step sizes

function [y]=problemG2(xslam)
    global HH;
    global mu;
    global n;
    global k;
    xs=xslam(1:mu*n);
    lam=xslam(mu*n+1:end);    
    fxs=evaluate_multiset(xs);
    xaux=vec2mat(xs,length(xs)/mu);
    faux=vec2mat(fxs,length(fxs)/mu);
    [faux,ind]=sortrows(faux);
    xaux=xaux(ind,:);
    xs=reshape(xaux',1,n*mu);
    fxs=reshape(faux',1,k*mu);
    lam=lam(ind);
    [nxs,~] = gradient_hypervolume(xs,fxs,0);
    matrix_H=jacobian_equality(xs);
    for i=1:mu
        nxs=nxs+matrix_H(i,:)*lam(i);                
    end
    for i=1:mu
        nxs=[nxs, HH(xaux(i,:)')];
    end
   
    y=norm(nxs);
end


%This function draw the arrows from multiset to multiset in parameter and
%objective space
function quiver_multiset(auxmultiset,multiset,num_fig,colorline,siz)
    global mu;
    newmultiset=auxmultiset-multiset;
    newmultiset=vec2mat(newmultiset,length(newmultiset)/mu);
    auxmultiset=vec2mat(auxmultiset,length(auxmultiset)/mu);
    
    
    figure(num_fig);    
    for i=1:mu
       h= quiver(auxmultiset(i,1),auxmultiset(i,2),-newmultiset(i,1),-newmultiset(i,2),'color',colorline);
       %if nargin ==5
       %    adjust_quiver_arrowhead_size(h, siz);
       %end
    end
        
end

function quiver_multiset3(auxmultiset,multiset,num_fig,colorline,siz)
    global mu;
    newmultiset=auxmultiset-multiset;
    newmultiset=vec2mat(newmultiset,length(newmultiset)/mu);
    auxmultiset=vec2mat(auxmultiset,length(auxmultiset)/mu);
    
    
    figure(num_fig);    
    for i=1:mu
       h= quiver3(auxmultiset(i,1),auxmultiset(i,2),auxmultiset(i,3),-newmultiset(i,1),-newmultiset(i,2),-newmultiset(i,3),'color',colorline);
       %if nargin ==5
       %    adjust_quiver_arrowhead_size(h, siz);
       %end
    end
        
end

%This function plot a multiset over a given space only for 2 objectives
function plot_multiset2D(multiset,num_fig,color,type_fig,size_fig)    
    global mu;
    multiset=vec2mat(multiset,length(multiset)/mu);
    figure(num_fig);
    if nargin < 2
        plot(multiset(:,1),multiset(:,2),'k.');
    else
        plot(multiset(:,1),multiset(:,2),type_fig,'MarkerFaceColor',color,'MarkerEdgeColor','k','MarkerSize',size_fig);
    end    
end

function plot_multiset3D(multiset,num_fig,color,type_fig,size_fig)    
    global mu;
    multiset=vec2mat(multiset,length(multiset)/mu);
    figure(num_fig);
    if nargin < 2
        plot3(multiset(:,1),multiset(:,2),multiset(:,3),'k.');
    else
        plot3(multiset(:,1),multiset(:,2),multiset(:,3),type_fig,'MarkerFaceColor',color,'MarkerEdgeColor','k','MarkerSize',size_fig);
    end    
end


%This function evaluates a multiset from R^(mu*n)->R^(mu*k)
function [multiset_fx]=evaluate_multiset(multiset_x)
    global mu;
    global k;
    global FH;
    aux=vec2mat(multiset_x,length(multiset_x)/mu);
    for i=1:mu
        faux(i,:)=FH(aux(i,:)')';
    end
    multiset_fx=reshape(faux',1,k*mu);
end


%This function compute HG and sort the population 
%INPUT:
%    multiset_x: a vector that contains all the decision variables of each
%    element in the population R^(mu*n)
%   multiset_fx: a vector that contains the image of all the elements in the
%    population in only one row R^(mu*k)
%     flag_norm: if this flag is 1 the HG is divided each element by each
%     2-norm
%OUTPUT:
%            dh: HG is a vector R^(mu*n)
%         new_x: return the population sorted in decision variable space R^(mu*n)
%        new_fx: return the population sorted in objective space R^(mu*k)
%            LH: return internal reference points information

function [dh,LH] = gradient_hypervolume(multiset_x,multiset_fx,flag_norm)    
    global R;  %Reference point
    global mu; %Number of elements in the population
    global DFH;%Function handle of the Jacobian of the MOP 
    global k;  %Number of objectives
    global n;  %Number of decision variables
       
    %Build the new reference points
    xaux=vec2mat(multiset_x,length(multiset_x)/mu);
    faux=vec2mat(multiset_fx,length(multiset_fx)/mu);    
    %Compute the HG
    dh=zeros(1,mu*n);
    LH=zeros(mu,k);
    for i=1:mu
        index1=i+1;
        index2=i-1;
        if i==1
            if mu==1
                r=[R(1), R(2)];
            else
                r=[faux(index1,1), R(2)];
            end
        elseif i==mu
            r=[R(1),faux(index2,2)];            
        else
            r=[faux(index1,1) ,faux(index2,2)]; 
        end
        dx=DFH(xaux(i,:)');
                
        LH(i,:)=fliplr(-(r-faux(i,:)));
        
        
        aux=(i-1)*n;
        for c1=1:n
           for c2=1:k
               dh(aux+c1)=dh(aux+c1)+LH(i,c2)*dx(c2,c1);
           end
        end       
        
    end    
    if flag_norm==1
        dh = dh/norm(dh,2);
    end
end

function [dde] = hessian_equality(multiset_x,var)    
    global mu; %Number of elements in the population
    global DDHH;%Function handle of the gradient of the equality 
    global n;  %Number of decision variables
       
    %Build the new reference points
    xaux=vec2mat(multiset_x,length(multiset_x)/mu);    
    %Compute the HG
    dde=zeros(mu*n,mu*n);
    con=1;
    for i=1:mu
        aux_con=i*n;
        if i==var
            dde(con:aux_con,con:aux_con)=DDHH(xaux(i,:)');
        end
        con=aux_con+1;
    end 
       
end

function [de] = gradient_equality(multiset_x,var)    
    global mu; %Number of elements in the population
    global DHH;%Function handle of the gradient of the equality 
    global n;  %Number of decision variables
       
    %Build the new reference points
    xaux=vec2mat(multiset_x,length(multiset_x)/mu);    
    %Compute the HG
    de=zeros(mu,n);
    for i=1:mu
        if i==var
            de(i,:)=DHH(xaux(i,:)');
        end
    end 
    de=reshape(de',1,mu*n);    
end

function [H] = jacobian_equality(multiset_x)    
    global mu; %Number of elements in the population
    global n;  %Number of decision variables
    
    H=zeros(mu,mu*n);
    for i=1:mu
        H(i,:)=gradient_equality(multiset_x,i);
    end 
    
end

%Compute the hypervolume of a multiset
function [multiset_h]=evaluate_hypervolume(multiset_fx)
    global mu;
    global R;
    faux=vec2mat(multiset_fx,length(multiset_fx)/mu)';
    multiset_h=lebesgue_measure(faux,R');    
end




%Compute the Hypervolume Hessian using the mathematical formulation
%INPUT:
%    multiset_x: a vector that contains all the decision variables of each
%    LH: return internal reference points information
%OUTPUT:
%         ddh: HH is a matrix R^(mu*n)x(mu*n)

function ddh=hessian_hypervolume(multiset_x,LH)
    global DFH;%Function handle of the Jacobian of the MOP 
    global DDFH;%Function handle of the Hessian of the MOP 
    global k;
    global mu; %Number of elements in the population
    global n;
    xaux=vec2mat(multiset_x,length(multiset_x)/mu);
    lim=mu*n;
    ddh=zeros(lim,lim);
    
    lim=mu*n;
    for a=1:lim
        i=ceil(a/n);
        j=mod(a,n);
        if j==0
            j=n;
        end
        dx=DFH(xaux(i,:)');
        ddx=DDFH(xaux(i,:)');
        for b=1:lim
            l=ceil(b/n);
            p=mod(b,n);
            if p==0
                p=n;
            end
            for z=1:k
                %Compute component 2
                if l==i
                    ddh(a,b)=ddh(a,b)+LH(i,z)*ddx(k*(j-1)+z,p);
                else
                    ddh(a,b)=ddh(a,b)+0;
                end
                
                
                %Compute component 1
                if l==i && z==1
                    ddh(a,b)=ddh(a,b)+dx(z,j)*dx(2,p);
                elseif l==(i-1) && z==1
                    dxau=DFH(xaux(i-1,:)');
                    ddh(a,b)=ddh(a,b)-dx(z,j)*dxau(2,p);
                elseif l==i && z==2
                    ddh(a,b)=ddh(a,b)+dx(z,j)*dx(1,p);
                elseif l==(i+1) && z==2
                    dxau=DFH(xaux(i+1,:)');
                    ddh(a,b)=ddh(a,b)-dx(z,j)*dxau(1,p);
                else
                    ddh(a,b)=ddh(a,b)+0;
                end
                
            end  
            %fprintf(1,'i=%d, j=%d, l=%d, p=%d , a=%d, b=%d \n',i,j,l,p, a, b);
        end
    end
    
end