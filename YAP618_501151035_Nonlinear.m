clc; clear; close all;
dof=3;                          %degree of freedom of nodes
Acc=load('Darfield6953X.txt');
LinSDisp=load('Linear_StoryDisp.txt');
LinBShear=load('Linear_StoryDisp-BaseShear.txt');

%%
%MODEL PROPERTIES
Dynamic=true;                       %SELECT BETWEEN STATIC OR DYNAMIC SOLVER
Example=false;                       %SELECT EXAMPLE FRAME FOR CROSSCHECK  (ONLY IT CHANGES AXES NUMBER AND STATIC LATERAL FORCES)
SpringModel=true;                   %SELECT BETWEEN MODEL WITH SPRING OR NOT
Nonlinear= true;                    %SELECT BETWEEN LINEAR OR NONLINEAR SOLVER (ONLY WITH DYNAMIC SOLVER)
NewtonRaphson=true;                 %SELECT BETWEEN NEWTON RAPHSON OR UNBALANCED FORCE CORRECTION (UFC) METHOD (ONLY WITH DYNAMIC & SPRING & NONLINEAR SOLVER)
EarthQuake=true;                    %SELECT BETWEEN EQ LOAD OR SINUSOIDAL LOAD (ONLY WITH DYNAMIC SOLVER)
LoadShape=1;                        %SELECT STATIC LOAD SHAPE 1=triangular, 2=rectangular (ONLY WITH STATIC SOLVER)

%%
%GRID PROPERTIES
x_axis=2;
y_axis=10;

x_length=7;
y_length=3.5;

if Example
    x_axis=4;
    y_axis=6;
end
%%
%FRAME PROPERTIES
beam_b=0.4;                         %m
beam_h=0.6;                         %m
column_b=0.6;                       %m
column_h=0.6;                       %m
%%
%MATERIAL PROPERTIES
E_Concrete=30*10^6;                 %kN/m2
%%
%DYNAMIC PROPERTIES
ksi=0.05;
%%
%LOAD PROPERTIES
SDead=30;                           %kN/m
load=1*9.806;                       %kN/m^2

%%
%BEAM NONLINEAR PROPERTIES
beamLp=beam_h/2;
beam_k1=5.94*10^4/beamLp;	        %kNm/rad
beam_k2=beam_k1*0.075;            %kNm/rad
beam_My=326.2;                      %kNm
beam_uy=beam_My/beam_k1;
beam_EI=beam_k1;


%COLUMN NONLINEAR PROPERTIES
columnLp=column_h/2;
column_k1=8.35*10^4/columnLp;       %kNm/rad
column_k2=column_k1*0.01;             %kNm/rad
column_My=565.1;                    %kNm
column_uy=column_My/column_k1;
column_EI=column_k1;


%%
%CONSTRAINTS PROPERTIES
for i=1:x_axis*dof
    specdofs(i,1)=i;
end
%%
if SpringModel
    ColumnObj=x_axis*(y_axis-1);    %NUMBER OF COLUMN OBJECT
    BeamObj=(x_axis-1)*(y_axis-1);  %NUMBER OF BEAM OBJECT
    ColumnElm=x_axis*(y_axis-1)*5;  %NUMBER OF COLUMN OBJECT
    BeamElm=(x_axis-1)*(y_axis-1)*5;%NUMBER OF BEAM OBJECT
    SpringElm=2*(ColumnObj+BeamObj);%NUMBER OF SPRING OBJECT
else
    ColumnObj=x_axis*(y_axis-1);
    BeamObj=(x_axis-1)*(y_axis-1);
    ColumnElm=ColumnObj;
    BeamElm=BeamObj;
    SpringElm=0;
end

if SpringModel
    node=x_axis*y_axis+4*ColumnObj+4*BeamObj;
else
    node=x_axis*y_axis;
end
%%
%NODES
if SpringModel
    for i=1:x_axis*y_axis*dof+(2*ColumnObj+2*BeamObj)*dof+(2*ColumnObj+2*BeamObj)
        alldofs(i,1)=i;
    end
else
    for i=1:x_axis*y_axis*dof
        alldofs(i,1)=i;
    end
end
freedofs=alldofs;
freedofs(specdofs)=[];
nofdofs=size(freedofs,1);           %NUMBER OF FREE DEGREE OF FREEDOMS
%%
%COORDINATE MATRIX
if SpringModel
    for i=1:x_axis*y_axis
        coord(i,1)=x_length*(mod((i-1),x_axis));
        coord(i,2)=y_length*floor((i-1)/x_axis);
    end
    i=i+1;
    for ii=1:y_axis-1
        for iii=1:x_axis
            for j=1:2
                coord(i,1)=x_length*(mod((iii-1),x_axis));
                coord(i,2)=y_length*(ii-1)+0.5+columnLp*(j-1);
                i=i+1;
            end
            for j=1:2
                coord(i,1)=x_length*(mod((iii-1),x_axis));
                coord(i,2)=y_length*ii-0.5+columnLp*(j-2);
                i=i+1;
            end
        end
    end
    
    for ii=1:y_axis-1
        for iii=1:x_axis-1
            for j=1:2
                coord(i,1)=x_length*(mod((iii-1),x_axis))+0.5+beamLp*(j-1);
                coord(i,2)=y_length*ii;
                i=i+1;
            end
            for j=1:2
                coord(i,1)=x_length*(mod((iii),x_axis))-0.5+beamLp*(j-2);
                coord(i,2)=y_length*ii;
                i=i+1;
            end
        end
    end
else
    for i=1:node
        coord(i,1)=x_length*(mod((i-1),x_axis));
        coord(i,2)=y_length*floor((i-1)/x_axis);
    end
end
%%
Columnid=[];
Beamid=[];
Springid=[];
Frameid=[];
iii=1;
l=x_axis*y_axis+1;
%CONNECTIVITY MATRIX
if SpringModel
    for i=1:y_axis-1
        for ii=1:x_axis
            connect(iii,1)=(i-1)*x_axis+ii;
            connect(iii,2)=l;
            Columnid=[Columnid;iii];
            iii=iii+1;
            for j=1:3
                connect(iii,1)=l;
                connect(iii,2)=l+1;
                if j==2
                    Columnid=[Columnid;iii];
                else
                    Springid=[Springid;iii];
                end
                iii=iii+1;
                l=l+1;
            end
            connect(iii,1)=l;
            connect(iii,2)=(i-1)*x_axis+ii+x_axis;
            Columnid=[Columnid;iii];
            iii=iii+1;
            l=l+1;
        end
    end
    
    for i=1:y_axis-1
        for ii=1:x_axis-1
            connect(iii,1)=i*x_axis+ii;
            connect(iii,2)=l;
            Beamid=[Beamid;iii];
            iii=iii+1;
            for j=1:3
                connect(iii,1)=l;
                connect(iii,2)=l+1;
                if j==2
                    Beamid=[Beamid;iii];
                else
                    Springid=[Springid;iii];
                end
                iii=iii+1;
                l=l+1;
            end
            connect(iii,1)=l;
            connect(iii,2)=i*x_axis+ii+1;
            Beamid=[Beamid;iii];
            iii=iii+1;
            l=l+1;
        end
    end
    Frameid=[Frameid;Columnid];
    Frameid=[Frameid;Beamid];
else
    for i=1:ColumnObj
        connect(i,1)=i;
        connect(i,2)=i+x_axis;
        Columnid=[Columnid;i];
    end
    
    i=i+1;
    for ii=1:y_axis-1
        for iii=1:x_axis-1
            connect(i,1)=ii*x_axis+iii;
            connect(i,2)=connect(i,1)+1;
            Beamid=[Beamid;i];
            i=i+1;
        end
    end
end
%%
%DEGREE OF FREEDOMS OF NODES MATRIX
if SpringModel
    for i=1:x_axis*y_axis
        dofs(i,1)=1+dof*(i-1);
        dofs(i,2)=2+dof*(i-1);
        dofs(i,3)=3+dof*(i-1);
    end
    l=3+dof*(i-1);
    for i=x_axis*y_axis+1:4:node
        dofs(i,1)=1+l;
        dofs(i,2)=2+l;
        dofs(i,3)=3+l;
        
        dofs(i+1,1)=dofs(i,1);
        dofs(i+1,2)=dofs(i,2);
        dofs(i+1,3)=4+l;
        
        dofs(i+2,1)=5+l;
        dofs(i+2,2)=6+l;
        dofs(i+2,3)=7+l;
        
        dofs(i+3,1)=5+l;
        dofs(i+3,2)=6+l;
        dofs(i+3,3)=8+l;
        l=l+8;
    end
else
    for i=1:node
        dofs(i,1)=1+dof*(i-1);
        dofs(i,2)=2+dof*(i-1);
        dofs(i,3)=3+dof*(i-1);
    end
end
%%
ii=1;
%DEGREE OF FREEDOMS OF FRAMES MATRIX
for i=1:ColumnElm+BeamElm
    framedofs(i,1:dof)=dofs(connect(i,1),:);
    framedofs(i,dof+1:dof+3)=dofs(connect(i,2),:);
    if   framedofs(i,1)==framedofs(i,4)
        if i<=ColumnElm
            SpringLoc(ii,1)=1;
        else
            SpringLoc(ii,1)=2;
        end
        springdofs(ii,1)=framedofs(i,3);
        springdofs(ii,2)=framedofs(i,6);
        ii=ii+1;
    end
end
%%

k{ColumnElm+BeamElm,7}=[];
for i=1:ColumnElm+BeamElm
    if any(i==Columnid) || any(i==Beamid)
        k{i,1}=E_Concrete;                              %E
        
        if i<=ColumnElm
            k{i,2}=(1/12)*column_b*column_h^3;          %I
        else
            k{i,2}=(1/12)*beam_b*beam_h^3;
        end
        
        if i<=ColumnElm
            k{i,3}=column_b*column_h;                   %A
        else
            k{i,3}=beam_b*beam_h;
        end
        
        k{i,4}=sqrt((coord(connect(i,2),2)-coord(connect(i,1),2))^2+(coord(connect(i,2),1)-coord(connect(i,1),1))^2); %L
        
        k{i,5}=stiff2D(k{i,1},k{i,2},k{i,3},k{i,4});    %local stiffness matrix
        
        if (coord(connect(i,2),1)-coord(connect(i,1),1))==0 %DEGREE
            k{i,6}=90;
        else
            k{i,6}=atan((coord(connect(i,2),2)-coord(connect(i,1),2))/(coord(connect(i,2),1)-coord(connect(i,1),1)));
        end
        
        k{i,7}=T_matrix(k{i,6}).'*k{i,5}*T_matrix(k{i,6}); %global stiffness matrix
    elseif any(i==Springid)
        if i<=ColumnElm
            k{i,1}=column_EI;
            k{i,5}=spring2D(k{i,1});
            k{i,7}=k{i,5};
        else
            k{i,1}=beam_EI;
            k{i,5}=spring2D(k{i,1});
            k{i,7}=k{i,5};
        end
    end
end
%%
%GLOBAL STIFFNESS MATRIX
K=zeros(size(alldofs,1));
for i=1:ColumnElm+BeamElm
    K(framedofs(i,:),framedofs(i,:))=K(framedofs(i,:),framedofs(i,:))+k{i,7};
end
%%
%MASS MATRIX
% if Example
%     Mass=[13.77	24.78	24.78	13.77   13.77	24.78	24.78	13.77   13.77	24.78	24.78	13.77   13.77	24.78	24.78	13.77   12.39	23.40	23.40	12.39];
%
%     for i=size(Mass,2)
%        M(
%     end
% else
for i=1:ColumnObj
    connectlin(i,1)=i;
    connectlin(i,2)=i+x_axis;
end
i=i+1;
for ii=1:y_axis-1
    for iii=1:x_axis-1
        connectlin(i,1)=ii*x_axis+iii;
        connectlin(i,2)=connectlin(i,1)+1;
        i=i+1;
    end
end

for i=1:size(alldofs,1)
    M(i,i)=10^-6;
end

% M=zeros(size(alldofs,1));
for i=1:ColumnObj
    for j=1:2
        M(alldofs(dofs(connectlin(i,j),1)),alldofs(dofs(connectlin(i,j),1)))=M(alldofs(dofs(connectlin(i,j),1)),alldofs(dofs(connectlin(i,j),1)))+(column_b*column_h*25)*y_length/2/9.81;
        M(alldofs(dofs(connectlin(i,j),2)),alldofs(dofs(connectlin(i,j),2)))=M(alldofs(dofs(connectlin(i,j),2)),alldofs(dofs(connectlin(i,j),2)))+10^-6;
        M(alldofs(dofs(connectlin(i,j),3)),alldofs(dofs(connectlin(i,j),3)))=M(alldofs(dofs(connectlin(i,j),3)),alldofs(dofs(connectlin(i,j),3)))+10^-6;
    end
end
for i=ColumnObj+1:BeamObj+ColumnObj
    for j=1:2
        M(alldofs(dofs(connectlin(i,j),1)),alldofs(dofs(connectlin(i,j),1)))=M(alldofs(dofs(connectlin(i,j),1)),alldofs(dofs(connectlin(i,j),1)))+(beam_b*beam_h*25+SDead)*x_length/2/9.81;
        M(alldofs(dofs(connectlin(i,j),2)),alldofs(dofs(connectlin(i,j),2)))=M(alldofs(dofs(connectlin(i,j),2)),alldofs(dofs(connectlin(i,j),2)))+10^-6;
        M(alldofs(dofs(connectlin(i,j),3)),alldofs(dofs(connectlin(i,j),3)))=M(alldofs(dofs(connectlin(i,j),3)),alldofs(dofs(connectlin(i,j),3)))+10^-6;
    end
end
% end
%%
%NATURAL FREQUENCIES AND MODE SHAPES
[Evec,Eval]=eig(K(freedofs,freedofs),M(freedofs,freedofs));
for i=1:nofdofs
    omega(i)=sqrt(Eval(i,i));
end
Periyod=2*pi./omega;
Periyod(1,1)

%%
%RAYLEIGH DAMPING
RD_a0=ksi*((2*omega(1)*omega(end))/(omega(1)+omega(end)));
RD_a1=ksi*(2/(omega(1)+omega(end)));
C=RD_a0*M(freedofs,freedofs)+RD_a1*K(freedofs,freedofs);

%%
%CONTROL
for i=1:nofdofs-1
    Control(1,i)=transpose(Evec(:,i))*K(freedofs,freedofs)*Evec(:,i+1);
    Control(2,i)=transpose(Evec(:,i))*M(freedofs,freedofs)*Evec(:,i+1);
    Control(3,i)=transpose(Evec(:,i))*C*Evec(:,i+1);
end
%%
%LOAD MATRIX
if Dynamic
    if EarthQuake
        %%
        %ACCELERATION TIME HISTORY
        
        %Acc=detrend(Acc);      %Base line correlation
        ga=9.8065*Acc(:,2);     %m/s^2                  ground acceleration
        dt=Acc(2,1)-Acc(1,1);   %s                      time interval
        TT=dt*(length(Acc)-1);  %s                      total time
        t=0:dt:TT;              %                       time vector
        N=length(Acc);          %                       total time steps
        
        gv = cumtrapz(ga)*dt;
        gd = cumtrapz(gv)*dt;
        
        F=(-M(freedofs,freedofs)*ones(1,nofdofs)'*ga')';
        
        figure
        subplot(3,1,1)
        plot(t,ga);
        xlabel('time [s]');
        ylabel('Ground Acceleration [m/s^2]');
        %title('D1 X RSN1160 KOCAELI FAT000');
        
        subplot(3,1,2)
        plot(t,gv);
        xlabel('time [s]');
        ylabel('Ground Velocity [m/s]');
        
        subplot(3,1,3)
        plot(t,gd);
        xlabel('time [s]');
        ylabel('Ground Displacement [m]');
        
    else
        %EXTERNAL HARMONIC FORCE
        Td=10;                    %s harmonic force duration
        dt=0.001;                %s time step
        t=0:dt:Td;               % time vector
        N=length(t);             % time steps number
        w=12.5538168;%2*pi;       %rad/s forcing frequency
        po=20;%193;              %N force amplitude
        
        for i=1:N
            F(i,:)=-M(freedofs,freedofs)*ones(1,nofdofs)'*po*sin(w*t(i)); %Period = 2pi/w
        end
        
        figure
        plot(t,F);
        xlabel('time [s]');
        ylabel('Sinusoidal Force [kN]');
        close
        %%
    end
else
    %%
    F=zeros(node*dof,1);
    if Example
        ii=1;
        FEQ=[44.365	79.8378	79.8378 44.365 88.7301	159.6755 159.6755 88.7301 133.0951 239.5133 239.5133 133.0951 177.4602 319.351 319.351	177.4602 214.9864 406.0276 406.0276 214.9864];
        for i=1:(y_axis-1)
            for j=1:x_axis
                F(dofs(i*x_axis+j,1),1)=FEQ(1,ii);
                ii=ii+1;
            end
        end
    else
        for i=1:(y_axis-1)
            for j=1:x_axis
                if LoadShape==1
                    F(dofs(i*x_axis+j,1),1)=i*(load/(y_axis-1))*M(dofs(i*x_axis+j,1),dofs(i*x_axis+j,1),1);
                elseif LoadShape==2
                    F(dofs(i*x_axis+j,1),1)=load*M(dofs(i*x_axis+j,1),dofs(i*x_axis+j,1),1);
                end
            end
        end
    end
end
%%
if Dynamic
    %INTEGRATION PARAMETERS
    alfa=1/2;
    beta=1/4;
    
    %INTEGRATION CONSTANTS
    a1=1/(beta*dt^2);
    a2=1/(beta*dt);
    a3=1/(2*beta);
    a4=alfa/(beta*dt);
    a5=alfa/beta;
    a6=dt*((alfa/(2*beta))-1);
    a7=dt*(1-(alfa/(2*beta)));
    
    %EFFECTIVE STIFFNESS
    if Nonlinear && SpringModel
    else
        kT=K(freedofs,freedofs)+(a1*M(freedofs,freedofs)+a4*C);
    end
    
    %INITIAL CONDITIONS
    u(freedofs,1)=zeros(size(freedofs,1),1);
    ud(freedofs,1)=zeros(size(freedofs,1),1);
    udd(freedofs,1)=inv(M(freedofs,freedofs))*(F(1,:)'-C*ud(freedofs,1)-K(freedofs,freedofs)*u(freedofs,1));
    
    Fs(freedofs,1)=zeros(size(freedofs,1),1);
    du(freedofs,1)=zeros(size(freedofs,1),1);
    springTeta=zeros(SpringElm,1);
    Moment=zeros(SpringElm,1);
    Ehys=zeros(SpringElm,1);
    ue=zeros(SpringElm,1);
    ER=10^-6;
    a=[-1 1];
    jmax=0;
    RRroof(:,1)=zeros(size(freedofs,1),1);
    if Nonlinear && SpringModel
        if NewtonRaphson
            
            %CALCULATIONS FOR EACH TIME STEP
            for i=2:N
                
                %INCREMENTAL EFFECTIVE LOAD AT TIME STEP i
                dP(i,:)=F(i,:)-F(i-1,:);
                dProof(i,:)=dP(i,:)'+(a2*M(freedofs,freedofs)+a5*C)*ud(freedofs,i-1)+(a3*M(freedofs,freedofs)+a6*C)*udd(freedofs,i-1);
                
                %INITIALIZE
                springTeta(:,i)=springTeta(:,i-1);
                MomentStart(:,i)=Moment(:,i-1);
                ueStart(:,i)=ue(:,i-1);
                duTotal(freedofs,i)=zeros(size(freedofs,1),1);
                
                du(freedofs,i)=dProof(i,:)/((a1*M(freedofs,freedofs)+a4*C)+K(freedofs,freedofs));
                Fs(freedofs,i)=zeros(size(freedofs,1),1);
                Rroof(:,i)=K(freedofs,freedofs)*du(freedofs,i);
                
                %DISPLACEMENT AT TIME STEP i+1
                for j=1:25
                    Fs(freedofs,i)=zeros(size(freedofs,1),1);
                    
                    if j>jmax
                        jmax=j
                    end
                    
                    %LINEAR INTERNAL FORCES
                    for ii=1:ColumnElm+BeamElm
                        if any(ii==Columnid) || any(ii==Beamid)
                            fl{ii,i}=k{ii,7}*du(framedofs(ii,:),i);
                            Fs(framedofs(ii,:),i)=Fs(framedofs(ii,:),i)+ fl{ii,i};
                        end
                    end
                    
                    %BILINEAR
                    for ii=1:SpringElm
                        springDTeta(ii,i)=a*du(springdofs(ii,:),i);
                        
                        if SpringLoc(ii,1)==1
                            [Moment(ii,i),ue(ii,i),Ehys(ii,i),k{Springid(ii,1),1}]=bilin(column_k1,column_k2,column_My,springTeta(ii,i),ueStart(ii,i),MomentStart(ii,i),springDTeta(ii,i),Ehys(ii,i-1));
                        elseif SpringLoc(ii,1)==2
                            [Moment(ii,i),ue(ii,i),Ehys(ii,i),k{Springid(ii,1),1}]=bilin(beam_k1,beam_k2,beam_My,springTeta(ii,i),ueStart(ii,i),MomentStart(ii,i),springDTeta(ii,i),Ehys(ii,i-1));
                        end
                        due(ii,i)=ue(ii,i)-ueStart(ii,i);
                        ueStart(ii,i)=ue(ii,i);
                        springTeta(ii,i)=springTeta(ii,i)+springDTeta(ii,i);
                        dMoment(ii,i)=Moment(ii,i)-MomentStart(ii,i);
                        MomentStart(ii,i)=Moment(ii,i);
                        k{Springid(ii,1),5}=spring2D(k{Springid(ii,1),1});
                        k{Springid(ii,1),7}=k{Springid(ii,1),5};
                        
                        %NONLINEAR INTERNAL FORCES
                        Fs(springdofs(ii,:),i)=Fs(springdofs(ii,:),i)+(transpose(a)*dMoment(ii,i));
                    end
                    
                    %GLOBAL STIFFNESS MATRIX
                    K=zeros(size(alldofs,1));
                    for ij=1:ColumnElm+BeamElm
                        K(framedofs(ij,:),framedofs(ij,:))=K(framedofs(ij,:),framedofs(ij,:))+k{ij,7};
                    end
                    
                    %RESIDUAL FORCE
                    Rroof(:,i)=Rroof(:,i)-Fs(freedofs,i);
                    duTotal(freedofs,i)=duTotal(freedofs,i)+du(freedofs,i);
                    
                    %CHECKING CONVERGENCE
                    if abs(Rroof(:,i))<ER
                        break;
                    else
                        du(freedofs,i)=(Rroof(:,i)'/((a1*M(freedofs,freedofs)+a4*C)+K(freedofs,freedofs)))';
                        Rroof(:,i)=K(freedofs,freedofs)*du(freedofs,i);
                    end
                end
                
                %INCREMENTAL VELOCITY AT TIME STEP i
                dud(freedofs,i)=a4*duTotal(freedofs,i)-a5*ud(freedofs,i-1)-a7*udd(freedofs,i-1);
                
                %INCREMENTAL ACCELERATION AT TIME STEP i
                dudd(freedofs,i)=a1*duTotal(freedofs,i)-a2*ud(freedofs,i-1)-a3*udd(freedofs,i-1);
                
                %STATE OF MOTION AT TIME STEP i+1
                u(freedofs,i)=u(freedofs,i-1)+duTotal(freedofs,i);
                ud(freedofs,i)=ud(freedofs,i-1)+dud(freedofs,i);
                udd(freedofs,i)=udd(freedofs,i-1)+dudd(freedofs,i);
            end
            %%
            
        else
            %CALCULATIONS FOR EACH TIME STEP
            for i=2:N
                
                %INCREMENTAL EFFECTIVE LOAD AT TIME STEP i
                dP(i,:)=F(i,:)-F(i-1,:);
                dProof(i,:)=dP(i,:)'+(a2*M(freedofs,freedofs)+a5*C)*ud(freedofs,i-1)+(a3*M(freedofs,freedofs)+a6*C)*udd(freedofs,i-1);
                
                %INITIALIZE
                du(freedofs,i)=dProof(i,:)/((a1*M(freedofs,freedofs)+a4*C)+K(freedofs,freedofs));
                Fs(freedofs,i)=zeros(size(freedofs,1),1);
                Rroof(:,i)=K(freedofs,freedofs)*du(freedofs,i)+RRroof(:,i-1);
                
                %LINEAR INTERNAL FORCES
                for ii=1:ColumnElm+BeamElm
                    if any(ii==Columnid) || any(ii==Beamid)
                        fl{ii,i}=k{ii,7}*du(framedofs(ii,:),i);
                        Fs(framedofs(ii,:),i)=Fs(framedofs(ii,:),i)+ fl{ii,i};
                    end
                end
                
                %BILINEAR
                for ii=1:SpringElm
                    springDTeta(ii,i)=a*du(springdofs(ii,:),i);
                    
                    if SpringLoc(ii,1)==1
                        [Moment(ii,i),ue(ii,i),Ehys(ii,i),k{Springid(ii,1),1}]=bilin(column_k1,column_k2,column_My,springTeta(ii,i-1),ue(ii,i-1),Moment(ii,i-1),springDTeta(ii,i),Ehys(ii,i-1));
                    elseif SpringLoc(ii,1)==2
                        [Moment(ii,i),ue(ii,i),Ehys(ii,i),k{Springid(ii,1),1}]=bilin(beam_k1,beam_k2,beam_My,springTeta(ii,i-1),ue(ii,i-1),Moment(ii,i-1),springDTeta(ii,i),Ehys(ii,i-1));
                    end
                    springTeta(ii,i)=springTeta(ii,i-1)+springDTeta(ii,i);
                    dMoment(ii,i)=Moment(ii,i)-Moment(ii,i-1);
                    k{Springid(ii,1),5}=spring2D(k{Springid(ii,1),1});
                    k{Springid(ii,1),7}=k{Springid(ii,1),5};
                    
                    %NONLINEAR INTERNAL FORCES
                    Fs(springdofs(ii,:),i)=Fs(springdofs(ii,:),i)+(transpose(a)*dMoment(ii,i));
                end
                
                %GLOBAL STIFFNESS MATRIX
                K=zeros(size(alldofs,1));
                for ij=1:ColumnElm+BeamElm
                    K(framedofs(ij,:),framedofs(ij,:))=K(framedofs(ij,:),framedofs(ij,:))+k{ij,7};
                end
                
                %RESIDUAL FORCE
                RRroof(:,i)=Rroof(:,i)-Fs(freedofs,i);
                
                %INCREMENTAL VELOCITY AT TIME STEP i
                dud(freedofs,i)=a4*du(freedofs,i)-a5*ud(freedofs,i-1)-a7*udd(freedofs,i-1);
                
                %INCREMENTAL ACCELERATION AT TIME STEP i
                dudd(freedofs,i)=a1*du(freedofs,i)-a2*ud(freedofs,i-1)-a3*udd(freedofs,i-1);
                
                %STATE OF MOTION AT TIME STEP i+1
                u(freedofs,i)=u(freedofs,i-1)+du(freedofs,i);
                ud(freedofs,i)=ud(freedofs,i-1)+dud(freedofs,i);
                udd(freedofs,i)=udd(freedofs,i-1)+dudd(freedofs,i);
            end
            %%
        end
    else
        %%
        for i=2:N
            dP(i,:)=F(i,:)-F(i-1,:);
            Proof(i,:)=dP(i,:)'+(a2*M(freedofs,freedofs)+a5*C)*ud(freedofs,i-1)+(a3*M(freedofs,freedofs)+a6*C)*udd(freedofs,i-1);
            
            %SOLUTION FOR du AT TIME STEP i
            du(freedofs,i)=Proof(i,:)/kT;
            
            %INCREMENTAL VELOCITY AT TIME STEP i
            dud(freedofs,i)=a4*du(freedofs,i)-a5*ud(freedofs,i-1)-a7*udd(freedofs,i-1);
            
            %INCREMENTAL ACCELERATION AT TIME STEP i
            dudd(freedofs,i)=a1*du(freedofs,i)-a2*ud(freedofs,i-1)-a3*udd(freedofs,i-1);
            
            %STATE OF MOTION AT TIME STEP i+1
            u(freedofs,i)=u(freedofs,i-1)+du(freedofs,i);
            ud(freedofs,i)=ud(freedofs,i-1)+dud(freedofs,i);
            udd(freedofs,i)=udd(freedofs,i-1)+dudd(freedofs,i);
        end
        %%
    end
else
    %%
    u=zeros(size(alldofs,1),1);
    u(freedofs)=inv(K(freedofs,freedofs))*(F(freedofs)-K(freedofs,specdofs)*u(specdofs));
    f(specdofs)=K(specdofs,:)*u;
    
end
%%

if Dynamic
    Base(specdofs,:)=K(specdofs,:)*u;
    BaseShear(1,:)=sum(Base(1:3:end,:));
end

%%
%%%%%---PLOTS---%%%%%

%Scale factor for deformed shape display
if Dynamic
    sf=50;%5/(max(abs(u(:,x_axis*(y_axis-1)*dof+1)))+coord(x_axis*(y_axis-1)+1,1));
else
    sf=50;
end

if Dynamic && Nonlinear
    %MOMENT-ROTATION
    for i=1:2:SpringElm
        figure
        subplot(1,2,1)
        plot(springTeta(i,:),Moment(i,:))
        title(sprintf('Frame-%d i End Moment-Rotation',(i+1)/2));
        xlabel('Rotation [rad]');
        ylabel('Moment [kNm]');
        
        subplot(1,2,2)
        plot(springTeta(i+1,:),Moment(i+1,:))
        title(sprintf('Frame-%d j End Moment-Rotation',(i+1)/2));
        xlabel('Rotation [rad]');
        ylabel('Moment [kNm]');
    end
    
    %     %         Ed(1)=0;
    %     %         for i=1:N-1
    %     %             Ed(i+1)=Ed(i)+(abs((fd(i+1)+fd(i))/2)*abs(u(i+1)-u(i)));
    %     %         end
    %
    %     EhysTotal=sum(Ehys,2);
    %     figure
    %     hold on
    %     %area(t,Ek);
    %     area(t,EhysTotal','FaceColor',[0.5 0.9 0.6],'EdgeColor',[0 0.5 0.1]);
    %     %         area(t,Ed,'FaceColor',[0.7 0.7 0.7],'EdgeColor','k');
    %     legend('Hysteresis');%,'Damping');%,'Kinetic');
    %     hold off
end

if Dynamic
    %ROOF DISPLACEMENT & BASE SHEAR
    figure
    plot(u(x_axis*y_axis*dof-2,:),BaseShear(1,:))
    title('Roof Displacement & Base Shear');
    xlabel('Roof Displacement [m]');
    ylabel('Base Shear [kN]');
    
    if Nonlinear && SpringModel
    else
        hold on
        plot(LinBShear(:,1),LinBShear(:,2),'R');
        legend('Matlab','Sap2000')
        hold off
    end
    
    %STORY DISPLACEMENTS
    figure
    hold on
    plot(t(1,:),u(x_axis*y_axis*dof-2,:))
    hold off
    title('Story Displacements');
    xlabel('Time [s]');
    ylabel('Displacement [m]');
    
    if Nonlinear && SpringModel
    else
        hold on
        plot(LinSDisp(:,1),LinSDisp(:,2),'R');
        legend('Matlab','Sap2000')
        hold off
    end
end

figure
hold on
title('Displacement');
axis([-5,(x_axis-1)*x_length+5,0,(y_axis-1)*y_length+1])

for i=1:(BeamElm+ColumnElm)
    plot(coord(connect(i,:),1),coord(connect(i,:),2),'k--')
end

if Dynamic
    
    h= animatedline('Color','r','LineWidth',2);
    for j = 1:N
        
        ndivs=1;
        for i=1:(BeamElm+ColumnElm)
            elnodes=connect(i,1:2);
            E1=[(coord(elnodes(2),1)-coord(elnodes(1),1))...
                (coord(elnodes(2),2)-coord(elnodes(1),2))];
            le=norm(E1);
            E1=E1/le;
            E2=[-E1(2) E1(1)];
            
            eldofs=framedofs(i,:);
            eldisp=u(eldofs,j);
            
            Qrot=[E1;E2];   %Transforms global to element d_E = Q d_G
            Qrot(dof,dof)=1;
            Tmatrix=[Qrot zeros(dof); zeros(dof) Qrot];
            eldispLOC=Tmatrix*eldisp;
            
            for ii=1:ndivs+1
                xi=(ii-1)/ndivs;
                xdispLOC=eldispLOC(1)*(1-xi)+eldispLOC(4)*xi;
                ydispLOC=eldispLOC(2)*(1-3*xi^2+2*xi^3)+eldispLOC(5)*(3*xi^2-2*xi^3)...
                    +eldispLOC(3)*le*(xi-2*xi^2+xi^3)+eldispLOC(6)*le*(-xi^2+xi^3);
                
                xydisp=(Qrot([1,2],[1,2]))'*[xdispLOC;ydispLOC];
                plotpts(ii,1)=coord(elnodes(1),1)+xi*le*E1(1)+sf*xydisp(1);
                plotpts(ii,2)=coord(elnodes(1),2)+xi*le*E1(2)+sf*xydisp(2);
            end
            %h(i) = animatedline(plotpts(:,1),plotpts(:,2),'Color','r','LineWidth',2);
            addpoints(h,plotpts(:,1),plotpts(:,2))
        end
        drawnow
        %delete(h(:,:))
        %pause(dt/1)
        clearpoints(h)
    end
else
    ndivs=50;
    for i=1:(BeamElm+ColumnElm)
        elnodes=connect(i,1:2);
        E1=[(coord(elnodes(2),1)-coord(elnodes(1),1))...
            (coord(elnodes(2),2)-coord(elnodes(1),2))];
        le=norm(E1);
        E1=E1/le;
        E2=[-E1(2) E1(1)];
        
        eldofs=framedofs(i,:);
        eldisp=u(eldofs);
        
        Qrot=[E1;E2];   %Transforms global to element d_E = Q d_G
        Qrot(dof,dof)=1;
        Tmatrix=[Qrot zeros(dof); zeros(dof) Qrot];
        eldispLOC=Tmatrix*eldisp;
        
        for ii=1:ndivs+1
            xi=(ii-1)/ndivs;
            xdispLOC=eldispLOC(1)*(1-xi)+eldispLOC(4)*xi;
            ydispLOC=eldispLOC(2)*(1-3*xi^2+2*xi^3)+eldispLOC(5)*(3*xi^2-2*xi^3)...
                +eldispLOC(3)*le*(xi-2*xi^2+xi^3)+eldispLOC(6)*le*(-xi^2+xi^3);
            
            xydisp=(Qrot([1,2],[1,2]))'*[xdispLOC;ydispLOC];
            plotpts(ii,1)=coord(elnodes(1),1)+xi*le*E1(1)+sf*xydisp(1);
            plotpts(ii,2)=coord(elnodes(1),2)+xi*le*E1(2)+sf*xydisp(2);
        end
        plot(plotpts(:,1),plotpts(:,2),'r.-','LineWidth',2)
    end
end

