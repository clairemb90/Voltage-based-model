
function S=model_neocortex1(y,set)

dt=0.1;
Ed=0;
theta_moins_ltp=y(3);
theta_moins_ltd=y(4);
[trace,filetoread,raw]=xlsread('traces_neocortex.xlsx',set);

if set==1 %Letzkus
    data=9;
    inter_dt=1000; %interval %ms
    nb_pairs(1:data)=150;
elseif set==2 %Sjostrom
    data=12;
    nb_pairs=[50;50;50;15;15;15;15;50;15;15;15;15];
    intervals=[0;10;25;10;0;10;25;-10;-10;-10;-10;10]; % %time delay (ms)
    freq=[0.1;0.1;0.1;10;20;20;20;0.1;10;20;40;40]; %frequency (Hz)
    interval=(1./freq)*1000; %ms
    inter_dt=floor(interval./dt);
end

dw=zeros(data,1);
norm_weight=zeros(data,1);
ww=0.5*ones(data,1);

for expdata=1:data
    
    if set==1
        u(50:50+length(trace)-6)=(trace(6:length(trace),expdata)-(trace(6,expdata)))*1000;
        u(isnan(u)) = 0;
        titi=zeros(length(u),1);
        if  expdata==1 || expdata==3 || expdata==5 || expdata==7 || expdata==9   
            titi(50)=1; %pre-post
        elseif expdata==2 || expdata==4 || expdata==6 || expdata==8 
            titi(50+10/dt)=1; %post-pre 10 ms
        end
    elseif set==2
        v(50:50+length(trace(:,expdata))-1)=trace(1:length(trace),expdata)-(trace(1,expdata));
        u=v(~isnan(v));
        titi=zeros(length(u),1);
        if expdata==8
            titi(50+10/dt)=1;
        elseif expdata==9
            for ll=1:5
                titi(50+10/dt+(ll-1)*100/dt)=1;
            end
        elseif expdata==10
            for ll=1:5
                titi(50+10/dt+(ll-1)*50/dt)=1;
            end
        elseif expdata==4
            for ll=1:5
                titi(50+(ll-1)*100/dt)=1;
            end
        elseif expdata==5 || expdata==6 || expdata==7
            for ll=1:5
                titi(50+(ll-1)*50/dt)=1;
            end
        else
            titi(50)=1;
        end
        
    end
    
    toto=zeros(length(u),1);
    wi=zeros(length(u),1); %weight
    x_barre=zeros(length(u),1);%a presynaptic spike increases the trace x of some biophysical quantity, then decays exp.
    ltd_veto=zeros(length(u),1);
    ltd_veto(1)=0;
    u_bare=zeros(length(u),1); %low-pass filtered version of the postsynaptic membrane potential
    u_bare_ltd=zeros(length(u),1);
    theta_moins_ltd=zeros(length(u),1);
    
    x_barre(1)=0;
    u_bare(1)=Ed;
    u_bare_ltd(1)=Ed;
    theta_moins_ltd(1)=y(4);
    
    wi(1)=0.5;
    
    tsim=2;
    
    timeexp=[];
    while tsim<=length(u)
       
        if titi(tsim)==1
            timeexp=[timeexp tsim];
        end
        sum_exp=0;
        for kti=1:length(timeexp)
            sum_exp=sum_exp+exp((timeexp(kti)-tsim)*dt/y(1));
        end
        
        x_barre(tsim)=sum_exp;
        
        u_bare(tsim)=u_bare(tsim-1)+dt*(-u_bare(tsim-1)+u(tsim-1))/y(2);
        u_bare_ltd(tsim)=u_bare_ltd(tsim-1)+dt*(-u_bare_ltd(tsim-1)+u(tsim-1))/y(7);
        BB=u_bare_ltd(tsim-1)>theta_moins_ltd(tsim-1);
        AA=u_bare(tsim-1)>theta_moins_ltp;
        toto(tsim-1)=y(5)*x_barre(tsim-1)*AA*(u_bare(tsim-1)-theta_moins_ltp);
        
        ltd_veto(tsim)=ltd_veto(tsim-1)+dt*(-ltd_veto(tsim-1)+y(8)*(toto(tsim-1)))/y(9);
        theta_moins_ltd(tsim)=y(4)+ltd_veto(tsim);
        
        wi(tsim)=wi(tsim-1)+y(5)*x_barre(tsim-1)*AA*(u_bare(tsim-1)-theta_moins_ltp)*dt-y(6)*x_barre(tsim-1)*BB*(u_bare_ltd(tsim-1)-theta_moins_ltd(tsim-1))*dt;
        
        tsim=tsim+1;
    end
    
    dw(expdata)=wi(length(wi));
    norm_weight(expdata)=(0.5+nb_pairs(expdata).*(dw(expdata)-ww(1)))./ww(1);
    
     clear u v
end

if set==1
    norm_w_exp=[92 129 81 99 118 100 137 85 98]/100;
    sd=[1 1 1 1 1 1 1 1 1];
elseif set==2
    norm_w_exp=[90 97 93 115 75 129 62 69 57 65 163 154]/100;
    sd=[1 1 1 1 1 1 1 1 1 1 1 1];
end

S=sum((norm_weight'-norm_w_exp).^2./sd.^2)

end

