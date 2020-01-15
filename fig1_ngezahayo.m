%%Ngezahao

int1=input('Stimulation frequency? 2, 20 or 40 Hz?');

dt=0.1;
if int1==2%interval %ms
    int=500;
    nb_pair=2;
elseif int1==20 %interval %ms    
    int=50;
    nb_pair=100;
elseif int1==40 %interval %ms
    int=25;
    nb_pair=100;
end

nb_tot=nb_pair*int;
inter=length([dt:dt:int]);
Ed=0;

volta=[0:0.2:35];
data=length(volta);
dw=zeros(data,1);
dw_ltd=zeros(data,1);
dw_ltp=zeros(data,1);
norm_weight=zeros(data,1);

theta_moins_ltp=y(3);

for expdata=1:data
    
    tspan=[dt:dt:nb_tot+500]; % ms
    tproto=length(tspan);
    wi=zeros(tproto,1);
    wi_sum=zeros(tproto,1); %weight
    w_ltd=zeros(tproto,1); %weight
    w_ltp=zeros(tproto,1);
    
    u=zeros(tproto,1);
    titi=zeros(tproto,1);
    toto=zeros(tproto,1);
    
    u(1:tproto)=volta(expdata)*ones(tproto,1);
    
    for pp=0:nb_pair-1    
        titi(5000+pp*inter)=1; %pre-post
    end
    
    x_barre=zeros(tproto,1);%a presynaptic spike increases the trace x of some biophysical quantity, then decays exp.
    x_barre_euler=zeros(tproto,1);
    
    ltd_veto=zeros(tproto,1);
    theta_moins_ltd=zeros(tproto,1);
    u_bare=zeros(tproto,1); %low-pass filtered version of the postsynaptic membrane potential
    u_bare_ltd=zeros(tproto,1);
    
    x_barre(1)=0;
    ltd_veto(1)=0;
    theta_moins_ltd(1)=y(4);
    u_bare(1)=Ed;
    u_bare_ltd(1)=Ed;
    
    wi(1)=0.5;
    w_ltd(1)=0;
    w_ltp(1)=0;
    wi_sum(1)=0.5;
    ltp=0;
    ltd=0;
    
    tsim=2;
    
    timeexp=[];
    while tsim<=tproto
        
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
        toto(tsim-1)=y(5)*x_barre(tsim-1)*AA*(u_bare(tsim-1)-theta_moins_ltp); %was tsim=tsim-1...
        
        ltd_veto(tsim)=ltd_veto(tsim-1)+dt*(-ltd_veto(tsim-1)+y(8)*(toto(tsim-1)))/y(9);
        theta_moins_ltd(tsim)=y(4)+ltd_veto(tsim);
        
        wi(tsim)=wi(tsim-1)+y(5)*x_barre(tsim-1)*AA*(u_bare(tsim-1)-theta_moins_ltp)*dt-y(6)*x_barre(tsim-1)*BB*(u_bare_ltd(tsim-1)-theta_moins_ltd(tsim-1))*dt;
        
        tsim=tsim+1;
        
    end
    
    dw(expdata,1)=wi(length(wi));
    x_barre(length(x_barre));
    max(theta_moins_ltd)-y(4);
    
end

ww=0.5*ones(data,1);

if nb_pair==2
    norm_weight=(0.5+50.*(dw-ww))./ww;
else
    norm_weight=(dw./ww);
end

figure,
plot(volta,norm_weight*100);



