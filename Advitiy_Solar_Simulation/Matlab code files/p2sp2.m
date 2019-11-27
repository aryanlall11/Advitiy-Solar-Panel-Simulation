function [m_p,V_m,I_m] = p2sp2(Iph,n_STC,N_Cell,Vt,Isat,Rs,G,Isc,Voc)
%Iph=[0,0,0.4414,0.5558];n_STC=2.6492;N_Cell=2;Vt=0.0278;Isat=3.159e-15;Rs=0.19;G=[0,0,776,1070];Isc=[0,0,0.4414,0.5558];Voc=[-5.8932e-20,-5.8932e-20,4.8034,4.8374];
pi=0;
pow=[];         % (1||2)s(3||4)
V1=0;
V2=0;
k=0;
if((sum(G(:) == 0))==4)
  I_m=0;  %non zero pair element
  V_m=4.2;
  m_p=0;
  k=1;
elseif(G(1)*G(4)==0 && G(2)*G(3)==0)
    if(G(1)==0)
        m1=Isc(4);
    else
        m1=Isc(1);
    end
    if(G(2)==0)
        m2=Isc(3);
    else
        m2=Isc(2);
    end
  I=min(m1,m2);
elseif(G(1)*G(4)==0)
    if(G(1)+G(4)==0)
        I=Isc(2)+Isc(3);
    elseif(G(1)==0)
        I=Isc(4);
    else
        I=Isc(1);
    end
 elseif(G(2)*G(3)==0)
    if(G(2)+G(3)==0)
        I=Isc(1)+Isc(4);
    elseif(G(2)==0)
        I=Isc(3);
    else
        I=Isc(2);
    end
end
%--------------------------------------------------------------------------
tic
while k==0
  if(toc<10)
      d=0.0001;
  else
      d=0.001;
  end
  cur1=0;
  if(G(1)==0 && G(4)==0)
   V1=0;
   m1=0;
  elseif(G(1)==0)
   V1=n_STC*N_Cell*Vt*log(((Iph(4)-I)/Isat)+1)-I*Rs;
   m1=Voc(4);
  elseif(G(4)==0)
   V1=n_STC*N_Cell*Vt*log(((Iph(1)-I)/Isat)+1)-I*Rs;
   m1=Voc(1);
  else
   syms I1;
   eq=n_STC*N_Cell*Vt*log(((Iph(1)-I1)/Isat)+1)-I1*Rs==n_STC*N_Cell*Vt*log(((Iph(4)-I+I1)/Isat)+1)-(I-I1)*Rs==0;
   cur1=double(vpasolve(eq,I1));
   V1=n_STC*N_Cell*Vt*log(((Iph(1)-cur1)/Isat)+1)-cur1*Rs;
   m1=min(Voc(1),Voc(4));
  end
  
  cur2=0;
  if(G(2)==0 && G(3)==0)
   V2=0;
   m2=0;
  elseif(G(2)==0)
   V2=n_STC*N_Cell*Vt*log(((Iph(3)-I)/Isat)+1)-I*Rs;
   m2=Voc(3);
  elseif(G(3)==0)
   V2=n_STC*N_Cell*Vt*log(((Iph(2)-I)/Isat)+1)-I*Rs;
   m2=Voc(2);
  else
   syms I2;
   eq=n_STC*N_Cell*Vt*log(((Iph(2)-I2)/Isat)+1)-I2*Rs==n_STC*N_Cell*Vt*log(((Iph(3)-I+I2)/Isat)+1)-(I-I2)*Rs==0;
   cur2=double(vpasolve(eq,I2));
   V2=n_STC*N_Cell*Vt*log(((Iph(2)-cur2)/Isat)+1)-cur2*Rs;
   m2=min(Voc(2),Voc(3));
  end
  %---------------------------------------------------------------------------------------
  
  if ~(isempty(cur1) || isempty(cur2))
   if(V1<0)
       V1=0;
   end
   if(V2<0)
       V2=0;
   end
   net_V=double(V1+V2);
%---------------------------------------------------------------------------------
 if(isreal(net_V)==1 && net_V<=(m1+m2))
   pow=[pow;net_V*I];
%     plot(net_V,I*net_V,'-x');
%     hold on
  if(net_V*I>=pi)
      pi=net_V*I;
      I=I-d;
      if(I<0.2)
          m_p=0;
          V_m=4.2;
          I_m=0;
          k=1;
      end
  else
      k=1;
      m_p=max(pow);
      V_m=net_V;
      I_m=I;
  end
  else
     I=I-d;
  end
 else
    I=I-d;
  end
end
end