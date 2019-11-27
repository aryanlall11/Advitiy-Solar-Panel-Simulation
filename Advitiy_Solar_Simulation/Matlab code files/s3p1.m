function [m_p,V_m,I_m] = s3p1(Iph,n_STC,N_Cell,Vt,Isat,Rs,G,Isc,Voc)
%Iph=[0.4791,0,0.4414,0.3985];n_STC=2.6492;N_Cell=2;Vt=0.0278;Isat=3.159e-15;Rs=0.19;G=[870,0,776,673];Isc=[0.4791,0,0.4414,0.3985];Voc=[4.8155,-5.89322759461873e-20,4.8034,4.7883];
pi=0;
pow=[];
cur=[];
vol=[];
k=0;
V=5.3;
net_I=0;
if(G(1)==0)              %parallel panel
    l=2;
else
    l=Isc(1);
end
while k==0 && net_I<=l
  syms I;
  if(G(2)==0)
   eq=I==0;
   m1=0;
  else
   eq=V-n_STC*N_Cell*Vt*log(((Iph(2)-I)/Isat)+1)+I*Rs==0;
   m1=Isc(2);
  end
  cur1=vpasolve(eq,I);
  
  syms I;
  if(G(3)==0)
   eq=I==0;
   m2=0;
  else
   eq=V-n_STC*N_Cell*Vt*log(((Iph(3)-I)/Isat)+1)+I*Rs==0;
   m2=Isc(3);
  end
  cur2=vpasolve(eq,I);
  
  syms I;
  if(G(4)==0)
   eq=I==0;
   m3=0;
  else
   eq=V-n_STC*N_Cell*Vt*log(((Iph(4)-I)/Isat)+1)+I*Rs==0;
   m3=Isc(4);
  end
  cur3=vpasolve(eq,I);
  
  if(cur1<0)
      cur1=0;
  end
  if(cur2<0)
      cur2=0;
  end
   if(cur3<0)
      cur3=0;
   end
  
  if isempty(cur2) 
    cur2=0;
  end
  if isempty(cur1) 
    cur1=0;
  end
  if isempty(cur3) 
    cur3=0;
  end
 net_I=cur1+cur2+cur3;     %net current
 
%---------------------------------------------------------------------------------
 if(isreal(net_I)==1 && net_I<=(m1+m2+m3))
   net_I=double(net_I);
   syms V1;
  if(G(1)==0)
   eq=V1==0;
   m4=0;
  else
   eq=V1-n_STC*N_Cell*Vt*log(((Iph(1)-net_I)/Isat)+1)+net_I*Rs==0;
   m4=Voc(1);
  end
  vol1=vpasolve(eq,V1);
  
   if(vol1<0)
      vol1=0;
  end
  
  if isempty(vol1) 
    vol1=0;
  end
%---------------------------------------------------------------------------------  
  if(isreal(vol1)==1 && vol1<=abs(m4))
      
  net_V=V+double(vol1);
   pow=[pow;net_V*net_I];      %power
   cur=[cur;net_I];
   vol=[vol,net_V];
    %plot(net_V,net_I,'-x');
    %hold on
  if(net_V*net_I>=pi)
      pi=net_V*net_I;
      V=V-0.05;
      if(V<1)
          m_p=0;
          V_m=4.2;
          I_m=0;
          k=1;
      end
  else
      k=1;
      m_p=max(pow);
      V_m=net_V;
      I_m=net_I;
      V=V-0.05;
  end
  else
     V=V-0.05;
  end
 else
    V=V-0.05;
 end
end

 if(k==0)
   m_p=max(pow);
   V_m=vol(end);
   I_m=cur(end);  
 end
end