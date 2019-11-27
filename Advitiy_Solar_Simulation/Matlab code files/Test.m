%------------------------------------------------2s-p-2s-------------------------------------------%
function [m_p,V_m,I_m] = Test(Iph,n_STC,N_Cell,Vt,Isat,Rs,G,Isc)
%Iph1=0.5065;Iph2=0.0037;Iph3=0.4414;Iph4=0.0037;n_STC=2.6492;N_Cell=2;Vt=0.0278;Isat=3.159e-15;Rs=0.19;G=[940,0,740,0];Isc=[0.5065,0.0037,0.4414,0.0037];
pi=0;
pow=[];
k=0;
V=9.5;
while k==0
  syms I;
  if(G(1)+G(4)==0)
   eq=I==0;
   m1=0.5;
  elseif(G(4)==0)
   eq=V-n_STC*N_Cell*Vt*log(((Iph(1)-I)/Isat)+1)+I*Rs==0;
   m1=Isc(1);
  elseif(G(1)==0)
   eq=V-n_STC*N_Cell*Vt*log(((Iph(4)-I)/Isat)+1)+I*Rs==0;
   m1=Isc(4);
  elseif(G(1)~=0 && G(4)~=0)
   eq=V-n_STC*N_Cell*Vt*log((((Iph(1)-I)/Isat)+1)*(((Iph(4)-I)/Isat)+1))+2*I*Rs==0;
   m1=min(Isc(1),Isc(4));
  end
  cur1=vpasolve(eq,I);
  
  syms I;
  if(G(2)+G(3)==0)
   eq=I==0;
   m2=0.5;
  elseif(G(2)==0)
   eq=V-n_STC*N_Cell*Vt*log(((Iph(3)-I)/Isat)+1)+I*Rs==0;
   m2=Isc(3);
  elseif(G(3)==0)
   eq=V-n_STC*N_Cell*Vt*log(((Iph(2)-I)/Isat)+1)+I*Rs==0;
   m2=Isc(2);
  elseif(G(2)~=0 && G(3)~=0)
   eq=V-n_STC*N_Cell*Vt*log((((Iph(2)-I)/Isat)+1)*(((Iph(3)-I)/Isat)+1))+2*I*Rs==0;
   m2=min(Isc(2),Isc(3));
  end
  cur2=vpasolve(eq,I);
  if(cur1<0)
      cur1=0;
  end
  if(cur2<0)
      cur2=0;
  end
  if isempty(cur2) 
    cur2=0;
  end
  if isempty(cur1) 
    cur1=0;
  end
  I1=cur1+cur2;
  if(isreal(I1)==1 && I1<(m1+m2))
   pow=[pow;V*I1];
%    plot(V,V*I1,'-x');
%    hold on
  if(V*I1>=pi)
      pi=V*I1;
      V=V-0.1;
      if(V<1)
          m_p=0;
          V_m=4.2;
          I_m=0;
          k=1;
      end
  else
      k=1;
      m_p=max(pow);
      V_m=V;
      I_m=I1;
      V=V-0.1;
  end
  else
     V=V-0.1;
  end
end
end
