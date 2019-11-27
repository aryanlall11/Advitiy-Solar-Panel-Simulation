function [m_p,V_m,I_m] = allp(Iph,n_STC,N_Cell,Vt,Isat,Rs,G,Isc)
%Iph=[0,0,0,0];n_STC=2.6492;N_Cell=2;Vt=0.0278;Isat=3.159e-15;Rs=0.19;G=[0,0,0,0];Isc=[0,0,0,0];
pi=0;
pow=[];
k=0;
V=5.2;
while k==0
  syms I;
  if(G(1)==0)
   eq=I==0;
   m1=0;
  else
   eq=V-n_STC*N_Cell*Vt*log(((Iph(1)-I)/Isat)+1)+I*Rs==0;
   m1=Isc(1);
  end
  cur1=vpasolve(eq,I);
  
  syms I;
  if(G(2)==0)
   eq=I==0;
   m2=0;
  else
   eq=V-n_STC*N_Cell*Vt*log(((Iph(2)-I)/Isat)+1)+I*Rs==0;
   m2=Isc(2);
  end
  cur2=vpasolve(eq,I);
  
  syms I;
  if(G(3)==0)
   eq=I==0;
   m3=0;
  else
   eq=V-n_STC*N_Cell*Vt*log(((Iph(3)-I)/Isat)+1)+I*Rs==0;
   m3=Isc(3);
  end
  cur3=vpasolve(eq,I);
  
  syms I;
  if(G(4)==0)
   eq=I==0;
   m4=0;
  else
   eq=V-n_STC*N_Cell*Vt*log(((Iph(4)-I)/Isat)+1)+I*Rs==0;
   m4=Isc(4);
  end
  cur4=vpasolve(eq,I);
  if(cur1<0)
      cur1=0;
  end
  if(cur2<0)
      cur2=0;
  end
   if(cur3<0)
      cur3=0;
   end
   if(cur4<0)
      cur4=0;
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
  if isempty(cur4) 
    cur4=0;
  end
  I1=cur1+cur2+cur3+cur4;
  if(isreal(I1)==1 && I1<=(m1+m2+m3+m4))
   pow=[pow;V*I1];
%    plot(V,V*I1,'-x');
%    hold on
  if(V*I1>=pi)
      pi=V*I1;
      V=V-0.1;
      if(V<2)
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
% 
% function m_p = Test(Iph1,Iph2,Iph3,Iph4,m,n_STC,N_Cell,Vt,Isat,Rs,I)
% pi=0;
% pow=[];
% k=0;
% V=2;
% while k==0
%   syms I;
%   eq=V-n_STC*N_Cell*Vt*log((((Iph1-I)/Isat)+1)*(((Iph2-I)/Isat)+1))+2*I*Rs==0;
%   cur1=vpasolve(eq,I);
%   syms I;
%   eq=V-n_STC*N_Cell*Vt*log((((Iph3-I)/Isat)+1)*(((Iph4-I)/Isat)+1))+2*I*Rs==0;
%   cur2=vpasolve(eq,I);
%   if(cur1<0)
%       cur1=0;
%   end
%   if(cur2<0)
%       cur2=0;
%   end
%   I=cur1+cur2;
%   if(isreal(I)==1 && I<m)
%    pow=[pow;V*I];
%   if(V*I>=pi)
%       pi=V*I;
%       V=V+0.1;
%   else
%       k=1;
%       m_p=max(pow);
%   end
%   else
%      V=V+0.3;
%   end
% end
