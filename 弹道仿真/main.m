
clc
clear
global a b c m_0 d L p g y_2 y_og t_1 a_0 
i=1;
n=1;
    prompt = {'  a :------(o=at+a_0)--  ';'  a_0 :--';'  b:------(F=bt)--   ';'  c :------(m=m_0-ct)--   ';'  d :------(M=dm_0)--   ';'  m_0 :--';'  p :-----(p=m/l)--   ';'  L :------';'  g:--';'  pause(?) :';' hold off ( y / n ):';'   y_2 :---'};
    dlg_title = 'missile......';
    num_lines= 1;
    def     = {'-.015','.001','30','.2','6','60','1','1100','9.81','.05','y','5000'};
    inp = inputdlg(prompt,dlg_title,num_lines,def);

while(strcmp(inp,'')==0)
    hh  = strcmp(inp(11,1),'y');
    inp = str2mat( inp );
    a   = str2num(inp(1,:));
    a_0 = str2num(inp(2,:));
    b   = str2num(inp(3,:));
    c   = str2num(inp(4,:));
    d   = str2num(inp(5,:));
    m_0 = str2num(inp(6,:));
    p   = str2num(inp(7,:));
    L   = str2num(inp(8,:));
    g   = str2num(inp(9,:));
    pp  = str2num(inp(10,:));
    y_2 = str2num(inp(12,:));
    
    if(L < (m_0)/p)
       msgbox('L < l')
      break
    end
    
    % #0# 
    x_0=0 ;      y_0 = L/2-((m_0)*(L-(m_0)/p))/(2*((1+d)*m_0));    fi_0=90;   t_0=0;  m_0;

    % #1# 
    y_og1=((m_0-c*t_1)*(L-(m_0-c*t_1)/p))/(2*((1+d)*m_0-c*t_1));
    x_1=0;   y_1=L/2-y_og1;    fi_1=90;   t_1=m_0*(1+d)/(b/g+c);  m_1=m_0-c*t_1;
  
   v_22=-13/10*(1175/2+13/40*t_1)/(600-13/5*t_1)+13/40*(50-13/10*t_1)/(600-13/5*t_1)+13/5*(50-13/10*t_1)*(1175/2+13/40*t_1)/(600-13/5*t_1)^2;
   
    
    % #2#
    tt=1;
    t_2 = t_1 + tt;
    [T_2,Y_2] = ode45(@state2,[t_1 t_2],[L/2-y_og1 v_22]);
    while(Y_2(length(Y_2(:,1)),1) < (y_2 + y_1))
        tt = tt + 1;t_2 = t_1 + tt;
        [T_2,Y_2] = ode45(@state2,[t_1 t_2],[L/2-y_og1 v_22]);
        
    end
    
    x_2 = 0;       fi_2 = 90;   m_2 = m_0-c*t_2;  y_2 = Y_2(length(Y_2(:,1)),1);
    
    v_y2 = Y_2(length(Y_2(:,2)),2);
  
    if(m_2 < 0)
         msgbox('fule is low (1)')
    break
    end

    % #3#
    t_3 = m_2/c;mm=0;
tt_3=1;
    [T_3,Y_3]=ode45(@state3,[t_2 t_2+tt_3],[0 y_2 90 0 v_y2 0]);
    while((Y_3(length(Y_3(:,5)),1) > 0))
        tt_3 = tt_3 + 1;t_2 = t_1 + tt_3;
        [T_3,Y_3]=ode45(@state3,[t_2 t_2+tt_3],[0 y_2 90 0 v_y2 0]);
        mm=mm+1;
   end
    
    if(t_3 < tt_3)
         msgbox('fule is low (2)')
    break
    end
       
    x_3 = Y_3(length(Y_3(:,1)),1);
    y_3 = Y_3(length(Y_3(:,2)),2);
    fi_3 = Y_3(length(Y_3(:,3)),3);
    v_x3 = Y_3(length(Y_3(:,4)),4);
    v_y3 = Y_3(length(Y_3(:,5)),5);
    w_fi3 = Y_3(length(Y_3(:,6)),6);

    % #4#
    t_4 = m_0/c;
    [T_4,Y_4]=ode45(@state4,[t_2+tt_3 t_4],[x_3 y_3 fi_3 v_x3 v_y3 w_fi3]);
    x_4 = Y_4(length(Y_4(:,1)),1);
    y_4 = Y_4(length(Y_4(:,2)),2);
    fi_4 = Y_4(length(Y_4(:,3)),3);
    v_x4 = Y_4(length(Y_4(:,4)),4);
    v_y4 = Y_4(length(Y_4(:,5)),5);
    w_fi4 = Y_4(length(Y_4(:,6)),6);

    
    
    % #5#

    tt=1;
    [T_5,Y_5]=ode45(@state5,[t_4 t_4+tt],[x_4 y_4 fi_4 v_x4 v_y4 w_fi4]);
    while(Y_5(length(Y_5(:,2)),2) > 0)
        [T_5,Y_5]=ode45(@state5,[t_4 t_4+tt],[x_4 y_4 fi_4 v_x4 v_y4 w_fi4]);
        tt = tt +1;
        
    end
 

   plot(zeros(length(Y_2),1),Y_2(:,1)  ,Y_3(:,1),Y_3(:,2) ,'-' ,Y_4(:,1),Y_4(:,2),'--' ,Y_5(:,1),Y_5(:,2),'b--' ),xlabel('x'),ylabel('y'),grid
   pause
   
   T_1=[0:1:t_1]';
   y_og2=L/2-((m_0-c*T_1).*(L-(m_0-c*T_1)/p))./(2*((1+d)*m_0-c*T_1)) ;
   v_2=-13./10*(1175./2+13./40.*T_1)./(600-13./5.*T_1)+13./40.*(50-13./10.*T_1)./(600-13./5.*T_1)+13./5.*(50-13./10.*T_1).*(1175./2+13./40.*T_1)./(600-13./5.*T_1).^2;
   GG=zeros(length(T_1),1);
   
   subplot(2,3,1),plot(T_1,GG,T_2,zeros(length(Y_2),1),T_3,Y_3(:,1),'-g',T_4,Y_4(:,1),'--'),ylabel('x'),xlabel('t'),grid
   subplot(2,3,2),plot(T_1,  y_og2 ,T_2,Y_2(:,1),T_3,Y_3(:,2) ,'-g',T_4,Y_4(:,2),'--'),ylabel('y'),xlabel('t'),grid
   subplot(2,3,3),plot(T_1,GG,T_2,90*ones(length(Y_2),1),T_3,Y_3(:,3),'-g',T_4,Y_4(:,3),'--'),ylabel('fi'),xlabel('t'),grid   
   subplot(2,3,4),plot(T_1,GG,T_2,zeros(length(Y_2),1),T_3,Y_3(:,4),'-g',T_4,Y_4(:,4),'--'),ylabel('dx'),xlabel('t'),grid
   subplot(2,3,5),plot(T_1,v_2,T_2,Y_2(:,2),T_3,Y_3(:,5),'-g',T_4,Y_4(:,5),'--'),ylabel('dy'),xlabel('t'),grid
   subplot(2,3,6),plot(T_1,GG,T_2,zeros(length(Y_2),1),T_3,Y_3(:,6),'-g',T_4,Y_4(:,6),'--'),ylabel('dfi'),xlabel('t'),grid


   pause
   cc=Y_3(1,2);
   i=1;
   while(abs((Y_3(i,2)-cc)) < 2)
       i=i+1;
       
   end
 
   X=[GG;zeros(length(Y_2(:,1)),1);Y_3(i:length(Y_3(:,2)),1);Y_4(:,1);Y_5(:,1)];
   Y=[y_og2;Y_2(:,1);Y_3(i:length(Y_3(:,2)),2);Y_4(:,2);Y_5(:,2)];
   FI=[90*ones(length(Y_2(:,1))+length(T_1),1) ; Y_3(i:length(Y_3(:,2)),3) ; Y_4(:,3);Y_5(:,3)];
   V_x=[GG;zeros(length(Y_2(:,1)),1);Y_3(i:length(Y_3(:,2)),4);Y_4(:,4);Y_5(:,4)];
   V_y=[v_2;Y_2(:,2);Y_3(i:length(Y_3(:,2)),5);Y_4(:,5);Y_5(:,5)];
   T=[T_1;T_2;T_3(i:length(T_3),1);T_4;T_5];
   X_max=X(1,1);   Y_max=Y(1,1); X_min=X(1,1);   Y_min=Y(1,1);
   
   for i=1:1:length(X)
       if(X(i,1) > X_max)
           X_max=X(i,1);
       end
       if(X(i,1) < X_min)
           X_min=X(i,1);
       end
       if(Y(i,1) > Y_max)
           Y_max=Y(i,1);
       end
       if(Y(i,1) < Y_min)
           Y_min=Y(i,1);
       end
   end
%-------------------------1---------------------------   
subplot(1,1,1)
for i=1:1:length(Y)
      if(T(i,1) <= (t_2 + t_3)) 
          y_og = ((m_0-c*T(i,1))*(L-(m_0-c*T(i,1))/p))/(2*((1+d)*m_0-c*T(i,1))) ;     
          L_m = (m_0-c*T(i,1))/p;
          aa = 7/200;
          L_w = (1-c*T(i,1)/(m_0*(1+d)))*L/5;
          L_F = L*b*T(i,1)/(5*m_0*(1+d)*g);
      else
          y_og=0;
          L_m=0;aa=0;
          L_w=L/(5*(1+1/d));
          L_F=0;               
      end
      if((T(i,1) <= (t_2))|(T(i,1) >= (t_2+t_3))) 
          cc=0;
      else 
          cc=a*(T(i,1)-t_2);
      end
     plot([X(i,1)-(L/2-y_og)*cos( FI(i,1)*pi/180 ),X(i,1)+(L/2+y_og)*cos( FI(i,1)*pi/180 )],[Y(i,1)-(L/2-y_og)*sin(FI(i,1)*pi/180),Y(i,1)+(L/2+y_og)*sin(FI(i,1)*pi/180)]);
     xlim([X_min - L , X_max + L])
     ylim([Y_min - L , Y_max + L])
     hold on
     XX_0 = X(i,1) - (L/2-y_og)*cos(FI(i,1)*pi/180);
     YY_0 = Y(i,1) - (L/2-y_og)*sin(FI(i,1)*pi/180);
     XX_1 = XX_0 - aa*L*sin( FI(i,1)*pi/180 );
     XX_2 = XX_1 + L_m*cos( FI(i,1)*pi/180 );
     XX_3 = XX_2 + aa*2*L*sin( FI(i,1)*pi/180 );
     XX_4 = XX_3 - L_m*cos( FI(i,1)*pi/180 );
     YY_1 = YY_0 + aa*L*cos( FI(i,1)*pi/180 );
     YY_2 = YY_1 + L_m*sin( FI(i,1)*pi/180 );
     YY_3 = YY_2 - aa*2*L*cos( FI(i,1)*pi/180 );
     YY_4 = YY_3 - L_m*sin( FI(i,1)*pi/180 );
     plot(X(i,1),Y(i,1),'+g')
     plot([XX_1,XX_2,XX_3,XX_4,XX_1],[YY_1,YY_2,YY_3,YY_4,YY_1],'r');
%fill([XX_1,XX_2,XX_3,XX_4,XX_1],[YY_1,YY_2,YY_3,YY_4,YY_1],'r');
     bb=2/15;     
%----------w
     plot([X(i,1),X(i,1)    , X(i,1)-bb*L_w/((2)^.5)     , X(i,1)-bb*L_w/((2)^.5)+(2)^.5*bb*L_w , X(i,1)],...
          [Y(i,1),Y(i,1)-L_w, Y(i,1)-L_w+bb*L_w/((2)^.5) , Y(i,1)-L_w+bb*L_w/((2)^.5)           , Y(i,1)-L_w],'k')
%----------F
    plot([XX_0 , XX_0-bb*L_F*cos((cc+FI(i,1)-45)*pi/180) , XX_0+bb*L_F*sin((cc+FI(i,1)-45)*pi/180) , XX_0 , XX_0-L_F*cos((cc+FI(i,1))*pi/180)],...
         [YY_0 , YY_0-bb*L_F*sin((cc+FI(i,1)-45)*pi/180) , YY_0-bb*L_F*cos((cc+FI(i,1)-45)*pi/180) , YY_0 , YY_0-L_F*sin((cc+FI(i,1))*pi/180)],'k')
%----------f_x
     L_vx = L/5*(sign(V_x(i,1)).*0.1.*V_x(i,1).^2)/(m_0*g*(1+d));
    plot([X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2 , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2                      , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) + L_vx],...
         [Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2 , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2 - bb*((2)^(.5))*L_vx , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 )       ],'r')
%----------f_y
     L_vy = L/5*(sign(V_y(i,1)).*0.1.*V_y(i,1).^2)/(m_0*g*(1+d));
     plot([X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) - L_vy*bb*((2)^(.5))/2 , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) - L_vy*bb*((2)^(.5))/2 + bb*((2)^(.5))*L_vy , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 )       ]...
         ,[Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vy*bb*((2)^(.5))/2 , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vy*bb*((2)^(.5))/2                      , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vy],'r')
     if(hh==1)
         hold off
     end
 pause(pp+.03)
end
% --------------2-----------------------------   

 for i=1:1:length(Y)
      if(T(i,1) <= (t_2 + t_3)) 
          y_og=((m_0-c*T(i,1))*(L-(m_0-c*T(i,1))/p))/(2*((1+d)*m_0-c*T(i,1))) ;     
          L_m=(m_0-c*T(i,1))/p; aa=7/200;
          L_w=(1-c*T(i,1)/(m_0*(1+d)))*L/5;
          L_F=L*b*T(i,1)/(5*m_0*(1+d)*g);
      else
          y_og=0;
          L_m=0;aa=0;
          L_w=L/(5*(1+1/d));
          L_F=0;               
      end
      if((T(i,1) <= (t_2))|(T(i,1) >= (t_2+t_3))) 
          cc=0;
      else 
          cc=a*(T(i,1)-t_2);
      end
     
    plot([X(i,1)-(L/2-y_og)*cos( FI(i,1)*pi/180 ),X(i,1)+(L/2+y_og)*cos( FI(i,1)*pi/180 )],[Y(i,1)-(L/2-y_og)*sin(FI(i,1)*pi/180),Y(i,1)+(L/2+y_og)*sin(FI(i,1)*pi/180)]);
    hold on
    xlim([X(i,1)-L X(i,1)+L])
    ylim([Y(i,1)-L Y(i,1)+L])
    XX_0=X(i,1)-(L/2-y_og)*cos(FI(i,1)*pi/180);
    YY_0=Y(i,1)-(L/2-y_og)*sin(FI(i,1)*pi/180);
    XX_1 = XX_0 - aa*L*sin( FI(i,1)*pi/180 );
    XX_2 = XX_1 + L_m*cos( FI(i,1)*pi/180 );
    XX_3 = XX_2 + aa*2*L*sin( FI(i,1)*pi/180 );
    XX_4 = XX_3 - L_m*cos( FI(i,1)*pi/180 );
    YY_1 = YY_0 + aa*L*cos( FI(i,1)*pi/180 );
    YY_2 = YY_1 + L_m*sin( FI(i,1)*pi/180 );
    YY_3 = YY_2 - aa*2*L*cos( FI(i,1)*pi/180 );
    YY_4 = YY_3 - L_m*sin( FI(i,1)*pi/180 );
    plot(X(i,1),Y(i,1),'+g')
    plot([XX_1,XX_2,XX_3,XX_4,XX_1],[YY_1,YY_2,YY_3,YY_4,YY_1],'r');
%fill([XX_1,XX_2,XX_3,XX_4,XX_1],[YY_1,YY_2,YY_3,YY_4,YY_1] ,'r');
  bb=2/15;
%----------w
     plot([X(i,1),X(i,1)    , X(i,1)-bb*L_w/((2)^.5)     , X(i,1)-bb*L_w/((2)^.5)+(2)^.5*bb*L_w , X(i,1)],...
          [Y(i,1),Y(i,1)-L_w, Y(i,1)-L_w+bb*L_w/((2)^.5) , Y(i,1)-L_w+bb*L_w/((2)^.5)           , Y(i,1)-L_w],'k')
%----------F
     plot([XX_0 , XX_0-bb*L_F*cos((cc+FI(i,1)-45)*pi/180) , XX_0+bb*L_F*sin((cc+FI(i,1)-45)*pi/180) , XX_0 , XX_0-L_F*cos((cc+FI(i,1))*pi/180)],...
          [YY_0 , YY_0-bb*L_F*sin((cc+FI(i,1)-45)*pi/180) , YY_0-bb*L_F*cos((cc+FI(i,1)-45)*pi/180) , YY_0 , YY_0-L_F*sin((cc+FI(i,1))*pi/180)],'k')
%----------f_x
     L_vx = L/5*(sign(V_x(i,1)).*0.1.*V_x(i,1).^2)/(m_0*g*(1+d));
     plot([X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2 , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2                      , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) + L_vx],...
          [Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2 , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vx*bb*((2)^(.5))/2 - bb*((2)^(.5))*L_vx , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 )       ],'r')
%----------f_y
     L_vy = L/5*(sign(V_y(i,1)).*0.1.*V_y(i,1).^2)/(m_0*g*(1+d));
     plot([X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) - L_vy*bb*((2)^(.5))/2 , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) - L_vy*bb*((2)^(.5))/2 + bb*((2)^(.5))*L_vy , X(i,1) + y_og*cos( FI(i,1)*pi/180 ) , X(i,1) + y_og*cos( FI(i,1)*pi/180 )       ]...
         ,[Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vy*bb*((2)^(.5))/2 , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vy*bb*((2)^(.5))/2                      , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) , Y(i,1) + y_og*sin( FI(i,1)*pi/180 ) + L_vy],'r')
     hold off
  pause(pp)
 end
pause

close
inp = inputdlg(prompt,dlg_title,num_lines,def);
n=n+1;
end %while
