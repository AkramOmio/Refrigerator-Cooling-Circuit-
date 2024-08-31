import CoolProp.CoolProp as CP
import numpy as np

'''Input Value'''
'''Compressor Input'''
#mass_flow_1=float(input('Mass Flow Rate (Kg/h):'))
mass_flow_1=1.45
mass_flow_2=mass_flow_1/3600
R_air_con=0.17646
R_air_eva=0.08
'''Condenser Input'''
# T1=float(input('D2 Temperature (C):'))
T1=91.2
T1=T1+273
#Tcond=float(input('Condenser Temperature (C):'))
Tcond=44.447
Tcond=Tcond+273
Pcond=CP.PropsSI('P','T',Tcond,'Q',1,'IsoButane')
print(Pcond)
#Tstainer=float(input('Degree of Strainer (C):'))
Tstainer=43.4
Tstainer=Tstainer+273
Tsub=Tcond-Tstainer
#Tamb=float(input('Ambient Temperature (C):'))
Tamb=38
Tamb=Tamb+273

#D_cond=float(input('Condenser Diameter (mm):'))
D_cond=3.2
D_cond=D_cond*1e-3
T2=CP.PropsSI('T','P',Pcond,'Q',1,'IsoButane')
T3=Tcond-Tsub
'''Capillary Tube Input'''
#D_cap=float(input('Capillary Tube Dia (mm):'))
D_cap=0.8
D_cap=D_cap*1e-3

'''Evaporator Input'''
#Teva=float(input('Evaporator Temperature (C):'))
Teva=-30.2
Teva=Teva+273
Peva=CP.PropsSI('P','T',Teva,'Q',1,'IsoButane')
#D_eva=float(input('Evaporattor Tube Diameter (mm):'))
D_eva=5.4
D_eva=D_eva*1e-3
#T_FC=float(input('FC Temperature:'))
T_FC=-26.9
T_FC=T_FC+273
#T_PC=float(input('PC Temperature:'))
T_PC=-1.4
T_PC=T_PC+273
#Q_FC=float(input('Insulation Heat in FC:'))
Q_FC=43.2
#Q_PC=float(input('Insulation Heat in PC:'))
Q_PC=25.9
Vfc=0.4
Vpc=0.6

'''  '''
'''Condenser Calculation'''
def condenser_length(Pcond,T1,Tcond,T3,mass_flow_2):
    '''Heat Calculation'''
    H1=CP.PropsSI('H','P',Pcond,'T',T1,'IsoButane')
    H2=CP.PropsSI('H','P',Pcond,'Q',1,'IsoButane')
    H3=CP.PropsSI('H','P',Pcond,'Q',0,'IsoButane')
    H4=CP.PropsSI('H','P',Pcond,'T',T3,'IsoButane')
    Q1=mass_flow_2*(H1-H2)
    Q2=mass_flow_2*(H2-H3)
    Q3=mass_flow_2*(H3-H4)
    
    '''Super Heated Zone'''
    LMTD=((T1-Tamb)-(Tcond-Tamb))/(np.log(((T1-Tamb)/(Tcond-Tamb))))
    k_vapor=CP.PropsSI('L','P',Pcond,'T',T1,'IsoButane')
    h_vapor=(4.36*k_vapor)/D_cond
    U_vapor=(1/(R_air_con+(2.58/h_vapor)))
    L_vapor=((Q1)/(U_vapor*np.pi*D_cond*LMTD))*0.2
    print('Super Heated Length %0.3f'%L_vapor)
    
    '''Two Phase Length'''
    viscosity_liquid=CP.PropsSI('V','P',Pcond,'Q',0,'IsoButane')
    Cp_liquid=CP.PropsSI('C','P',Pcond,'Q',0,'IsoButane')
    k_liquid=CP.PropsSI('L','P',Pcond,'Q',0,'IsoButane')
    Pr_liquid=(viscosity_liquid*Cp_liquid)/k_liquid
    Re_liquid=(4*mass_flow_2)/(D_cond*viscosity_liquid*np.pi)
    h_tp=0.023*(Re_liquid**0.8)*(Pr_liquid**0.4)*(k_liquid/D_cond)
    U_tp=h_tp/2.58
    dt=Tcond-Tamb
    L_tp=((Q2)/(U_tp*np.pi*D_cond*dt))
    print('Two Phase Length:%0.3f'%L_tp)
    
    '''Sub-cooled Length'''
    LMTD=((Tcond-Tamb)-(T3-Tamb))/(np.log(((Tcond-Tamb)/(T3-Tamb))))
    k_liquid=CP.PropsSI('L','P',Pcond,'Q',0,'IsoButane')
    h_liquid=(4.36*k_liquid)/D_cond
    U_liquid=h_liquid/2.58
    L_liquid=(Q3)/(U_liquid*np.pi*D_cond*LMTD)
    print('Sub-Coolde Length: %0.3f'%L_liquid)
    print('Total Length: %0.3f'%(L_vapor+L_tp+L_liquid))
print('\n') 
print('250 Circuit Report')    
print('\n') 
print('Condenser Result\n')    
condenser_length (Pcond,T1,Tcond,T3,mass_flow_2)


'''Capillary Tube Calculation'''
def capillary_length(P4,P2,dt,d,m):
    A=((np.pi*d**2)/4)
    G=m/A
    T1=CP.PropsSI('T','P',P4,'Q',0,'IsoButane')
    T2=T1-dt
    h3=CP.PropsSI('H','P',P4,'Q',0,'IsoButane')
    cp=CP.PropsSI('C','P',P4,'Q',0,'IsoButane')
    h1=h3-cp*dt
    P1=CP.PropsSI('P','T',T2,'Q',0,'IsoButane')
    D4=CP.PropsSI('D','P',P4,'Q',0,'IsoButane')
    v4=1/D4
    h2=h1=CP.PropsSI('H','P',P1,'Q',0,'IsoButane')
    x=CP.PropsSI('Q','P',P2,'H',h2,'IsoButane')
    vis_1=CP.PropsSI('V','P',P1,'Q',0,'IsoButane')
    vis_l=CP.PropsSI('V','P',P2,'Q',0,'IsoButane')
    vis_v=CP.PropsSI('V','P',P2,'Q',1,'IsoButane')
    vis_a_star=((x*vis_v)+((1-x)*vis_l))/vis_1
    vis_a=vis_a_star*vis_1
    Re=(G*d)/vis_a
    B=(37530/Re)**16
    j=3.27e-4
    i=(7/Re)**0.9+0.27*j
    A=(2.457*np.log(1/i))**16
    s=(8/Re)**12+(1/((A+B)**(3/2)))
    lamb=8*(s**(1/12))
    L41=(2*(P4-P1)*d)/(lamb*(G**2)*v4)
   
    k1=(2.62*(10**5))/(P1**0.75)
    D1=CP.PropsSI('D','P',P1,'Q',0,'IsoButane')
    v1=1/D1
    k2=(v1*(G**2))/P1
    p_star=P2/P1
    u=p_star-1-((k1/(1-k1))*np.log(k1+((1-k1)*p_star)))
    w=1/(k2*(1-k1))
    z=w*u
    y=np.log(p_star/(k1+(1-k1)*p_star))
    L_star=y-z
    L12=(L_star*d*2)/lamb
    print('Adibatic Length: %0.3f'%L41)
    L=L12+L41
    print('Non Adibatic Length: %0.3f'%L12)
    print('Total Length (m):%0.3f'%L)
    return L41,L12,L,x
print('\n')  
print('Capillary Tube Result\n')    
capillary_lsp,capillary_ltp,capillary_l,cap_x=capillary_length(Pcond,Peva,Tsub,D_cap,mass_flow_2)

'''Evaporator Calculation'''
'''Heat Calculation'''
def evaporator_length (Peva,T_FC,T_PC,Q_FC,Q_PC,mass_flow_2,D_eva):
    
    viscosity_liquid=CP.PropsSI('V','P',Peva,'Q',0,'IsoButane')
    viscosity_vapor=CP.PropsSI('V','P',Peva,'Q',1,'IsoButane')
    density_liquid=CP.PropsSI('D','P',Peva,'Q',0,'IsoButane')
    density_vapor=CP.PropsSI('D','P',Peva,'Q',1,'IsoButane')
    x=cap_x
    X_tt=(density_vapor/density_liquid)**0.5*(viscosity_liquid/viscosity_vapor)**0.1*((1-x)/x)**0.9
    h1=CP.PropsSI('H','P',Peva,'Q',0,'IsoButane')
    h2=CP.PropsSI('H','P',Peva,'Q',1,'IsoButane')
    hfg=h2-h1
    k_liquid=CP.PropsSI('L','P',Peva,'Q',0,'IsoButane')
    
    '''FC Calculation'''
    dt_FC=Vfc*(T_FC-Teva)
    h_tp=X_tt*((density_liquid*(density_liquid-density_vapor)*9.81*hfg*k_liquid**3)/(D_eva*dt_FC*viscosity_liquid))**0.25
    U_tp=1/(R_air_eva+5/h_tp)
    L_FC=((Q_FC)/(U_tp*np.pi*D_eva*dt_FC))*0.1
    print('FC Length (m): %0.3f'%L_FC)
    
    '''PC Calculation'''
    dt_PC=Vpc*(T_PC-Teva)
    h_tp=X_tt*((density_liquid*(density_liquid-density_vapor)*9.81*hfg*k_liquid**3)/(D_eva*dt_PC*viscosity_liquid))**0.25
    U_tp=1/(R_air_eva+5/h_tp)
    L_PC=(Q_FC)/(U_tp*np.pi*D_eva*dt_PC)*0.1
    print('PC length (m):%0.3f'%L_PC)
    print('Total length (m):%0.3f'%(L_PC+L_FC))
    
print('\n')  
print('Evaporator Result\n')  
evaporator_length (Peva,T_FC,T_PC,Q_FC,Q_PC,mass_flow_2,D_eva)