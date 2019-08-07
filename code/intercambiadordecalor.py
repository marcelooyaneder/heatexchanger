# -*- coding: utf-8 -*-
#################MARCELO OYANEDER LABARCA##################
#############CONTACTO:MARCELO.OYANEDER@USACH.CL############

import numpy as np
import pandas as pd
from iapws import IAPWS95
import math
from CoolProp.CoolProp import PropsSI
import time

###########################################################
###################DECLARACION DE INPUTS###################
###########################################################
s=time.time()
archivo=pd.read_excel("inputs.xlsx",'inputs',index_col=None)
print 'las variables ingresadas al sistema son: '
print'****************************************'
print archivo
print'****************************************'
T_ev=archivo['T_ev[K]'].to_numpy()
T_sv=archivo['T_sv[K]'].to_numpy()
T_ea=archivo['T_ea[K]'].to_numpy()
T_sa=archivo['T_sa[K]'].to_numpy()
T_inf=archivo['T_inf[K]'].to_numpy()
Work_pressure=archivo['P[psig]'].to_numpy()
dif_altura_man=archivo['h[cm]'].to_numpy()
steam_mass_flow=archivo['Steam_mass_flow[kg/s]'].to_numpy()
T_prom_shell=archivo['T_prom_shell[K]'].to_numpy()

###########################################################
###################DECLARACION DE CONSTANTES###############
###########################################################
pi=np.pi
Inner_diameter_tube=0.01
Outer_diameter_tube=0.013
Inner_diameter_casquet=0.107
Outer_diameter_casquet=0.115
Tube_long=0.825
Casquet_long=0.8
#!son 7 tubos en total, pero son 2 etapas
Total_tubes=14.0
#!Areas, REVISAR!!!****
A_cil=pi*Outer_diameter_casquet*Casquet_long
A_tubo_int=pi*Inner_diameter_tube*Tube_long
A_tubo_ext=pi*Outer_diameter_tube*Tube_long

###########################################################
###################DECLARACION DE FUNCIONES################
###########################################################
def Air_physical_properties(T):
    Prandtl_air=PropsSI('Prandtl','T',T,'P',101325,'Air')
    Density_air=PropsSI('D','T',T,'P',101325,'Air')
    Viscosity_air=PropsSI('V','T',T,'P',101325,'Air')
    Thermal_conductivity_air=PropsSI('L','T',T,'P',101325,'Air')
    return np.array([Prandtl_air,Density_air,Viscosity_air,Thermal_conductivity_air])

def water_mass_flow(dif_altura_man):
    water_mass_flow=0.039333333*dif_altura_man+0.147888889
    #resultado en [kg/s]
    return water_mass_flow

def Grashof(T_prom_corasa, T_inf, T_film,Diameter,Density,Viscosity):
    g=9.8
    Dif_temp=T_prom_corasa-T_inf
    betha=(T_film+273.15)**(-1.0)
    Grashof=g*betha*Dif_temp*(Diameter**3.0)*(Density**2.0)*(Viscosity**(-2.0))
    return Grashof

def ChurchillChu(Gr,Pr):
    #!Se obtiene el numero de Nusselt^-1
    ChurchillChu=(0.6+0.387*((Gr*Pr)/((1.0+(0.559/Pr)**(9.0/16.0))**(16.0/9.0)))**(1.0/6.0))**2.0
    return ChurchillChu

def Reynolds(Mass_flow,Viscosity,Diameter):
    Reynolds=4*Mass_flow/(Viscosity*Diameter*pi*7.)
    return Reynolds

def Heat_transfer_water_iteration(Re_tuberia,Prandtl):
    Heat_transfer_water_iteration=0.023*pow(Re_tuberia,0.8)*pow(Prandtl,0.334)
    return Heat_transfer_water_iteration

def Chen_correlation(Densidad_l,Densidad_g,Entalpia,Conductividad,Viscosidad,T_sat,x,Outer_diameter_tube):
    g=9.8
    a=g*(Densidad_l*(Densidad_l-Densidad_g))*Entalpia*1000.*(Conductividad**3)
    b=Viscosidad*(T_sat-x)*4.*Outer_diameter_tube
    c=(a/b)**0.25
    Chen_correlation=0.725*c
    return Chen_correlation

def T_wall(Chen_correlation,h_correlacion_agua,T_delta_liq,T_delta_vapor,A_tubo_ext,A_tubo_int):
    a=A_tubo_ext*Chen_correlation*T_delta_vapor+A_tubo_int*h_correlacion_agua*T_delta_liq
    b=A_tubo_int*h_correlacion_agua+A_tubo_ext*Chen_correlation
    T_wall=a/b
    return T_wall

def Coef_calor_limpio(Outer_diameter_tube,Inner_diameter_tube,h_correlacion_agua,h_correlacion_vapor):
    #1/U_c
    a=Outer_diameter_tube/(Inner_diameter_tube*h_correlacion_agua)
    b=1/h_correlacion_vapor
    Coef_calor_limpio=a+b
    return Coef_calor_limpio

def DMLT(T_ea,T_sa,T_ev,T_sv):
    a=T_sv-T_ea
    b=T_ev-T_sa
    DMLT=(a-b)/math.log(a/b)
    return DMLT

def Coef_calor_sucio(A_tubo_int,DMLT,Calor_absorbido):
    Coef_calor_sucio=14.*A_tubo_int*DMLT*0.98/(Calor_absorbido*1.) 
    #Coef_calor_sucio=1.*A_tubo_int*DMLT*0.98/(Calor_absorbido*1000.)
    return Coef_calor_sucio


###########################################################
###################BLOQUE DE DESARROLLO####################
###########################################################

#Ingreso de la presión de trabajo y transformacion a [MPa]
Work_pressure=Work_pressure[0]
Work_pressure=Work_pressure+14.6959
Work_pressure_MPa=Work_pressure*0.00689476
print 'La presión de trabajo, corresponde a: ',Work_pressure_MPa,'[MPa]'

#Obtencion calor cedido por el vapor
T_ev=T_ev[0]
T_sv=T_sv[0]
T_prom_vapor=(T_ev+T_sv)/2.
steam_mass_flow=steam_mass_flow[0]
sat_liquid=IAPWS95(P=Work_pressure_MPa, x=0)
sat_steam=IAPWS95(P=Work_pressure_MPa, x=1)
T_sat=sat_steam.T
Vaporization_Enthalpy=sat_steam.h-sat_liquid.h
Specific_heat_liquid=sat_liquid.cp
Q_ced=(Vaporization_Enthalpy+Specific_heat_liquid*(T_ev-T_sv))*steam_mass_flow
print 'el calor cedido por el vapor es: ',Q_ced,'[KW]'

#Obtencion calor absorbido por el agua
T_ea=T_ea[0]
T_sa=T_sa[0]
dif_altura_man=dif_altura_man[0]
water_mass_flow=water_mass_flow(dif_altura_man)
T_prom_liq=(T_ea+T_sa)/2.
sat_liquid=IAPWS95(T=T_prom_liq, x=0)
Specific_heat_liquid=sat_liquid.cp
Q_abs=Specific_heat_liquid*(T_sa-T_ea)*water_mass_flow
print 'el calor absorbido por el liquido es: ',Q_abs,'[KW]'

#Obtencion calor perdido del sistema
Q_lost=Q_ced-Q_abs
print 'el calor perdido por el sistema es: ',Q_lost,'[KW]'

#calculos por correlaciones 
T_prom_shell=T_prom_shell[0]
T_inf=T_inf[0]
T_film=(T_prom_shell+T_inf)/2.
Air_physical_properties=Air_physical_properties(T_film)
#np.array([Prandtl_air,Density_air,Viscosity_air,Thermal_Conductivity_air])
Gr=Grashof(T_prom_shell,T_inf,T_film,Outer_diameter_casquet,Air_physical_properties[[1]],Air_physical_properties[2])
Nu_casquet=ChurchillChu(Gr,Air_physical_properties[0])
h_casquet=Nu_casquet*Air_physical_properties[3]/Outer_diameter_casquet
Q_lost_churchillchu=h_casquet*A_cil*(T_prom_shell-T_inf)
print 'el calor perdido calculado por la correlacion de Churchill & Chu es: ',Q_lost_churchillchu,'[W]'

#calculo de error.
error_calor_perdido=abs(1000*Q_lost-Q_lost_churchillchu)
print 'La diferencia entre energias perdidas calculada por ambos metodos es: ',error_calor_perdido
print 'La eficiencia del intercambiador de calor es: ', Q_abs*100/Q_ced,'%'

#Iteracion temperatura de pared
error=float(input('ingresar la diferencia a calcular para la iteracion: '))
x=T_prom_shell
dif_t_wall=10.
n=1
while (dif_t_wall>=error):
    print 'iteración numero: ',n
    #Agua
    t_film_ite=(T_prom_liq+x)/2.
    t_film_ite_steam=(T_sat+x)/2.
    print 'temperatura de film: ', t_film_ite,'[K]'
    print 'temperatura de pared: ', x,'[K]'
    sat_liquid=IAPWS95(T=t_film_ite, x=0)
    #vapor a P cte, por lo tanto temperatura constante.-
    a=sat_liquid.Prandt
    e=sat_liquid.k
    #Vaporization_Enthalpy=sat_steam.h-sat_liquid.h
    Re_tuberia=Reynolds(water_mass_flow,sat_liquid.mu,Inner_diameter_tube)
    Nussel_correlacion_agua=Heat_transfer_water_iteration(Re_tuberia,a)
    h_correlacion_agua=Nussel_correlacion_agua*e/Inner_diameter_tube
    #function Chen_correlation (Density_water, density_steam,Vaporization_enthalpy,Thermal_conductivity_water,Viscosity_water,T_sat,T_prom_corasa,Outer_diameter_tube)
    #propiedades del condensado
    sat_liquid=IAPWS95(T=t_film_ite_steam, x=0)
    b=sat_liquid.mu
    d=sat_liquid.rho
    e=sat_liquid.k
    h_correlacion_vapor=Chen_correlation(d,sat_steam.rho,Vaporization_Enthalpy,e,b,T_sat,x,Outer_diameter_tube)
    #debería ser la diferencia de temperatura y no el promedio de estas.
    T_pared=T_wall(h_correlacion_vapor,h_correlacion_agua,T_prom_liq,T_prom_vapor,A_tubo_ext,A_tubo_int)
    dif_t_wall=abs(x-T_pared)
    x=T_pared
    print'valor iteración ',n,'es: ',x,'[K]'
    n=n+1

print 'La temperatura de pared de la tuberia, obtenida por iteraciones es:',x-273.15,'[°C]'
print 'Coeficiente de transferencia calor convectivo por vapor: ',h_correlacion_vapor,'[W/m^2 K]'
print 'Coeficiente de transferencia calor convectivo por liquido: ',h_correlacion_agua,'[W/m^2 K]'
print 'Q_abs: ', -h_correlacion_agua*(T_prom_liq-x)*A_tubo_int, '[W]'
print 'Q_ced: ', -h_correlacion_vapor*(T_prom_vapor-x)*A_tubo_ext, '[W]'
Q_abs=h_correlacion_vapor*(T_prom_vapor-x)*A_tubo_ext
#CALCULO COEFICIENTE GLOBAL DE TRANSFERENCIA DE CALOR LIMPIO (LA FUNCION ES 1/U_c)
U_c=(Coef_calor_limpio(Outer_diameter_tube,Inner_diameter_tube,h_correlacion_agua,h_correlacion_vapor))**(-1.)
print 'el coeficiente global de transferencia de calor limpio es:', U_c,'[W/m^2 K]'
    
#CALCULO COEFICIENTE GLOBAL DE TRANSFERENCIA DE CALOR SUCIO (LA FUNCION ES 1/U_d)
DMLT=DMLT(T_ea,T_sa,T_ev,T_sv)
print'el valor de DMLT es:', DMLT, '[°C]'
U_d=pow(Coef_calor_sucio(A_tubo_int,DMLT,Q_abs),-1.)
print 'el coeficiente global de transferencia de calor sucio es:', U_d,'[W/m^2 K]'

#CALCULO FACTOR DE INCRUSTACION
R_d=pow(U_d,-1.)+pow(U_c,-1.)
print 'el coeficiente del factor de incrustación es:',R_d,'[m^2 °C/W]'
e=time.time()
print '***************************************************************'
print 'Tiempo de calculo: ',e-s,'[s]'