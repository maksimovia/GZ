from CoolProp.CoolProp import PropsSI
import json, CoolProp.CoolProp as CP
CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, 'C:\\Program Files (x86)\\REFPROP\\')

#Gas_composition_0={'N2':78.03,'O2':12.37,'CO2':3.01,'H2O':5.94,'Ar':0.65}
Gas_method='REFPROP'
Water_method='HEOS'


class Composition_properties(object):
    pass


def Composition_string(Gas_composition):
    Composition_properties = Gas_method + '::' + '&'.join(
        map(lambda x: x + '[' + str(Gas_composition[x] / 100) + ']', Gas_composition))
    return Composition_properties

fluid_water = Water_method+'::Water'

#Gas properties

def gas_H_PT(P,T,Gas_composition):
    return PropsSI('H','T',T+273.15,'P',P*1000000,Composition_string(Gas_composition))/1000

def gas_T_HP(H,P,Gas_composition):
    return float(PropsSI('T','H',H*1000,'P',P*1000000,Composition_string(Gas_composition)))-273.15

def gas_DVIS_PT(P,T,Gas_composition):
    Gas_composition_local=dict(Gas_composition)
    if Gas_composition_local['H2O']>5:
        Gas_composition_local['CO2']+=(Gas_composition_local['H2O']-5)
        Gas_composition_local['H2O']=5
    return PropsSI('V','T',T+273.15,'P',P*1000000,Composition_string(Gas_composition_local))

def gas_KINVIS_PT(P,T,Gas_composition):
    Gas_composition_local=dict(Gas_composition)
    if Gas_composition_local['H2O']>5:
        Gas_composition_local['CO2']+=(Gas_composition_local['H2O']-5)
        Gas_composition_local['H2O']=5
    return PropsSI('V','T',T+273.15,'P',P*1000000,Composition_string(Gas_composition_local))/PropsSI('D','T',T+273.15,'P',P*1000000,Composition_string(Gas_composition_local))

def gas_Prandtl_PT(P,T,Gas_composition):
    Gas_composition_local=dict(Gas_composition)
    if Gas_composition_local['H2O']>5:
        Gas_composition_local['CO2']+=(Gas_composition_local['H2O']-5)
        Gas_composition_local['H2O']=5
    return PropsSI('Prandtl','T',T+273.15,'P',P*1000000,Composition_string(Gas_composition_local))

def gas_D_PT(P,T,Gas_composition):

    return PropsSI('D','T',T+273.15,'P',P*1000000,Composition_string(Gas_composition))

def gas_L_PT(P,T,Gas_composition):
    Gas_composition_local=dict(Gas_composition)
    if Gas_composition_local['H2O']>5:
        Gas_composition_local['CO2']+=(Gas_composition_local['H2O']-5)
        Gas_composition_local['H2O']=5
    return PropsSI('L','T',T+273.15,'P',P*1000000,Composition_string(Gas_composition_local))

#print(gas_Prandtl_PT(0.1013,200,Gas_composition_0))

#Water properties

def water_H_PT(P,T):
    return PropsSI('H','T',T+273.15,'P',P*1000000,fluid_water)/1000

def water_T_HP(H,P):
    return float(PropsSI('T','H',H*1000,'P',P*1000000,fluid_water))-273.15

def water_HSS_T(T):
    return PropsSI('H','T',T+273.15,'Q',1,fluid_water)/1000

def water_HSS_P(P):
    return PropsSI('H','P',P*1000000,'Q',1,fluid_water)/1000

def water_TSS_P(P):
    return PropsSI('T','P',P*1000000,'Q',0,fluid_water)-273.15

def water_HSW_P(P):
    return PropsSI('H','P',P*1000000,'Q',0,fluid_water)/1000

def water_HSW_T(T):
    return PropsSI('H','T',T+273.15,'Q',0,fluid_water)/1000

def water_D_PT(P,T):
    return PropsSI('D','T',T+273.15,'P',P*1000000,fluid_water)

def water_DSS_P (P):
    return PropsSI('D','Q',1,'P',P*1000000,fluid_water)

def water_Q_HP(H,P):
    return float(PropsSI('Q','H',H*1000,'P',P*1000000,fluid_water))

def water_DSW_P(P):
    return PropsSI('D','Q',0,'P',P*1000000,fluid_water)

def water_DVIS_PT(P,T):
    return PropsSI('V','T',T+273.15,'P',P*1000000,fluid_water)

def water_KINVIS_PT(P,T):
    return PropsSI('V','T',T+273.15,'P',P*1000000,fluid_water)/PropsSI('D','T',T+273.15,'P',P*1000000,fluid_water)

def water_Prandtl_PT(P,T):
    return PropsSI('PRANDTL','T',T+273.15,'P',P*1000000,fluid_water)

def water_THCOND_PT(P,T):
    return PropsSI('CONDUCTIVITY','T',T+273.15,'P',P*1000000,fluid_water)

#print(water_THCOND_PT(1,20))