import mat_properties as prop
import numpy as n
import pandas as pd


def mixing_gases_molar(stream1, stream2, stream3, working_table):
    fractions1 = working_table.loc[stream1, "N2":]
    fractions2 = working_table.loc[stream2, "N2":]
    [gas1_T, gas1_P, gas1_H, gas1_G] = working_table.loc[stream1, "T":"G"]
    [gas2_T, gas2_P, gas2_H, gas2_G] = working_table.loc[stream2, "T":"G"]
    molar_mases = {
        "N2": 0.0280134,
        "O2": 0.03199806,
        "CO2": 0.0440095,
        "H2O": 0.01801528,
        "Ar": 0.039948,
    }
    components = list(working_table.columns)[4:]
    Sr_mol_mass1 = sum(map(lambda x1, x2: molar_mases[x1] * x2, components, fractions1))
    Sr_mol_mass2 = sum(map(lambda x1, x2: molar_mases[x1] * x2, components, fractions2))
    Molar_flow1 = gas1_G / Sr_mol_mass1
    Molar_flow2 = gas2_G / Sr_mol_mass2
    Molar_flow_components = list(
        map(lambda x1, x2: x1 * Molar_flow1 + x2 * Molar_flow2, fractions1, fractions2)
    )
    fractions3 = list(map(lambda x: x / sum(Molar_flow_components),Molar_flow_components))
    gas3_G = gas1_G + gas2_G
    gas3_H = (gas1_H * gas1_G + gas2_H * gas2_G) / gas3_G
    
    gasmix = "Nitrogen*Oxygen*CO2*Water*Argon"
    RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")

    gas_s_prop = prop.Materials_prop(
        gasmix,
        fractions3,
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
    )

    gas3_P = min(gas2_P, gas1_P)
    gas3_T = gas_s_prop.p_h(gas3_P, gas3_H)["T"]
    gas3_properties = [gas3_T, gas3_P, gas3_H, gas3_G] + fractions3
    working_table.loc[stream3, :] = gas3_properties
    gas3_properties_return = dict(zip(working_table.columns, gas3_properties))
    return gas3_properties_return

class steam_transformer:
    def __init__(self,**kwargs):
        self.stream11 = kwargs['stream11']
        self.water = kwargs['water']
        self.water_streams = kwargs['water_streams']
        self.heaters = kwargs['heaters']
        self.Pdr1 = kwargs['Pdr1']
        self.Pdr2 = kwargs['Pdr2']
        self.P2 = kwargs['P2']
        self.dT = kwargs['dT']
        self.dTmin = kwargs['dTmin']
        self.Tdec = kwargs['Tdec']
        
    def calc(self):
        P11 = self.water_streams.at[self.stream11, 'P']
        T11 = self.water_streams.at[self.stream11, 'T']
        H11 = self.water_streams.at[self.stream11, 'H']
        G1  = self.water_streams.at[self.stream11, 'G']

        Q12 = 0
        P12 = P11
        
        H12 = self.water.p_q(P12, Q12)['h']
        T12 = self.water.p_q(P12, Q12)['T']

        P13 = self.Pdr1
        H13 = H12
                
        T13 = self.water.p_h(P13, H13)['T']
        Q13 = self.water.p_h(P13, H13)['Q']
        
        P24 = self.P2
        T24 = T11-self.dT
        H24 = self.water.p_t(P24, T24)['h']

        T23 = T13-self.dTmin
        P23 = P24
        H23 = self.water.p_t(P23, T23)['h']

        G2 = G1 * ((H11-H12)/(H24-H23))
                
        Q14 = 0
        P14 = P13
        H14 = self.water.p_q(P14, Q14)['h']
        T14 = self.water.p_q(P14, Q14)['T']

        P15 = self.Pdr2
        H15 = H14
        T15 = self.water.p_h(P15, H15)['T']
        Q15 = self.water.p_h(P15, H15)['Q']

        P22 = self.P2
        H22 = H23 - (G1*(H13-H14)/G2)
        T22 = self.water.p_h(P22, H22)['T']

        T21 = 15
        P21 = self.P2
        H21 = self.water.p_t(P21, T21)['h']

        H16 = H15 - G2*(H22-H21)/G1
        P16 = P15
        T16 = self.water.p_h(P16, H16)['T']

        Qtrans = G1*(H11-H16)
        
        
        T17 = 80
        P17 = P16
        H17 = self.water.p_t(P17, T17)['h']
        Qcool80 = G1*(H16-H17)
        self.heaters.loc['Strans','Qw'] = Qtrans
        self.heaters.loc['Strans_cool','Qw'] = Qcool80

        return {'H11': H11,'H12': H12,'H13': H13,'H14': H14,'H15': H15,'H16': H16,
               'H21': H21,'H22': H22,'H23': H23,'H24': H24, 'G1':G1, 'G2':G2,'Q':Qtrans, 'P2':self.P2, 'P17':P17,
               'T11': T11,'T12': T12,'T13': T13,'T14': T14,'T15': T15,'T16': T16, 'T17':T17,'H17':H17, 
               'T21': T21,'T22': T22,'T23': T23,'T24': T24, 'Qcool80':Qcool80}

class reformer:
    def __init__(self, **kwargs):
        self.stream11 = kwargs['stream11']
        self.water = kwargs['water']
        self.water_streams = kwargs['water_streams']
        self.heaters = kwargs['heaters']
        self.Methane = kwargs['Methane']
        self.waterMethane = kwargs['waterMethane']
        self.Tref = kwargs['Tref']
        self.Pref = kwargs['Pref']
        self.gas_KU = kwargs['gas_KU']
        self.T1gas = kwargs['T1gas']
        self.T2gas = kwargs['T2gas']
        
    def calc(self):
        Gsteam = self.water_streams.at[self.stream11,'G']
        Hsteam = self.water_streams.at[self.stream11,'H']
        Gmeth = (0.151140511695727/0.848859488304273)*Gsteam
        
        Gref = Gmeth + Gsteam
        Tmeth = 121.53615158005
        Hmeth = self.Methane.p_t(self.Pref, Tmeth)['h']
        H1r = (Gsteam*Hsteam + Gmeth*Hmeth) / Gref
        H2r = self.waterMethane.p_t(self.Pref, self.Tref)['h']
        Qdt = Gref*(H2r-H1r)
        Qreac = (53634.90756/44.6397313913235)*Gref
        Qref = Qdt+Qreac
        H1gas = 3191.2518095937544
        H2gas = self.gas_KU.p_t(0.1, self.T2gas)['h']
        
        Ggas = Qref/(H1gas-H2gas)
        
        Gair = Ggas*48.5966429877363/51.2672005145991
        Gch4 = Ggas*2.67055752686283/51.2672005145991
        Hsg = H2r + (Qreac/Gref)
                
        SGsost = ['N2','O2','CO2','Ar','H2O','CH4','H2','CO']
        SGfrac = [0,0,0.06356,0,0.53662,0.05004,0.32590,0.02388]
        # All_in_one=sum(SGfrac)
        # SGfrac= map(lambda x: x/All_in_one,SGfrac)
        SGfrac = dict(zip(SGsost,SGfrac))

        Gassost = ['N2','O2','CO2','H2O','Ar']
        Gasfrac = [0.710320591016015,0.00996710270335893,0.090538556815177,0.180531273012258,0.00864247645319178]
        Gasfrac = dict(zip(Gassost,Gasfrac))
        self.heaters.loc['Reformer','Qw'] = Qref
        return {'Q':Qref,'Ggas':Ggas,'Gair':Gair,'Gch4':Gch4, 'Tref':self.Tref,
                'Gref':Gref, 'Pref':self.Pref,'Hsg':Hsg,'SGfrac':SGfrac,'H2gas':H2gas,'Gasfrac':Gasfrac}

class PKM_cooler:
    def __init__(self, stream1, stream2, syngas_streams, heaters, Tout):
        self.stream1 = stream1
        self.stream2 = stream2
        self.syngas_streams = syngas_streams
        self.heaters = heaters
        self.Tout = Tout
    def calc(self):
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = (0,0,0.06356,0,0.53662,0.05004,0.32590,0.02388)
        Tin    = self.syngas_streams.at[self.stream1,'T']
        Pin    = self.syngas_streams.at[self.stream1,'P']
        G      = self.syngas_streams.at[self.stream1,'G']
        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
        SG = prop.Materials_prop(SGsost,SGfrac,prop.REFPROP_h_s,prop.REFPROP_p_t,prop.REFPROP_p_h,
                                     prop.REFPROP_p_s,prop.REFPROP_p_q,prop.REFPROP_t_q,prop.REFPROP_p_rho,prop.REFPROP_s_q,RP=RP)
        Hin  = SG.p_t(Pin,Tin)['h']
        Tout = self.Tout
        Hout = SG.p_t(Pin,Tout)['h']
        Qcooler = G*(Hin-Hout)
        self.syngas_streams.loc[self.stream2,'T':'G'] = [Tout, Pin, Hout, G]
        self.syngas_streams.loc['REF-COOL','H'] = Hin
        self.syngas_streams.loc[self.stream2,'N2':'CO'] = self.syngas_streams.loc[self.stream1,'N2':'CO']
        self.heaters.loc['Ref_cooler','Qw'] = Qcooler
        return {'Tout':Tout, 'Hout':Hout, 'Pout':Pin,'G':G, 'Qcooler':Qcooler}
    
class HTS:
    def __init__(self, stream1, stream2, syngas_streams, heaters, Tout):
        self.stream1 = stream1
        self.stream2 = stream2
        self.syngas_streams = syngas_streams
        self.heaters = heaters
        self.Tout = Tout
        
    def calc(self):
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = (0,0,0.06356,0,0.53662,0.05004,0.32590,0.02388)
        Tin    = self.syngas_streams.at[self.stream1,'T']
        Pin    = self.syngas_streams.at[self.stream1,'P']
        G      = self.syngas_streams.at[self.stream1,'G']
        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
        SG = prop.Materials_prop(SGsost,SGfrac,prop.REFPROP_h_s,prop.REFPROP_p_t,prop.REFPROP_p_h,
                                     prop.REFPROP_p_s,prop.REFPROP_p_q,prop.REFPROP_t_q,prop.REFPROP_p_rho,prop.REFPROP_s_q,RP=RP)
        Hin  = SG.p_t(Pin,Tin)['h']
#         Hdt = SG.p_t(Pin,self.Tout)['h']
#         Qdt = G*(Hin-Hdt)
#         Qreac = (2701.173347/44.6397313913235)*G
#         Qhts = Qdt+Qreac
        SGfracnew = (0,0,0.0864149892361543,0,0.513766606256586,0.0500421238434872,0.348747199710724,0.00102908095304842)
        SGnew = prop.Materials_prop(SGsost,SGfracnew,prop.REFPROP_h_s,prop.REFPROP_p_t,prop.REFPROP_p_h,
                                     prop.REFPROP_p_s,prop.REFPROP_p_q,prop.REFPROP_t_q,prop.REFPROP_p_rho,prop.REFPROP_s_q,RP=RP)
        Hout  = SGnew.p_t(Pin,self.Tout)['h']
        Qhts = G*(Hin-Hout)
        self.syngas_streams.loc[self.stream2,'T':'G'] = [self.Tout, Pin, Hout, G]
        self.syngas_streams.loc[self.stream2,'N2':'CO'] = list(SGfracnew)
        self.heaters.loc['Ref_HTS','Qw'] = Qhts
        return {'Qhts':Qhts, 'Tout': self.Tout, 'Pout':Pin, 'Hout':Hout, 'G':G}
    
    
class PKM_all:
    def calc(PKM_zaryad,gas_streams,syngas_streams,water_streams,water_streams0,heaters):
        
        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
        water = prop.Materials_prop(
            "water",
            [1.0, 0, 0, 0, 0],
            prop.REFPROP_h_s,
            prop.REFPROP_p_t,
            prop.REFPROP_p_h,
            prop.REFPROP_p_s,
            prop.REFPROP_p_q,
            prop.REFPROP_t_q,
            prop.REFPROP_p_rho,
            prop.REFPROP_s_q,
            RP=RP,
        )
        waterMethane = prop.Materials_prop(
            "Water*METHANE",
            (0.833372660622383, 0.166627339377617, 0, 0, 0),
            prop.REFPROP_h_s,
            prop.REFPROP_p_t,
            prop.REFPROP_p_h,
            prop.REFPROP_p_s,
            prop.REFPROP_p_q,
            prop.REFPROP_t_q,
            prop.REFPROP_p_rho,
            prop.REFPROP_s_q,
            RP=RP,
        )
        Methane = prop.Materials_prop(
            "METHANE",
            [1.0, 0, 0, 0, 0],
            prop.REFPROP_h_s,
            prop.REFPROP_p_t,
            prop.REFPROP_p_h,
            prop.REFPROP_p_s,
            prop.REFPROP_p_q,
            prop.REFPROP_t_q,
            prop.REFPROP_p_rho,
            prop.REFPROP_s_q,
            RP=RP,
        )
        gas_KU_PKM = prop.Materials_prop(
            "Nitrogen*Oxygen*CO2*Water*Argon",
            (0.710320591016015,0.00996710270335893,0.090538556815177,0.180531273012258,0.00864247645319178),
            prop.REFPROP_h_s,
            prop.REFPROP_p_t,
            prop.REFPROP_p_h,
            prop.REFPROP_p_s,
            prop.REFPROP_p_q,
            prop.REFPROP_t_q,
            prop.REFPROP_p_rho,
            prop.REFPROP_s_q,
            RP=RP,
        )
        
        
        
        
        