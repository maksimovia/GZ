import mat_properties as prop
import numpy as n
import pandas as pd
from scipy.optimize import root

RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")


def mixing_gases_molar(stream1, stream2, stream3, working_table):
    fractions1 = working_table.loc[stream1, "N2":]
    fractions2 = working_table.loc[stream2, "N2":]
    [gas1_T, gas1_P, gas1_H, gas1_G] = working_table.loc[stream1, "T":"G"]
    [gas2_T, gas2_P, gas2_H, gas2_G] = working_table.loc[stream2, "T":"G"]
    molar_mases = {"CO2": 0.04401,
                       "H2O": 0.01801528,
                       "Ar":  0.039948,
                       "CH4": 0.01604,
                       "H2":  0.002,
                       "CO":  0.02801,
                       "N2": 0.0280134,
                       "O2": 0.015999*2}
    components = list(working_table.columns)[4:]
    Sr_mol_mass1 = sum(
        map(lambda x1, x2: molar_mases[x1] * x2, components, fractions1))
    Sr_mol_mass2 = sum(
        map(lambda x1, x2: molar_mases[x1] * x2, components, fractions2))
    Molar_flow1 = gas1_G / Sr_mol_mass1
    Molar_flow2 = gas2_G / Sr_mol_mass2
    Molar_flow_components = list(
        map(lambda x1, x2: x1 * Molar_flow1 + x2 * Molar_flow2, fractions1, fractions2))
    fractions3 = list(
        map(lambda x: x / sum(Molar_flow_components), Molar_flow_components))
    gas3_G = gas1_G + gas2_G
    gas3_H = (gas1_H * gas1_G + gas2_H * gas2_G) / gas3_G

    gas_s_prop = prop.Materials_prop("Nitrogen*Oxygen*CO2*Water*Argon",
                                     fractions3,
                                     prop.REFPROP_h_s,
                                     prop.REFPROP_p_t,
                                     prop.REFPROP_p_h,
                                     prop.REFPROP_p_s,
                                     prop.REFPROP_p_q,
                                     prop.REFPROP_t_q,
                                     prop.REFPROP_p_rho,
                                     prop.REFPROP_s_q,
                                     RP=RP)

    gas3_P = min(gas2_P, gas1_P)
    gas3_T = gas_s_prop.p_h(gas3_P, gas3_H)["T"]
    gas3_properties = [gas3_T, gas3_P, gas3_H, gas3_G] + fractions3
    working_table.loc[stream3, :] = gas3_properties
    gas3_properties_return = dict(zip(working_table.columns, gas3_properties))
    return gas3_properties_return


class steam_transformer:
    def __init__(self, **kwargs):
        self.stream11 = kwargs['stream11']
        self.water_streams = kwargs['water_streams']
        self.water = kwargs['water']
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
        G1 = self.water_streams.at[self.stream11, 'G']
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
#         T24 = 500
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

        T17 = self.water_streams.loc['SWIN', 'T']+10
        P17 = P16
        H17 = self.water.p_t(P17, T17)['h']

        Qcool80 = max(G1*(H16-H17), 0)
        self.heaters.loc['Strans', 'Qw'] = Qtrans
        self.heaters.loc['Strans_cool', 'Qw'] = Qcool80
        return {'H11': H11, 'H16': H16, 'H21': H21, 'H24': H24, 'G1': G1, 'G2': G2, 'Q': Qtrans, 'P2': self.P2, 'P17': P17,
                'T11': T11, 'T17': T17, 'H17': H17, 'T21': T21, 'T24': T24, 'Qcool80': Qcool80}


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
        Gsteam = self.water_streams.at[self.stream11, 'G']
        Hsteam = self.water_streams.at[self.stream11, 'H']
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

        SGsost = ['N2', 'O2', 'CO2', 'Ar', 'H2O', 'CH4', 'H2', 'CO']
        SGfrac = [0, 0, 0.06356, 0, 0.53662, 0.05004, 0.32590, 0.02388]
        # All_in_one=sum(SGfrac)
        # SGfrac= map(lambda x: x/All_in_one,SGfrac)
        SGfrac = dict(zip(SGsost, SGfrac))

        Gassost = ['N2', 'O2', 'CO2', 'H2O', 'Ar']
        Gasfrac = [0.710320591016015, 0.00996710270335893,
                   0.090538556815177, 0.180531273012258, 0.00864247645319178]
        Gasfrac = dict(zip(Gassost, Gasfrac))
        self.heaters.loc['Reformer', 'Qw'] = Qref
        return {'Q': Qref, 'Ggas': Ggas, 'Gair': Gair, 'Gch4': Gch4, 'Tref': self.Tref,
                'Gref': Gref, 'Pref': self.Pref, 'Hsg': Hsg, 'SGfrac': SGfrac, 'H2gas': H2gas, 'Gasfrac': Gasfrac}


class PKM_cooler:
    def __init__(self, stream1, stream2, syngas_streams, heaters, Tout):
        self.stream1 = stream1
        self.stream2 = stream2
        self.syngas_streams = syngas_streams
        self.heaters = heaters
        self.Tout = Tout

    def calc(self):
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = (0, 0, 0.06356, 0, 0.53662, 0.05004, 0.32590, 0.02388)
        Tin = self.syngas_streams.at[self.stream1, 'T']
        Pin = self.syngas_streams.at[self.stream1, 'P']
        G = self.syngas_streams.at[self.stream1, 'G']
        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
        SG = prop.Materials_prop(SGsost,
                                 SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)
        Hin = SG.p_t(Pin, Tin)['h']
        Tout = self.Tout
        Hout = SG.p_t(Pin, Tout)['h']
        Qcooler = G*(Hin-Hout)
        self.syngas_streams.loc[self.stream2, 'T':'G'] = [Tout, Pin, Hout, G]
        self.syngas_streams.loc['REF-COOL', 'H'] = Hin
        self.syngas_streams.loc[self.stream2,
                                'N2':'CO'] = self.syngas_streams.loc[self.stream1, 'N2':'CO']
        self.heaters.loc['Ref_cooler', 'Qw'] = Qcooler
        return {'Tout': Tout, 'Hout': Hout, 'Pout': Pin, 'G': G, 'Qcooler': Qcooler}


class HTS:
    def __init__(self, stream1, stream2, syngas_streams, heaters, Tout):
        self.stream1 = stream1
        self.stream2 = stream2
        self.syngas_streams = syngas_streams
        self.heaters = heaters
        self.Tout = Tout

    def calc(self):
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = (0, 0, 0.06356, 0, 0.53662, 0.05004, 0.32590, 0.02388)
        Tin = self.syngas_streams.at[self.stream1, 'T']
        Pin = self.syngas_streams.at[self.stream1, 'P']
        G = self.syngas_streams.at[self.stream1, 'G']
        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
        SG = prop.Materials_prop(SGsost,
                                 SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)
        Hin = SG.p_t(Pin, Tin)['h']
        SGfracnew = (0, 0, 0.0864149892361543, 0, 0.513766606256586,
                     0.0500421238434872, 0.348747199710724, 0.00102908095304842)
        SGnew = prop.Materials_prop(SGsost,
                                    SGfracnew,
                                    prop.REFPROP_h_s,
                                    prop.REFPROP_p_t,
                                    prop.REFPROP_p_h,
                                    prop.REFPROP_p_s,
                                    prop.REFPROP_p_q,
                                    prop.REFPROP_t_q,
                                    prop.REFPROP_p_rho,
                                    prop.REFPROP_s_q,
                                    RP=RP)
        Hout = SGnew.p_t(Pin, self.Tout)['h']
        Qhts = G*(Hin-Hout)
        self.syngas_streams.loc[self.stream2, 'T':'G'] = [
            self.Tout, Pin, Hout, G]
        self.syngas_streams.loc[self.stream2, 'N2':'CO'] = list(SGfracnew)
        self.heaters.at['Ref_HTS','Qw'] = Qhts
        return {'Qhts': Qhts, 'Tout': self.Tout, 'Pout': Pin, 'Hout': Hout, 'G': G}


class HTS_cooler:
    def calc(stream1, stream2, syngas_streams, heaters, Tout):
        Hin = syngas_streams.loc[stream1, 'H']
        Pin = syngas_streams.loc[stream1, 'P']
        G = syngas_streams.loc[stream1, 'G']
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = syngas_streams.loc[stream1, 'N2':'CO']
        SG = prop.Materials_prop(SGsost,
                                 SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)
        Hout = SG.p_t(Pin, Tout)['h']
        syngas_streams.loc['HTSCOOL-Separ', 'T':'G'] = [Tout, Pin, Hout, G]
        Q = G*(Hin-Hout)
        qual = SG.p_t(Pin, Tout)['Q']
        heaters.loc['HTS_cooler', 'Qw'] = Q
        syngas_streams.loc[stream2, 'N2':'CO'] = list(SGfrac)
        return Q


class separator:
    def calc(stream1, stream2, syngas_streams, heaters):
        Hin = syngas_streams.loc[stream1, 'H']
        Pin = syngas_streams.loc[stream1, 'P']
        Tin = syngas_streams.loc[stream1, 'T']
        G = syngas_streams.loc[stream1, 'G']
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = syngas_streams.loc[stream1, 'N2':'CO']
        SG = prop.Materials_prop(SGsost,
                                 SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)
        qual = SG.p_t(Pin, Tin)['Q']
        Mass_liq_f = 1-qual
        Gliq = G*Mass_liq_f
        Gout = G - Gliq

        molar_mases = {"CO2": 0.04401,
                       "H2O": 0.01801528,
                       "Ar":  0.039948,
                       "CH4": 0.01604,
                       "H2":  0.002,
                       "CO":  0.02801,
                       "N2": 0.0280134,
                       "O2": 0.015999*2}

        Sr_mol_mass = sum(map(
            lambda x1, x2: molar_mases[x1] * x2, list(syngas_streams.columns)[4:], SGfrac))
        Mole_flow = G/Sr_mol_mass
        Mole_flow_CO2 = Mole_flow*SGfrac["CO2"]
        Mole_flow_H2O = Mole_flow*SGfrac["H2O"]
        Mole_flow_Ar = Mole_flow*SGfrac["Ar"]
        Mole_flow_CH4 = Mole_flow*SGfrac["CH4"]
        Mole_flow_H2 = Mole_flow*SGfrac["H2"]
        Mole_flow_CO = Mole_flow*SGfrac["CO"]
        Mole_flow_N2 = Mole_flow*SGfrac["N2"]
        Mole_flow_O2 = Mole_flow*SGfrac["O2"]
        Mole_sep = Gliq/molar_mases['H2O']
        Mole_flow_H2O = max(Mole_flow_H2O-Mole_sep, 0)
        Mole_flow_out = Mole_flow_CO2+Mole_flow_Ar + \
            (Mole_flow_H2O)+Mole_flow_CH4+Mole_flow_H2 + \
            Mole_flow_CO+Mole_flow_N2+Mole_flow_O2
        SGfrac_new = [Mole_flow_N2/Mole_flow_out, Mole_flow_O2/Mole_flow_out, Mole_flow_CO2/Mole_flow_out, Mole_flow_Ar/Mole_flow_out,
                      (Mole_flow_H2O)/Mole_flow_out, Mole_flow_CH4/Mole_flow_out, Mole_flow_H2/Mole_flow_out, Mole_flow_CO/Mole_flow_out]
        Sr_mol_mass_out = sum(map(
            lambda x1, x2: molar_mases[x1] * x2,  list(syngas_streams.columns)[4:], SGfrac_new))
        Mass_flow_out = Mole_flow_out*Sr_mol_mass_out
        SGnew = prop.Materials_prop(SGsost,
                                    SGfrac_new,
                                    prop.REFPROP_h_s,
                                    prop.REFPROP_p_t,
                                    prop.REFPROP_p_h,
                                    prop.REFPROP_p_s,
                                    prop.REFPROP_p_q,
                                    prop.REFPROP_t_q,
                                    prop.REFPROP_p_rho,
                                    prop.REFPROP_s_q,
                                    RP=RP)
        Hout = SGnew.p_t(Pin, Tin)['h']
        syngas_streams.loc[stream2, 'N2':'CO'] = SGfrac_new
        syngas_streams.loc[stream2, 'T':'G'] = [Tin, Pin, Hout, Mass_flow_out]


class PKM_all:
    def calc(PKM_zaryad, gas_streams, syngas_streams, water_streams, water_streams0, heaters):

        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
        water = prop.Materials_prop("water",
                                    [1.0, 0, 0, 0, 0],
                                    prop.REFPROP_h_s,
                                    prop.REFPROP_p_t,
                                    prop.REFPROP_p_h,
                                    prop.REFPROP_p_s,
                                    prop.REFPROP_p_q,
                                    prop.REFPROP_t_q,
                                    prop.REFPROP_p_rho,
                                    prop.REFPROP_s_q,
                                    RP=RP)
        waterMethane = prop.Materials_prop("Water*METHANE",
                                           (0.833372660622383,
                                            0.166627339377617, 0, 0, 0),
                                           prop.REFPROP_h_s,
                                           prop.REFPROP_p_t,
                                           prop.REFPROP_p_h,
                                           prop.REFPROP_p_s,
                                           prop.REFPROP_p_q,
                                           prop.REFPROP_t_q,
                                           prop.REFPROP_p_rho,
                                           prop.REFPROP_s_q,
                                           RP=RP)
        Methane = prop.Materials_prop("METHANE",
                                      [1.0, 0, 0, 0, 0],
                                      prop.REFPROP_h_s,
                                      prop.REFPROP_p_t,
                                      prop.REFPROP_p_h,
                                      prop.REFPROP_p_s,
                                      prop.REFPROP_p_q,
                                      prop.REFPROP_t_q,
                                      prop.REFPROP_p_rho,
                                      prop.REFPROP_s_q,
                                      RP=RP)
        gas_KU_PKM = prop.Materials_prop("Nitrogen*Oxygen*CO2*Water*Argon",
                                         (0.710320591016015, 0.00996710270335893, 0.090538556815177,
                                          0.180531273012258, 0.00864247645319178),
                                         prop.REFPROP_h_s,
                                         prop.REFPROP_p_t,
                                         prop.REFPROP_p_h,
                                         prop.REFPROP_p_s,
                                         prop.REFPROP_p_q,
                                         prop.REFPROP_t_q,
                                         prop.REFPROP_p_rho,
                                         prop.REFPROP_s_q,
                                         RP=RP)
        if PKM_zaryad == True:
            water_streams.loc["DROSVD-TURBVD", "G"] = water_streams.loc["PEVD-DROSVD",
                                                                        "G"]-water_streams.loc["DROSVD-ST", "G"]
            water_streams.loc['DROSVD-ST',
                              'P'] = water_streams.loc['PEVD-DROSVD', 'P']
            # паротрансформатор
            ST = steam_transformer(stream11="DROSVD-ST",
                                   water=water,
                                   water_streams=water_streams,
                                   heaters=heaters,
                                   Pdr1=2, Pdr2=0.8, P2=2, dT=15, dTmin=5, Tdec=10)
            steam_trans = ST.calc()

            # Ввод в табл выходов из паротрансформатора
            water_streams.loc["ST-GPK", "T":"G"] = [steam_trans["T17"],
                                                    steam_trans["P17"], steam_trans["H17"], steam_trans["G1"]]
            water_streams.loc["ST-PKM", "T":"G"] = [steam_trans["T24"],
                                                    steam_trans["P2"], steam_trans["H24"], steam_trans["G2"]]

            # реформер
            # from PKM import reformer
            ref = reformer(stream11="ST-PKM",
                           water=water,
                           gas_KU=gas_KU_PKM,
                           Methane=Methane,
                           waterMethane=waterMethane,
                           water_streams=water_streams,
                           heaters=heaters,
                           Tref=700, Pref=2, T1gas=1968.58395330148, T2gas=800)
            reform = ref.calc()

            # Синтезгаз из рефрмера
            syngas_streams.loc["REF-COOL", "T":"G"] = [reform["Tref"],
                                                       reform["Pref"], reform["Hsg"], reform["Gref"]]
            syngas_streams.loc["REF-COOL",
                               "N2":"CO"] = list(reform["SGfrac"].values())

            # Газы реформера
            gas_streams.loc["AIR-REF", "T":"G"] = [15,
                                                   0.1, 288.39, reform["Gair"]]
            gas_streams.loc["CH4-REF", "T":"G"] = [15,
                                                   0.7, 881.50, reform["Gch4"]]
            gas_streams.loc["REF-SMESH", "T":"G"] = [800,
                                                     0.1, reform["H2gas"], reform["Ggas"]]
            gas_streams.loc["REF-SMESH",
                            "N2":"Ar"] = list(reform["Gasfrac"].values())

            # Смешение
            gas_streams.loc["GTU-PEVD", "G"] = (
                gas_streams.at["REF-SMESH", "G"] + gas_streams.at["GTU-KU", "G"])
            gas_streams.loc["GTU-PEVD", "H"] = (gas_streams.at["REF-SMESH", "G"] * gas_streams.at["REF-SMESH", "H"] +
                                                gas_streams.at["GTU-KU", "G"] * gas_streams.at["GTU-KU", "H"]) / gas_streams.loc["GTU-PEVD", "G"]
            gas_streams.loc["GTU-PEVD", "P"] = 0.1
            # from PKM import mixing_gases_molar
            mixing_gases_molar("GTU-KU", "REF-SMESH", "GTU-PEVD", gas_streams)
            for stream in gas_streams.index[4:10]:
                gas_streams.loc[stream,
                                "N2":"Ar"] = gas_streams.loc["GTU-PEVD", "N2":"Ar"]
            #Cooler + HTS
            cool = PKM_cooler('REF-COOL', 'COOL-HTS',
                              syngas_streams, heaters, 450).calc()
            hts = HTS('COOL-HTS', 'HTS-HTSCOOL',
                      syngas_streams, heaters, 275).calc()
            # +HTS_cooler
            cool = HTS_cooler.calc(
                'HTS-HTSCOOL', 'HTSCOOL-Separ', syngas_streams, heaters, 100)
            # from PKM import separator
            sep = separator.calc(
                'HTSCOOL-Separ', 'Separ-SGaccum', syngas_streams, heaters)
            Qref_all = heaters.loc["Ref_HTS", "Qw"] + heaters.loc["Ref_cooler",
                                                                  "Qw"] + heaters.loc["Strans_cool", "Qw"]+heaters.loc["HTS_cooler", "Qw"]
            
            
        else:
            water_streams.loc["DROSVD-TURBVD",
                              "G"] = water_streams.loc["PEVD-DROSVD", "G"]
            gas_streams.loc["GTU-PEVD",
                            "T":"Ar"] = gas_streams.loc["GTU-KU", "T":"Ar"]
            water_streams.loc["ST-GPK", "T":"G"] = [80, 2, 336.57, 0]
            heaters.loc["Strans":"HTS_cooler", "Qw"] = 0
            Qref_all = 0
        heaters.loc["Ref_all", "Qw"] = Qref_all
        return {'Qref_all': Qref_all}


class syngas_GTU:
    def calc(water_streams,syngas_streams, P2, KPDcomp, KPDturb, heaters, electric,accumulation):
        Tsg = syngas_streams.loc["Sepout-COMB", "T"]
        Psg = syngas_streams.loc["Sepout-COMB", "P"]
        Hsg = syngas_streams.loc["Sepout-COMB", "H"]
        Gsg = syngas_streams.loc["Sepout-COMB", "G"]
        
        if int(Gsg) == 0:
            print('Раход на СГтурбину турбину =0, всё ушло на осн ГТУ')
            syngas_streams.loc['AIR-COMP', 'T':'G'] = [0, 0, 0, 0]
            syngas_streams.loc['COMP-COMB', 'T':'G'] = [0, 0, 0, 0]
            syngas_streams.loc['COMB-TURB', 'T':'G'] = [0, 0, 0, 0]
            syngas_streams.loc['TURB-COOL','T':'G'] = [0, 0, 0, 0]
            syngas_streams.loc['COOL-EX', 'T':'G'] = [0, 0, 0, 0]
            syngas_streams.loc['COMB-TURB', 'N2':'CO'] = 0
            syngas_streams.loc['TURB-COOL', 'N2':'CO'] = 0
            syngas_streams.loc['COOL-EX', 'N2':'CO'] = 0
            heaters.loc['SGgtu_Qsw', 'Qw'] = 0
            electric.loc['SGgtu_comp', 'Ni'] = 0
            electric.loc['SGgtu_turb', 'Ni'] = 0
            return 0        
        RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
        air = prop.Materials_prop("Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO",
                                  (0.78, 0.21, 0.01, 0, 0, 0, 0, 0),
                                  prop.REFPROP_h_s,
                                  prop.REFPROP_p_t,
                                  prop.REFPROP_p_h,
                                  prop.REFPROP_p_s,
                                  prop.REFPROP_p_q,
                                  prop.REFPROP_t_q,
                                  prop.REFPROP_p_rho,
                                  prop.REFPROP_s_q,
                                  RP=RP)
        # Компр
        Hair1 = air.p_t(0.1, 15)['h']
        Sair1 = air.p_t(0.1, 15)['s']
        Hair2t = air.p_s(P2, Sair1)['h']
        Hair2 = Hair1 + (Hair2t-Hair1)/KPDcomp
        Tair2 = air.p_h(P2, Hair2)['T']

        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = list(syngas_streams.loc['Sepout-COMB', 'N2':'CO'])
        SG = prop.Materials_prop(SGsost,
                                 SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)
        # ПЕРЕДЕЛАТЬ СГОРАНИЕ
        molar_mases = {"CO2": 0.04401,
                       "H2O": 0.01801528,
                       "Ar":  0.039948,
                       "CH4": 0.01604,
                       "H2":  0.002,
                       "CO":  0.02801,
                       "N2": 0.0280134,
                       "O2": 0.015999*2}

        fractions_syngas = SGfrac
        components = list(syngas_streams.columns)[4:]

        Sr_mol_mass_syngas = sum(
            map(lambda x1, x2: molar_mases[x1] * x2, components, fractions_syngas))
        Qnr_molar_elements = {"CO2": 0,
                              "H2O": 0,
                              "Ar":  0,
                              "CH4": 802.34,
                              "H2":  244,
                              "CO":  283.24,
                              "N2": 0,
                              "O2": 0}
        Qnr_molar = sum(
            map(lambda x1, x2: Qnr_molar_elements[x1]*x2, components, fractions_syngas))
        Qnr_mass = Qnr_molar/Sr_mol_mass_syngas
        Molar_flow_syngas = Gsg/Sr_mol_mass_syngas
        Molar_components_flows_syngas = list(
            map(lambda x: x*Molar_flow_syngas, fractions_syngas))
        #Writing LHV in accumulation
        accumulation.at["PKM", "Qnr"] = Qnr_mass/1000
        # сост ух газов
        air_composition = ["N2", "O2", "Ar"]
        air_fraction_composition = [0.78, 0.21, 0.01]
        combsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        combfrac=syngas_streams.loc['COMB-TURB', 'N2':'CO'] 
        if combfrac[0]==0:
            combfrac = (0.5, 0.2,0.2,
                    0, 0.1, 0, 0, 0)  #первое приближение

        max_iterations_gas_turbine = 20
        Gair_tol = 10**-4
        Tcomb = 750
        Gair_it = [0]
        Sr_mol_mass_air = sum(map(lambda x1, x2: x1*molar_mases[x2], air_fraction_composition, air_composition))

        Molar_composition_exhaust=combfrac
        

        for i in range(max_iterations_gas_turbine):
            
            COMB = prop.Materials_prop(combsost,
                                       Molar_composition_exhaust,
                                       prop.REFPROP_h_s,
                                       prop.REFPROP_p_t,
                                       prop.REFPROP_p_h,
                                       prop.REFPROP_p_s,
                                       prop.REFPROP_p_q,
                                       prop.REFPROP_t_q,
                                       prop.REFPROP_p_rho,
                                       prop.REFPROP_s_q,
                                       RP=RP)
            Hcomb = COMB.p_t(P2, Tcomb)['h']
            Gair = -Gsg*((Hsg-Hcomb+Qnr_mass)/(Hair2-Hcomb))
            Gcomb = Gair+Gsg
            Molar_flow_air = Gair/Sr_mol_mass_air
            Molar_components_flows_air = list(
                map(lambda x: x*Molar_flow_air, air_fraction_composition))
            
            # Hydrogen combustion
            Molar_flow_H2O_H2 = Molar_components_flows_syngas[components.index(
                "H2")]
            Molar_flow_O2used_H2 = Molar_components_flows_syngas[components.index(
                "H2")]/2

            # Methane combustion
            Molar_flow_H2O_CH4 = Molar_components_flows_syngas[components.index(
                "CH4")]*2
            Molar_flow_CO2_CH4 = Molar_components_flows_syngas[components.index(
                "CH4")]
            Molar_flow_O2used_CH4 = Molar_components_flows_syngas[components.index(
                "CH4")]*2

            # CO combustion
            Molar_flow_CO2_CO = Molar_components_flows_syngas[components.index(
                "CO")]
            Molar_flow_O2used_CO = Molar_components_flows_syngas[components.index("CO")]/2

            Molar_flow_components_exhaust = Molar_components_flows_syngas.copy()

            Molar_flow_components_exhaust[components.index("H2")] = 0
            Molar_flow_components_exhaust[components.index("CH4")] = 0
            Molar_flow_components_exhaust[components.index("CO")] = 0
            Molar_flow_components_exhaust[components.index(
                "H2O")] += Molar_flow_H2O_H2 + Molar_flow_H2O_CH4
            
            Molar_flow_components_exhaust[components.index(
                "CO2")] += Molar_flow_CO2_CO + Molar_flow_CO2_CH4
            
            Molar_flow_components_exhaust[components.index("O2")] +=Molar_components_flows_air[air_composition.index("O2")]-Molar_flow_O2used_H2-Molar_flow_O2used_CH4-Molar_flow_O2used_CO
            Molar_flow_components_exhaust[components.index("O2")]=max(0,Molar_flow_components_exhaust[components.index("O2")])
            
            Molar_flow_components_exhaust[components.index("N2")] +=Molar_components_flows_air[air_composition.index("N2")]
            Molar_flow_components_exhaust[components.index("Ar")] +=Molar_components_flows_air[air_composition.index("Ar")]
            
            Molar_composition_exhaust=list(map(lambda x: x/sum(Molar_flow_components_exhaust),Molar_flow_components_exhaust))
            
            Sr_molar_mass_exhaust=sum(map(lambda x1, x2: x1*molar_mases[x2], Molar_composition_exhaust, components))
            G_exhaust_calculated=Sr_molar_mass_exhaust*sum(Molar_flow_components_exhaust)
            Error_G=abs((G_exhaust_calculated-Gcomb)/Gcomb*100)
            
            Error_composition=abs((Molar_composition_exhaust[0]-combfrac[0])/combfrac[0]*100)
            
            Gair=Gair-Error_composition
            Gair_it.append(Gair)
            Error_Gair = abs((Gair_it[-1]-Gair_it[-2])/Gair_it[-1]*100)
            
            combfrac=Molar_composition_exhaust
            if Gair<0:
                print("Расход воздуха меньше 0")
                print(f"Gsg:{Gsg}, Hsg:{Hsg}, Hcomb:{Hcomb}, Qnr_mass:{Qnr_mass}, Hair2:{Hair2}")
            
            if Error_G > 0.1:
                print("Не сходится баланс камеры сгорания",G_exhaust_calculated,Gcomb, G_exhaust_calculated- Gcomb)
                       
            if i == max_iterations_gas_turbine-1:
                print("Достигнуто максимальное количество итераций расхода воздуха в пиковую ГТУ: ",
                      max_iterations_gas_turbine)
                print("Gair_it ", Gair_it)
            if Error_Gair < Gair_tol and Error_composition < Gair_tol:
                print("Рассчитана пиковая ГТУ, расход воздуха:",Gair_it[-1])
                break

        # Qnr = 23375.32317

        # турб
        Scomb = COMB.p_h(P2, Hcomb)['s']
        Hcombext = COMB.p_s(0.1, Scomb)['h']
        Hcombex = Hcomb - (Hcomb-Hcombext)*KPDturb
        Tcombex = COMB.p_h(0.1, Hcombex)['T']

        Tex = water_streams.loc['SWIN', 'T']+20
        Hex = COMB.p_t(0.1, Tex)['h']

        Ncomp = Gair*(Hair2-Hair1)
        Nturb = Gcomb*(Hcomb - Hcombex)
        Ngtu = Nturb-Ncomp
        Qgtu_ex = Gcomb*(Hcombex-Hex)

        syngas_streams.loc['AIR-COMP', 'T':'G'] = [15, 0.1, Hair1, Gair]
        syngas_streams.loc['COMP-COMB', 'T':'G'] = [Tair2, P2, Hair2, Gair]

        syngas_streams.loc['COMB-TURB', 'T':'G'] = [Tcomb, P2, Hcomb, Gcomb]
        syngas_streams.loc['TURB-COOL',
                           'T':'G'] = [Tcombex, 0.1, Hcombex, Gcomb]
        syngas_streams.loc['COOL-EX', 'T':'G'] = [Tex, 0.1, Hex, Gcomb]
        syngas_streams.loc['COMB-TURB', 'N2':'CO'] = combfrac
        syngas_streams.loc['TURB-COOL', 'N2':'CO'] = combfrac
        syngas_streams.loc['COOL-EX', 'N2':'CO'] = combfrac

        heaters.loc['SGgtu_Qsw', 'Qw'] = Qgtu_ex
        electric.loc['SGgtu_comp', 'Ni'] = Ncomp
        electric.loc['SGgtu_turb', 'Ni'] = Nturb


class accum:

    def __init__(self, water_streams, accumulation):
        self._V = 1
        self._D = 1
        self._F = 1
        self._H = 1
        self._D = 1
        self._kolichestvo = 1
        self._V = n.pi*self._D**3/4
        self._F = 1.5*n.pi*self._D**2
        self._khi = 1
        self._lambda_min_vata = 0.045
        self.delta_min_vata = 0.01
        self.water_streams = water_streams
        self.accumulation = accumulation
        self._T_nar_vozd = self.water_streams.at['AIR', 'T']
        self._f = 0
        self._T_accum = 0
        self._h_accum = 0
        self._P_accum = 0
        self._Mass = 0
        self._Q = 0

    def set_construct(self, constr):

        self._D = constr['Diametr']
        self._kolichestvo = constr['kolichestvo']
        self._H = constr['Visota']
        self._lambda_min_vata = constr['lambda_min_vata']
        self._delta_min_vata = constr['delta_min_vata']
        self._V = n.pi*self._D**2/4 * self._H
        self._F = n.pi*self._D*self._H + n.pi*self._D**2/2

        pass

    def zaryad(self, t, accumulation, gas_streams, syngas_streams, water_streams, water_streams0, heaters, electric):
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGfrac = (0, 0, 0.168738364456343, 0, 0.0503573198627571,
                  0.0975081144748292, 0.681387369772999, 0.00200883143307171)
        SG = prop.Materials_prop(SGsost, SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)
        V = self._V
        Tsg = 100
        Psg = 2
        rosg = SG.p_t(Psg, Tsg)['rho']
        Gsg = V*self._kolichestvo*rosg/(t*3600)
        def sved_Gst(Gst):
            water_streams.loc['DROSVD-ST', 'G'] = Gst
            PKM = PKM_all.calc(True, gas_streams, syngas_streams,
                               water_streams, water_streams0, heaters)
            Gsgnew = syngas_streams.loc["Separ-SGaccum", "G"]
            return ((Gsg-Gsgnew)/Gsg)*100
        Gst = float(
            root(sved_Gst, water_streams.loc['DROSVD-ST', 'G'], method="hybr", tol=10**-2).x)
        if Gst > water_streams.loc['PEVD-DROSVD', 'G']:
            print('!Расчетный расход пара в ПКМ больше расхода пара в КУ:',Gst, water_streams.loc['PEVD-DROSVD', 'G'])
            Gst = water_streams.loc['PEVD-DROSVD', 'G']
        water_streams.loc['DROSVD-ST', 'G'] = Gst
        water_streams.loc['DROSVD-ST','T':'H'] = water_streams.loc['PEVD-DROSVD', 'T':'H']
        PKM = PKM_all.calc(True, gas_streams, syngas_streams,water_streams, water_streams0, heaters)
        steamVD_to_turbine = water_streams.loc['PEVD-DROSVD','G'] - water_streams.loc['DROSVD-ST', 'G']

        Qteplofic = water_streams.loc['SWIN', 'G'] * (water_streams.loc['SWOUT', 'H']-water_streams.loc['SWIN', 'H'])
        Qref = heaters.loc["Ref_all", "Qw"]
        if Qref < Qteplofic:
            Gw_pkm = Qref / (water_streams.loc['SWOUT', 'H'] - water_streams.loc['SWIN', 'H'])
            water_streams.loc['SWIN-TURB','G'] = water_streams.loc['SWIN', 'G']-Gw_pkm
            print('Тепла от ПКМ не хватает на теплофикацию. Расход на СП:',
                  round(water_streams.loc['SWIN-TURB', 'G'], 1))
            Teplo = 1
        else:
            print('Тепла от ПКМ хватает на теплофикацию')
            Teplo = 0

        return {'steamVD_to_turbine': steamVD_to_turbine, 'Teplo': Teplo}

    def razryad(self, t, accumulation, gas_streams, syngas_streams, water_streams, water_streams0, heaters, electric):

        PKM = PKM_all.calc(False, gas_streams, syngas_streams,
                           water_streams, water_streams0, heaters)
        steamVD_to_turbine = water_streams.loc['PEVD-DROSVD', 'G']
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        SGmass = syngas_streams.loc['Separacc-Sepout', 'G']  # масса
        SGfrac = list(syngas_streams.loc['Separacc-Sepout', 'N2':'CO'])

        SG = prop.Materials_prop(SGsost, SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)

        Tsg = syngas_streams.loc['Separacc-Sepout', 'T']
        Psg = syngas_streams.loc['Separacc-Sepout', 'P']
        rosg = SG.p_t(Psg, Tsg)['rho']
        Hsg = SG.p_t(Psg, Tsg)['h']
        Gsg = SGmass/(t*3600)
        if Gsg <= 9.73*1/9:
            Gsg_GTUmain = Gsg
            Gsg_GTUsg = 0
        else:
            Gsg_GTUmain = 9.73*1/9
            Gsg_GTUsg = Gsg - Gsg_GTUmain
        syngas_streams.loc["Sepout-COMB", "T":"G"] = [Tsg, Psg, Hsg, Gsg_GTUsg]
        syngas_streams.loc["Sepout-GTU", "T":"H"] = [Tsg, Psg, Hsg]
        syngas_streams.loc["Sepout-GTU", "G"] = Gsg_GTUmain
        syngas_streams.loc["Sepout-GTU", "N2":"CO"] = SGfrac
        syngas_streams.loc["Sepout-COMB", "N2":"CO"] = SGfrac

        Qteplofic = water_streams.loc['SWIN', 'G'] * \
            (water_streams.loc['SWOUT', 'H']-water_streams.loc['SWIN', 'H'])
        syngas_GTU.calc(water_streams,syngas_streams, 1.2, 0.88, 0.9, heaters, electric,accumulation)
        Qgvto = heaters.loc['SGgtu_Qsw', 'Qw']
        print("Qgvto:",Qgvto)
        if Qgvto < Qteplofic:
            Gw_pkm = Qgvto / \
                (water_streams.loc['SWOUT', 'H'] -
                 water_streams.loc['SWIN', 'H'])
            water_streams.loc['SWIN-TURB',
                              'G'] = water_streams.loc['SWIN', 'G']-Gw_pkm
            Teplo = 1
            print('Тепла в ГВТО не хватает на теплофикацию. Расход на СП:',
                  water_streams.loc['SWIN-TURB', 'G'])

        else:
            print('Тепла в ГВТО хватает на теплофикацию')
            Teplo = 0

        return {'steamVD_to_turbine': steamVD_to_turbine, 'Teplo': Teplo}

    def jdat(self, time_jdat, accumulation, gas_streams, syngas_streams, water_streams, heaters, electric):
        SGsost = "Nitrogen*O2*CO2*Ar*H2O*Methane*H2*CO"
        ##########НАДО СЧИТАТЬ С ТАБЛИЦЫ##########
        SGfrac = list(syngas_streams.loc['Separ-SGaccum', 'N2':])
        if isinstance(SGfrac[-2], float):
            SGfrac = (0, 0, 0.168738364456343, 0, 0.0503573198627571,
                      0.0975081144748292, 0.681387369772999, 0.00200883143307171)

        SG = prop.Materials_prop(SGsost, SGfrac,
                                 prop.REFPROP_h_s,
                                 prop.REFPROP_p_t,
                                 prop.REFPROP_p_h,
                                 prop.REFPROP_p_s,
                                 prop.REFPROP_p_q,
                                 prop.REFPROP_t_q,
                                 prop.REFPROP_p_rho,
                                 prop.REFPROP_s_q,
                                 RP=RP)
        T_nar_vozd = water_streams.loc['AIR', 'T']
        T_accum = syngas_streams.loc['Separ-SGaccum', 'T']
        Psg = syngas_streams.loc['Separ-SGaccum', 'P']
        ######НАДО СЧИТАТЬ С ТАБЛИЦЫ######
        if isinstance(T_accum, float):
            T_accum = 100  # после сепаратора HTS
        if isinstance(Psg, float):
            Psg = 2

        rosg = SG.p_t(Psg, T_accum)['rho']

        n_step = 100
        time_jdat_step = time_jdat/n_step
        H1 = SG.p_t(Psg, T_accum)['h']
        H2 = SG.p_t(Psg, T_accum)['h']
        

            # # заменить на логарифм
            # poteri = self._khi*(self._lambda_min_vata/self.delta_min_vata)*self._F * \
            #     self._kolichestvo * (T_accum-T_nar_vozd) * \
            #     (time_jdat_step*3600)/1000
            # H2 = H2 - poteri/(self._V*self._kolichestvo*rosg)
            # T_accum = SG.p_h(Psg, H2)['T']
            
            
        for i in range(n_step):
            poteri = time_jdat_step*3600/1000*3.14*self._H*(T_accum-self._T_nar_vozd)/(1/2*self._lambda_min_vata*n.log((self._D+2*self._delta_min_vata)/self._D)+1/(100000*self._D)+1/(20*(self._D+2*self._delta_min_vata)))+ time_jdat_step*3600/1000*2*3.14*self._D**2/4*(self._T_accum-self._T_nar_vozd)/(1/20+1/100000+self._delta_min_vata/self._lambda_min_vata) #kJ
            H2 = H2 - poteri/(self._V*self._kolichestvo*rosg)
            T_accum = SG.p_h(Psg, H2)['T']

        self.accumulation.at["PKM", "Qw"] = (
            self._V*self._kolichestvo*rosg)*(H2)  # хз надо ли
        self.accumulation.at["PKM", "T"] = T_accum
        syngas_streams.loc["SGaccum-Separacc", "T"] = T_accum
        syngas_streams.loc["SGaccum-Separacc", "P"] = Psg
        syngas_streams.loc["SGaccum-Separacc", "H"] = H2

        syngas_streams.loc["SGaccum-Separacc",
                           "G"] = (self._V*self._kolichestvo*rosg)  # Масса
        syngas_streams.loc['SGaccum-Separacc', 'N2':'CO'] = SGfrac
        # from PKM import separator
        sep = separator.calc('SGaccum-Separacc',
                             'Separacc-Sepout', syngas_streams, heaters)

        syngas_streams.loc["Sepout-COMB",
                           "T":"CO"] = syngas_streams.loc["Separacc-Sepout", "T":"CO"]
        return {'T_accum': T_accum, 'poteri': poteri}
