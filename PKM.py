import mat_properties as prop

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
        
        return {'H11': H11,'H12': H12,'H13': H13,'H14': H14,'H15': H15,'H16': H16,
               'H21': H21,'H22': H22,'H23': H23,'H24': H24, 'G1':G1, 'G2':G2,'Q':Qtrans, 'P2':self.P2, 'P16':P16,
               'T11': T11,'T12': T12,'T13': T13,'T14': T14,'T15': T15,'T16': T16,
               'T21': T21,'T22': T22,'T23': T23,'T24': T24}

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
        
        Qreac = (52166.505091208484/44.6397313913235)*Gref
        Qref = Qdt+Qreac
#         H1gas = self.gas_KU.p_t(0.1, 1968.58395330148)['h']
        H1gas = 3191.2518095937544
        H2gas = self.gas_KU.p_t(0.1, self.T2gas)['h']
        
        Ggas = Qref/(H1gas-H2gas)
        Gair = Ggas*48.5966429877363/51.2672005145991
        Gch4 = Ggas*2.67055752686283/51.2672005145991
        Hsg = H2r + (Qreac/Gref)
                
        SGsost = ['N2','O2','CO2','Ar','H2O','CH4','H2','CO']
        SGfrac = [0,1.96010880218854E-22,0.06356655972348,0,0.536615036955725,0.0500421242206529,0.325898769048586,0.0238775100515557]
        SGfrac = dict(zip(SGsost,SGfrac))

        Gassost = ['N2','O2','CO2','H2O','Ar']
        Gasfrac = [0.710320591016015,0.00996710270335893,0.090538556815177,0.180531273012258,0.00864247645319178]
        Gasfrac = dict(zip(Gassost,Gasfrac))
        
        return {'Q':Qref,'Ggas':Ggas,'Gair':Gair,'Gch4':Gch4, 'Tref':self.Tref,
                'Gref':Gref, 'Pref':self.Pref,'Hsg':Hsg,'SGfrac':SGfrac,'H2gas':H2gas,'Gasfrac':Gasfrac}