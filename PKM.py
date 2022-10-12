import mat_properties as prop

class steam_transformer:
    def __init__(self,**kwargs):
        self.stream11 = kwargs['stream11']
        self.stream12 = kwargs['stream12']
        self.stream21 = kwargs['stream21']
        self.stream22 = kwargs['stream22']
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
        G1 = self.water_streams.at[self.stream11, 'G']

        Q12 = 0
        P12 = P11
        H12 = self.water.p_q(P12, Q12)['h']
        T12 = self.water.p_q(P12, Q12)['T']

        P13 = self.Pdr1
        H13 = H12
        T13 = water.p_h(P13, H13)['T']
        Q13 = water.p_h(P13, H13)['Q']
    
        P22 = self.P2
        T22 = T11-self.dT
        H22 = water.p_t(P22, T22)['h']

        T21 = T13-self.dTmin
        P21 = P22
        H21 = water.p_t(P21, T21)['h']

        G2 = G1 * ((H11-H12)/(H22-H21))
        
        Q14 = 0
        P14 = P13
        H14 = self.water.p_q(P14, Q14)['h']
        T14 = self.water.p_q(P14, Q14)['T']
        
        P15 = self.Pdr2
        H15 = H14
        T15 = water.p_h(P15, H15)['T']
        Q15 = water.p_h(P15, H15)['Q']
        
        T16 = T15 - self.Tdec
        P16 = P15
        H16 = water.p_t(P16, T16)['h']
        
        
        
        
        return {'T12': ,}
        
        
        














class reformer:
    def __init__(self, stream11, gas, water, Methane, waterMethane, gas_streams, water_streams, heaters, P11, P12, P13, P2, Tref):
        
        
        
        