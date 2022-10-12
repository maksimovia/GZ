import mat_properties as prop


# SR = steam_transformer(stream11 = 'stream11', stream12 = 'streams12', stream21 = 'stream21', stream22 = 'stream22', water = 'water', water_streams = 'water_streams', heaters = 'heaters', Pdr1 = 2, Pdr2 = 0.8, P2 = 2, dT = 15, dTmin = 5, Tdec = 10)

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

        T16 = T15 - self.Tdec
        P16 = P15
        H16 = self.water.p_t(P16, T16)['h']

        P22 = P2
        H22 = H23 - (G1*(H13-H14)/G2)
        T22 = self.water.p_h(P22, H22)['T']

        P21 = P2
        H21 = H22 - (G1*(H15-H16)/G2)
        T21 = self.water.p_h(P21, H21)['T']
        
        
        
        
        return {'T12': ,}
        
        
        














class reformer:
    def __init__(self, stream11, gas, water, Methane, waterMethane, gas_streams, water_streams, heaters, P11, P12, P13, P2, Tref):
        
        
        
        