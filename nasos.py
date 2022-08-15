class nasos:
    def __init__(self,stream11,stream12,water,KPDnasos,water_streams):
        self.KPDnasos = KPDnasos
        self.stream11=stream11
        self.stream12=stream12
        self.water=water
        self.water_streams=water_streams
    def calc(self):
        P0  = self.water_streams.at[self.stream11,'P']
        P1  = self.water_streams.at[self.stream12,'P']
        H0 = self.water_streams.at[self.stream11,'H']
        T0 = self.water.p_h(P0,H0)['T']
        G  = self.water_streams.at[self.stream11,'G']
        s1teor = self.water.p_t(P0,T0)['s']
        h1teornasos = self.water.p_s(P1,s1teor)['h']
        h1realnasos = self.water.p_t(P0,T0)['h']+(h1teornasos-self.water.p_t(P0,T0)['h'])/self.KPDnasos
        T1nasosCO2=self.water.p_h(P1,h1realnasos)['T']
        Rabota = G*(h1realnasos - self.water.p_t(P0,T0)['h'])
        return [T1nasosCO2,P1,h1realnasos,G,Rabota]