class nasos():
    def __init__(self,stream11,stream12,water, KPDnasos,  water_streams, water_streams0):
        self.KPDnasos = KPDnasos
        self.stream11=stream11
        self.stream12=stream12
        self.water=water
        self.water_streams=water_streams
        self.water_streams0=water_streams0
        self.G0=water_streams0.at[self.stream11,'G']
    
        
    def KPD(self,G):
        g_otn=G/self.G0
        KPD_otn =(603.31*g_otn**4 - 1902.8*g_otn**3 + 2241.5*g_otn**2 - 1168*g_otn + 227.13)
        if KPD_otn > 1:
            KPD = self.KPDnasos
        else:
            KPD = (603.31*g_otn**4 - 1902.8*g_otn**3 + 2241.5*g_otn**2 - 1168*g_otn + 227.13) * self.KPDnasos
        return KPD
    
    def calc(self):
        
        P1  = self.water_streams.at[self.stream11,'P']
        P2  = self.water_streams.at[self.stream12,'P']
        h1 = self.water_streams.at[self.stream11,'H']
        T1 = self.water.p_q(P1,0)['T']
        G1  = self.water_streams.at[self.stream11,'G']
        # print (G1,'G1')
        s1 = self.water.p_q(P1,0)['s']
        h2teor = self.water.p_s(P2,s1)['h']
        h2real = h1+(h2teor-h1)/self.KPD(G1)
        # print (h2real,'h2')
        T2=self.water.p_h(P2,h2real)['T']
        Rabota = G1*(h2real - h1)
        return [T2, P2, h2real, G1, Rabota, self.KPD(G1)]