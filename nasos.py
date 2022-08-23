class nasos():
    def __init__(self,stream11,stream12,water, KPDnasos,  water_streams, water_streams0):
        self.KPDnasos = KPDnasos
        self.stream11=stream11
        self.stream12=stream12
        self.water=water
        self.water_streams=water_streams
        self.water_streams0=water_streams0
        self.G0=water_streams0.at[self.stream11,'G']
        self.P01  = water_streams0.at[self.stream11,'P']
        self.P02  = water_streams0.at[self.stream12,'P']
        
    def KPD(self,G):
        g_otn=G/self.G0
        if g_otn<0.65:
            g_otn=0.65
        KPD_otn =(603.31*g_otn**4 - 1902.8*g_otn**3 + 2241.5*g_otn**2 - 1168*g_otn + 227.13)
        if KPD_otn > 1:
            KPD = self.KPDnasos
        else:
            KPD = (603.31*g_otn**4 - 1902.8*g_otn**3 + 2241.5*g_otn**2 - 1168*g_otn + 227.13) * self.KPDnasos
        return KPD
    
    
    def KPDm (self):
        P1  = self.water_streams.at[self.stream11,'P']
        P2  = self.water_streams.at[self.stream12,'P']
        G=self.water_streams.at[self.stream11,'G']
        H_otn = (P2-P1)/(self.P02-self.P01)
        Q_otn = G/self.G0
        H0 = -0.1957 * Q_otn - 0.1254 * Q_otn + 1.3138
        n = (H_otn/H0 -0.0233)/0.9402
        KN = 0.9402 * n + 0.0233
        if KN > 1:
            KPDgm = 0.9736 * KN**(-0.3358)
        if KN>0:
            KPDgm = 0.9736 * KN**(0.3358)
        else:
            KPDgm=0.1
        return KPDgm
    
    
    def calc(self):
        
        P1  = self.water_streams.at[self.stream11,'P']
        P2  = self.water_streams.at[self.stream12,'P']
        h1 = self.water_streams.at[self.stream11,'H']
        T1 = self.water.p_q(P1,0)['T']
        G1  = self.water_streams.at[self.stream11,'G']
        s1 = self.water.p_q(P1,0)['s']
        h2teor = self.water.p_s(P2,s1)['h']
        h2real = h1+(h2teor-h1)/self.KPD(G1)
        T2=self.water.p_h(P2,h2real)['T']
        Rabota = G1*(h2real - h1)/1000
        Ngm = Rabota/self.KPDm()
        return {'T2': T2, 'P2': P2, 'h2real': h2real, 'G1': G1, 'Ni': Rabota, 'Ngm': Ngm, 'KPD': self.KPD(G1), 'KPDm': self.KPDm()}
