import mat_properties as prop
import pandas as pd


RP = prop.init_REFPROP(r'C:\Program Files (x86)\REFPROP')



water=prop.Materials_prop('water',[1.0,0,0,0,0], prop.REFPROP_h_s,
                                              prop.REFPROP_p_t,
                                              prop.REFPROP_p_h,
                                              prop.REFPROP_p_s,
                                              prop.REFPROP_p_q,
                                              prop.REFPROP_t_q,
                                              RP=RP)



class nasos:
    def __init__(self,stream11,stream12,water,KPDnasos,water_streams):
        self.KPDnasos = KPDnasos
        self.P0  = water_streams.at[stream11,'P']
        self.P1  = water_streams.at[stream12,'P']
        self.T0 = water_streams.at[stream11,'T']
        self.G  = water_streams.at[stream11,'G']
        self.water=water
        
    def calc(self):
        self.s1teor = self.water.p_t(self.P0,self.T0)['s']
        self.h1teornasos = self.water.p_s(self.P1,self.s1teor)['h']
        self.h1realnasos = self.water.p_t(self.P0,self.T0)['h']+(self.h1teornasos-self.water.p_t(self.P0,self.T0)['h'])/self.KPDnasos
        self.T1nasosCO2=self.water.p_h(self.P1,self.h1realnasos)['T']
        T1 = self.T1nasosCO2
        Rabota = self.G*(self.h1realnasos - self.water.p_t(self.P0,self.T0)['h'])
        return [Rabota,T1,self.P0,self.P1,self.T0,self.G]