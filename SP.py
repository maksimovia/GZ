import mat_properties as prop
import numpy as n

class od:
    def __init__(self,stream11,stream12,stream21,stream22,KPD,calctolerance,water,calcmethod,water_streams0,water_streams):
        self.KPD=KPD
        self.water_streams0=water_streams0
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.stream11=stream11
        self.stream21=stream21
        self.stream12=stream12
        self.stream22=stream22
        self.water=water
        delta_0=10 #задаемся недогревом на холодном конце, чтобы оценить расход        
        self.Cp=4.187
        Qod_0=6.81*self.Cp/3.6 #нагрузка охладителя дренажа (расчетная точка)
        tdr_in_0 = water_streams0.at[self.stream11,'T'] #температура горячего потока на входе 
        G_dr_0=water_streams0.at[self.stream11,'G'] #расход греющего потока
        tw1_0= water_streams0.at[self.stream21,'T']#температура холодного потока на входе
        tod_out_0=tdr_in_0-delta_0 #температура холодного потока на выходе
        tdr2_0=tdr_in_0-Qod_0*1000/G_dr_0/self.Cp
        deltat=((tdr_in_0-tod_out_0)-(tdr2_0-tw1_0))/n.log((tdr_in_0-tod_out_0)/(tdr2_0-tw1_0))
        self.F=n.sqrt((tdr_in_0-tdr2_0)*(tod_out_0-tw1_0))/deltat
    def calc(self):
        Gdr = self.water_streams.at[self.stream11,'G']
        God = self.water_streams.at[self.stream21,'G']
        tdr_in = self.water_streams.at[self.stream11,'T']
        tw_in = self.water_streams.at[self.stream21,'T']
        Delta=tdr_in-tw_in
        epsilon=min(0.999, 1/(0.35*min(Gdr,God)/max(Gdr,God)+0.65+1/self.F*n.sqrt(min(Gdr,God)/max(Gdr,God))))
        Q=epsilon*min(Gdr,God)*self.Cp*Delta/1000
        tod_out=tw_in+Q*1000/God/self.Cp
        tdr_out=tdr_in-Q*1000/Gdr/self.Cp
        print(God)
        return [tdr_out,tod_out]
    


class sp2:
    def __init__(self,stream11,stream12,stream21,stream22,KPD,calctolerance,water,calcmethod,water_streams0,water_streams):
        self.KPD=KPD
        self.water_streams0=water_streams0
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.stream11=stream11
        self.stream21=stream21
        self.stream12=stream12
        self.stream22=stream22
        self.water=water
        self.Cp=4.187
        self.dP=5/100
        tsp2_0=water_streams0.at[self.stream12,'T'] #температура насыщения в номинале 
        tsp2_in_0 = water_streams0.at[self.stream21,'T'] #температура холодного поткоа на входе в номинале 
        tw2_0 = water_streams0.at[self.stream22,'T'] #температура воды на выходе
        Gsv_0 =water_streams0.at[self.stream21,'G'] #расход нагреваемого потока
        self.KF=Gsv_0*self.Cp*n.log((tsp2_0-tsp2_in_0)/(tsp2_0-tw2_0))

    def calc(self):
        Gsv=self.water_streams.at[self.stream21,'G'] #расход нагреваемого потока
        tw_out= self.water_streams.at[self.stream22,'T'] #температура воды на выходе
        tw_in = self.water_streams.at[self.stream21,'T'] #температура холодного поткоа на входе
        h_otb_sp= self.water_streams0.at[self.stream11,'H'] #температура воды на выходе
        m=n.exp(-self.KF/Gsv/self.Cp)
        
        tsp=(tw_out-tw_in*m)/(1-m)
        Nas_sp=self.water.t_q(tsp,0)
        hsp_nas=Nas_sp['h']
        p_sp=Nas_sp['P']
        p_otb=p_sp/(1-self.dP)
        Qsp=Gsv*self.Cp*(tw_out-tw_in)/self.KPD/1000
        Gotb=Qsp*1000/(h_otb_sp-hsp_nas)    
        print(h_otb_sp)
        return [p_otb,Gotb]
    



#делаю класс СП1 (верхний отбор по пару, второй по воде)

class sp1:
    def __init__(self,stream11,stream12,stream21,stream22,stream111, KPD,calctolerance,water,calcmethod,water_streams0,water_streams):
        self.KPD=KPD
        self.water_streams0=water_streams0
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.stream11=stream11
        self.stream21=stream21
        self.stream12=stream12
        self.stream22=stream22
        self.stream111=stream111
        self.water=water
        self.Cp=4.187
        self.dP=5/100
        tsp1_0=water_streams0.at[self.stream12,'T'] #температура насыщения в номинале 
        tsp1_in_0 = water_streams0.at[self.stream21,'T'] #температура холодного потока на входе в номинале 
        tw1_0 = water_streams0.at[self.stream22,'T'] #температура воды на выходе
        Gsv_0 =water_streams0.at[self.stream21,'G'] #расход нагреваемого потока
        self.KF=Gsv_0*self.Cp*n.log((tsp1_0-tsp1_in_0)/(tsp1_0-tw1_0))

    def calc(self):
        Gsv= self.water_streams.at[self.stream21,'G'] #расход нагреваемого потока
        tw_in = self.water_streams.at[self.stream21,'T'] #температура холодного поткоа на входе
        h_otb_sp= self.water_streams0.at[self.stream11,'H'] #температура воды на выходе
        p_otb_sp= self.water_streams0.at[self.stream11,'P'] #температура воды на выходе
        h_nas_sp2= self.water_streams0.at[self.stream111,'H'] #температура воды на выходе
        G_nas_sp2= self.water_streams0.at[self.stream111,'G'] #температура воды на выходе
        m=n.exp(-self.KF/Gsv/self.Cp)
        p_sp=p_otb_sp*(1-self.dP)
        Nas_sp=self.water.p_q(p_sp,0)
        t_sp=Nas_sp['T']
        hsp_nas=Nas_sp['h']
        tw_out=t_sp-(t_sp-tw_in)*m     
        Qsp=Gsv*self.Cp*(tw_out-tw_in)/self.KPD/1000
        Gotb=(Qsp*1000-G_nas_sp2*(h_nas_sp2-hsp_nas))/(h_otb_sp-hsp_nas)    
        print(h_otb_sp)
        return [tw_out,Gotb]





