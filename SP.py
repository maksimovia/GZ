import mat_properties as prop
import numpy as n

def Tset(Tnv):
    if Tnv <=-23.459064:
        Tpr=120
    if Tnv>-23.459064 and Tnv<=6.190058:
        Tpr=-1.01183432*Tnv+96.26331361
    if Tnv>6.190058:
        Tpr=90
    if Tnv <=-27.1:
        Tobr=70
    if Tnv>-27.1 and Tnv<=6.190058:
        Tobr=-1*Tnv+42.9
    if Tnv>6.190058:
        Tobr=0.1067856*Tnv+36.04893297
    return Tpr, Tobr

class od:
    def __init__(self,stream11,stream12,stream21,stream22,KPD,water,water_streams0,water_streams):
        self.KPD=KPD
        self.water_streams0=water_streams0
        self.water_streams=water_streams
        self.stream11=stream11
        self.stream21=stream21
        self.stream12=stream12
        self.stream22=stream22
        self.water=water
        delta_0=10 #задаемся недогревом на холодном конце, чтобы оценить расход        
        self.Cp=4.187
        # Qod_0=6.81*self.Cp/3.6 #нагрузка охладителя дренажа (расчетная точка)
        Qod_0=6.81*self.Cp*self.KPD/3.6 #нагрузка охладителя дренажа (расчетная точка)
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
        Hdr_in = self.water_streams.at[self.stream11,'H']
        tw_in = self.water_streams.at[self.stream21,'T']
        pdr_in = self.water_streams.at[self.stream11,'P']
        pdr_out=pdr_in
        Delta=tdr_in-tw_in
        epsilon=min(0.999, 1/(0.35*min(Gdr,God)/max(Gdr,God)+0.65+1/self.F*n.sqrt(min(Gdr,God)/max(Gdr,God))))
        Q=min(epsilon,1)*min(Gdr,God)*self.Cp*Delta/1000
        tod_out=tw_in+Q*1000*self.KPD/God/self.Cp
        Hdr_out=Hdr_in-Q*1000/Gdr
        tdr_out=self.water.p_h(pdr_out,Hdr_out)['T']

        self.water_streams.at[self.stream12,'G']= Gdr
        self.water_streams.at[self.stream22,'G']= God    
 
        self.water_streams.at[self.stream12,'H']= Hdr_out
        self.water_streams.at[self.stream22,'H']= self.Cp*tod_out   
        
        self.water_streams.at[self.stream12,'T']= tdr_out
        self.water_streams.at[self.stream22,'T']= tod_out 
        
        self.water_streams.at[self.stream12,'P']= self.water_streams.at[self.stream11,'P']
        self.water_streams.at[self.stream22,'P']= self.water_streams.at[self.stream21,'P']   
        
        res = dict()
        Qw=Q*self.KPD
        res['tdr_out']=tdr_out
        res['tod_out']=tod_out
        res['Qw']=Qw
        res['epsilon']=epsilon
      

        return res
    


class sp2:
    def __init__(self,stream11,stream12,stream21,stream22,KPD,water,water_streams0,water_streams):
        self.KPD=KPD
        self.water_streams0=water_streams0
        self.water_streams=water_streams
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
        self.Gsv_0 =water_streams0.at[self.stream21,'G'] #расход нагреваемого потока
        self.KF=self.Gsv_0*self.Cp*n.log((tsp2_0-tsp2_in_0)/(tsp2_0-tw2_0))

    
    def KF_otn(self, Gsv):
        Gsv_0=self.Gsv_0
        G_otn=Gsv/self.Gsv_0
        KF_ot=1.2285*(G_otn)**4-4.8973*(G_otn)**3+6.3132*(G_otn)**2-2.6031*(G_otn)**1+0.96
        # KF_ot=1
        return KF_ot
        
    def calc(self):
        Gsv=self.water_streams.at[self.stream21,'G'] #расход нагреваемого потока
        tw_out= self.water_streams.at[self.stream22,'T'] #температура воды на выходе
        tw_in = self.water_streams.at[self.stream21,'T'] #температура холодного поткоа на входе
        h_otb_sp= self.water_streams.at[self.stream11,'H'] 
        m=n.exp(-self.KF*self.KF_otn(Gsv)/Gsv/self.Cp)
        tsp=(tw_out-tw_in*m)/(1-m)
        Nas_sp=self.water.t_q(tsp,0)
        hsp_nas=Nas_sp['h']
        p_sp=Nas_sp['P']
        p_otb=p_sp/(1-self.dP)
        Qsp=Gsv*self.Cp*(tw_out-tw_in)/1000
        Gotb=Qsp*1000/(h_otb_sp-hsp_nas)/self.KPD

        Qw=Qsp
        
        self.water_streams.at[self.stream12,'G']= Gotb
        self.water_streams.at[self.stream11,'G']= Gotb
        self.water_streams.at[self.stream22,'G']= Gsv   
        
        self.water_streams.at[self.stream12,'H']= hsp_nas
        self.water_streams.at[self.stream22,'H']= self.Cp*tw_out   
        
        self.water_streams.at[self.stream12,'T']= tsp  
        
        self.water_streams.at[self.stream12,'P']= p_sp 
        self.water_streams.at[self.stream11,'P']= p_otb
        self.water_streams.at[self.stream22,'P']= self.water_streams.at[self.stream21,'P']   
        
        res = dict()
        res['p_otb']=p_otb
        res['Gotb']=Gotb
        res['Qw']=Qw

        return res
    



#делаю класс СП1 (верхний отбор по пару, второй по воде)

class sp1:
    def __init__(self,stream11,stream12,stream21,stream22,stream111, KPD,water,water_streams0,water_streams):
        self.KPD=KPD
        self.water_streams0=water_streams0
        self.water_streams=water_streams
        self.stream11=stream11
        self.stream21=stream21
        self.stream12=stream12
        self.stream22=stream22
        self.stream111=stream111
        self.Gsv_0 =water_streams0.at[self.stream21,'G'] #расход нагреваемого потока
        self.water=water
        self.Cp=4.187
        self.dP=5/100
        tsp1_0=water_streams0.at[self.stream12,'T'] #температура насыщения в номинале 
        tsp1_in_0 = water_streams0.at[self.stream21,'T'] #температура холодного потока на входе в номинале 
        tw1_0 = water_streams0.at[self.stream22,'T'] #температура воды на выходе
        Gsv_0 =water_streams0.at[self.stream21,'G'] #расход нагреваемого потока
        self.KF=Gsv_0*self.Cp*n.log((tsp1_0-tsp1_in_0)/(tsp1_0-tw1_0))

    
    def KF_otn(self, Gsv):
        Gsv_0=self.Gsv_0
        G_otn=Gsv/self.Gsv_0
        KF_ot=1.2285*(G_otn)**4-4.8973*(G_otn)**3+6.3132*(G_otn)**2-2.6031*(G_otn)**1+0.96
        # KF_ot=1
        return KF_ot
    
    
    def calc(self):
        Gsv= self.water_streams.at[self.stream21,'G'] #расход нагреваемого потока
        tw_in = self.water_streams.at[self.stream21,'T'] #температура холодного поткоа на входе
        h_otb_sp= self.water_streams.at[self.stream11,'H'] #температура воды на выходе
        p_otb_sp= self.water_streams.at[self.stream11,'P'] #температура воды на выходе
        h_nas_sp2= self.water_streams.at[self.stream111,'H'] #температура воды на выходе
        G_nas_sp2= self.water_streams.at[self.stream111,'G'] #температура воды на выходе
        m=n.exp(-self.KF*self.KF_otn(Gsv)/Gsv/self.Cp)
        p_sp=p_otb_sp*(1-self.dP)
        Nas_sp=self.water.p_q(p_sp,0)
        t_sp=Nas_sp['T']
        hsp_nas=Nas_sp['h']
        if t_sp<tw_in:
            Qw=0
            tw_out=tw_in+G_nas_sp2*(h_nas_sp2-hsp_nas)*self.KPD/Gsv/self.Cp
            Gotb=0
        else:
            tw_out=t_sp-(t_sp-tw_in)*m     
            Qw=Gsv*self.Cp*(tw_out-tw_in)/1000
            Gotb=(Qw*1000-G_nas_sp2*(h_nas_sp2-hsp_nas)*self.KPD)/(h_otb_sp-hsp_nas) / self.KPD
        
        self.water_streams.at[self.stream12,'G']= (G_nas_sp2+Gotb)
        self.water_streams.at[self.stream11,'G']= Gotb
        self.water_streams.at[self.stream22,'G']= Gsv   
        
        self.water_streams.loc[self.stream22, 'T'] = tw_out
        self.water_streams.loc[self.stream12, 'T'] = t_sp
        
        self.water_streams.at[self.stream12,'H']= hsp_nas
        self.water_streams.at[self.stream22,'H']= self.Cp*tw_out   
        
        self.water_streams.at[self.stream12,'P']= p_sp
        self.water_streams.at[self.stream22,'P']= self.water_streams.at[self.stream21,'P']   

        
        res = dict()
        res['tw_out']=tw_out
        res['Gotb']=Gotb
        res['Qw']=Qw

        return res


class teplofik_systema:
    def __init__(
        self, KPD_SP, water, water_streams0, water_streams
    ):

        self.OD = od(
            "SP1-OD",
            "OD-GPK",
            "SWIN-OD",
            "OD-SP1",
            KPD_SP,
            water,
            water_streams0,
            water_streams,
        )
        self.SP1 = sp1(
            "OTB1-SP1",
            "SP1-OD",
            "WIN-SP1",
            "SP1-SP2",
            "SP2-SP1",
            KPD_SP,
            water,
            water_streams0,
            water_streams,
        )
        self.SP2 = sp2(
            "OTB2-SP2",
            "SP2-SP1",
            "SP1-SP2",
            "SP2-WOUT",
            KPD_SP,
            water,
           water_streams0,
            water_streams,
        )

        self.water_streams = water_streams
        self.water_streams0 = water_streams0
        self.KPD_SP = KPD_SP

    def calculate(self, maxiterations, calctolerance):

        for i in range(maxiterations):

            # разделение потоков на ОД и байпас ОД
            Gsv_in = self.water_streams.at["SWIN-TURB", "G"]
            Gsv_in0 = self.water_streams0.at["SWIN-TURB", "G"]
            G_od0 = self.water_streams0.at["SWIN-OD", "G"]
            G_od = G_od0 * Gsv_in / Gsv_in0
            self.water_streams.at["SWIN-OD", "G"] = G_od
            t_ow = self.water_streams.at["SWIN-TURB", "T"]
            self.water_streams.at["SWIN-OD", "T"] = t_ow
            self.water_streams.at["SWIN-OD", "P"] = self.water_streams.at[
                "SWIN-TURB", "P"
            ]

            # Расчет ОД
            OD_res = self.OD.calc()
            Qw_od=OD_res['Qw']

            # Смешение
            t_wod = self.water_streams.at["OD-SP1", "T"]
            t_sp1in = (t_ow * (Gsv_in - G_od) + t_wod * G_od) / Gsv_in
            self.water_streams.at["WIN-SP1", "T"] = t_sp1in
            self.water_streams.at["WIN-SP1", "G"] = Gsv_in
            self.water_streams.at["WIN-SP1", "P"] = self.water_streams.at[
                "SWIN-TURB", "P"
            ]

            # Расчет СП
            SP1_res = self.SP1.calc()
            Qw_sp1=SP1_res['Qw']

            SP2_res = self.SP2.calc()
            Qw_sp2=SP2_res['Qw']

            # Баланс по сетевой воде
            Qw_summ = (
                (
                    self.water_streams.at["OTB2-SP2", "G"]
                    * (
                        self.water_streams.at["OTB2-SP2", "H"]
                        - self.water_streams.at["OD-GPK", "H"]
                    )
                    + self.water_streams.at["OTB1-SP1", "G"]
                    * (
                        self.water_streams.at["OTB1-SP1", "H"]
                        - self.water_streams.at["OD-GPK", "H"]
                    )
                )
                * self.KPD_SP/1000
            )
            Cp = 4.187
            t_pw = self.water_streams.at["SP2-WOUT", "T"]
            Qw_all = Cp * (t_pw - t_ow) * Gsv_in / 1000
            Error_all = (Qw_summ - Qw_all) / Qw_all * 100
#             Qw_sp1_sp2=Cp * (t_pw - t_sp1in) * Gsv_in / 1000
#             t_sp2in=self.water_streams.at["SP1-SP2", "T"]
#             Qw_sp2=Cp * (t_pw - t_sp2in) * Gsv_in / 1000
            Qw_summ_sp1 = ((self.water_streams.at["OTB1-SP1", "G"]* (self.water_streams.at["OTB1-SP1", "H"]- self.water_streams.at["SP1-OD", "H"]))+self.water_streams.at["SP2-SP1", "G"]* (self.water_streams.at["SP2-SP1", "H"]- self.water_streams.at["SP1-OD", "H"]))* self.KPD_SP/1000
    
            Qw_summ_sp2 = (self.water_streams.at["OTB2-SP2", "G"]* (self.water_streams.at["OTB2-SP2", "H"]- self.water_streams.at["SP2-SP1", "H"])) * self.KPD_SP/1000
            Qw_summ_od = (self.water_streams.at["SP1-OD", "G"]* (self.water_streams.at["SP1-OD", "H"]- self.water_streams.at["OD-GPK", "H"])) * self.KPD_SP/1000
#             print(self.water_streams.at["OTB2-SP2", "H"])
#             Qw_sp2q=Qw_sp2
                
                
#             Qw_summq=Qw_od+Qw_sp1+Qw_sp2

            # print('Qw_summ',Qw_summ)
            # print('Qw_all',Qw_all)
            # print('Error_all',Error_all)
#             print('ALL_delt',Qw_summ-Qw_all)
#             print('SP2_delt',Qw_summ_sp2-Qw_sp2)
#             print('SP1_delt',Qw_summ_sp1-Qw_sp1)
#             print('OD_delt',Qw_summ_od-Qw_od)

            if abs(Error_all) < calctolerance:
                # print('Мощность теплофикационной установки',Qw_summ)
                # print('Погрешность определения тепловой мощности теплофикационной установки',Error_all)
                
                break
            if i==maxiterations-1:
                print('Достигнуто максимальное количество итераций')

        return {'SP1':SP1_res, 'SP2':SP2_res, 'OD':OD_res}



