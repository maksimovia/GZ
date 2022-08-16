import numpy as n
from scipy.optimize import root

class heatex:
    def __init__(self,stream11,stream12,stream21,stream22,KPD,calctolerance,gas0,gas1,water,calcmethod,gas_streams0,water_streams0,gas_streams,water_streams):
        self.KPD=KPD
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas_streams=gas_streams
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.stream11=stream11
        self.stream12=stream12
        self.stream21=stream21
        self.stream22=stream22
        self.gas=gas1
        self.water=water
        self.gas0=gas0
        self.H011 = gas_streams0.at[stream11,'H']
        self.H012 = gas_streams0.at[stream12,'H']
        self.H021 = water_streams0.at[stream21,'H']
        self.H022 = water_streams0.at[stream22,'H']
        self.G01  = gas_streams0.at[stream11,'G']
        self.G02  = water_streams0.at[stream21,'G']
        self.P01  = gas_streams0.at[stream11,'P']
        self.P021 = water_streams0.at[stream21,'P']
        self.P022 = water_streams0.at[stream22,'P']        
        self.Q0     = self.G01*(self.H011-self.H012)*self.KPD
        T011   = self.gas0.p_h(self.P01,self.H011)['T']
        T012   = self.gas0.p_h(self.P01,self.H012)['T']
        T021   = self.water.p_h(self.P021,self.H021)['T']
        T022   = self.water.p_h(self.P022,self.H022)['T']
        dTmin0 = min(T011-T022,T012-T021)
        dTmax0 = max(T011-T022,T012-T021)
        self.LMTD0  = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av  = (T011+T012)/2
        T02av  = (T021+T022)/2
        P02av  = (self.P021+self.P022)/2
        self.lambda01av = self.gas0.p_t(self.P01,T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01,T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01,T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01,T01av)['rho']
        self.ro02av = self.water.p_t(P02av,T02av)['rho']
        self.ro021  = self.water.p_q(self.P021,1)['rho']        

    def calc(self):
        H11  = self.gas_streams.at[self.stream11,'H']
        H21  = self.water_streams.at[self.stream21,'H']
        G1   = self.gas_streams.at[self.stream11,'G']
        G2   = self.water_streams.at[self.stream21,'G']
        P1   = self.gas_streams.at[self.stream11,'P']
        P21  = self.water_streams.at[self.stream21,'P']
        ro21   = self.water.p_q(P21,1)['rho']
        ddp    = (ro21/self.ro021)*((self.G02/G2)**2)
        P22    = P21 - ((self.P021-self.P022)/ddp)
        P12 = P1
        T21 = self.water.p_h(P21,H21)['T']
        T11 = self.gas.p_h(P1,H11)['T']
        def T12sved(T12):
            if T12<T21 or T12>T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1,T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = H21 + (Q/G2)
                T22 = self.water.p_h(P22,H22)['T']
                dTmin= min(T11-T22,T12-T21)
                dTmax= max(T11-T22,T12-T21)
                if dTmin<0 or dTmax<0 or dTmin==dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin)) 
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1,T1av)['L']
                Pr1av = self.gas.p_t(P1,T1av)['Prandtl']
                nu1av = self.gas.p_t(P1,T1av)['nu']
                ro1av = self.gas.p_t(P1,T1av)['rho']
                
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33)*(((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst=T11
        try:
            Tfirst=T12
        except:
            Tfirst=T11*0.9
          #  max(Tfirst,T21+5)
        sol = root(T12sved,max(Tfirst,T21+5), method=self.calcmethod, tol=self.calctolerance)
        T12=float(sol.x)
        H12 = self.gas.p_t(P1,T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = H21 + (Q/G2)
        T22 = self.water.p_h(P22,H22)['T']
        return [T12,P12,H12,G1,T22,P22,H22,G2,Q]
    def calculate_and_write(self):
        calculation_result=self.calc()
        
        

class heatexPEND:
    def __init__(self,stream11,stream12,stream21,stream22,KPD,calctolerance,gas0,gas1,water,calcmethod,gas_streams0,water_streams0,gas_streams,water_streams):
        self.KPD=KPD
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas_streams=gas_streams
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.stream11=stream11
        self.stream12=stream12
        self.stream21=stream21
        self.stream22=stream22
        self.gas=gas1
        self.water=water
        self.gas0=gas0
        self.H011 = gas_streams0.at[stream11,'H']
        self.H012 = gas_streams0.at[stream12,'H']
        self.H021 = water_streams0.at[stream21,'H']
        self.H022 = water_streams0.at[stream22,'H']
        self.G01  = gas_streams0.at[stream11,'G']
        self.G02  = water_streams0.at[stream21,'G']
        self.P01  = gas_streams0.at[stream11,'P']
        self.P021 = water_streams0.at[stream21,'P']
        self.P022 = water_streams0.at[stream22,'P']        
        self.Q0     = self.G01*(self.H011-self.H012)*self.KPD
        T011   = self.gas0.p_h(self.P01,self.H011)['T']
        T012   = self.gas0.p_h(self.P01,self.H012)['T']
        T021   = self.water.p_h(self.P021,self.H021)['T']
        T022   = self.water.p_h(self.P022,self.H022)['T']
        dTmin0 = min(T011-T022,T012-T021)
        dTmax0 = max(T011-T022,T012-T021)
        self.LMTD0  = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av  = (T011+T012)/2
        T02av  = (T021+T022)/2
        P02av  = (self.P021+self.P022)/2
        self.lambda01av = self.gas0.p_t(self.P01,T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01,T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01,T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01,T01av)['rho']
        self.ro02av = self.water.p_t(P02av,T02av)['rho']
        self.ro021  = self.water.p_q(self.P021,1)['rho']        

    def calc(self):
        H11  = self.gas_streams.at[self.stream11,'H']
        H21  = self.water_streams.at[self.stream21,'H']
        G1   = self.gas_streams.at[self.stream11,'G']
        G2   = self.water_streams.at[self.stream21,'G']
        P1   = self.gas_streams.at[self.stream11,'P']
        P21  = self.water_streams.at[self.stream21,'P']
        ro21   = self.water.p_q(P21,1)['rho']
        ddp    = (ro21/self.ro021)*((self.G02/G2)**2)
        P22    = P21 - ((self.P021-self.P022)/ddp)
        P12 = P1
        T21 = self.water.p_h(P21,H21)['T']
        T11 = self.gas.p_h(P1,H11)['T']
        def T12sved(T12):
            if T12<T21 or T12>T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1,T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = H21 + (Q/G2)
                T22 = self.water.p_h(P22,H22)['T']
                dTmin= min(T11-T22,T12-T21)
                dTmax= max(T11-T22,T12-T21)
                if dTmin<0 or dTmax<0 or dTmin==dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin)) 
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1,T1av)['L']
                Pr1av = self.gas.p_t(P1,T1av)['Prandtl']
                nu1av = self.gas.p_t(P1,T1av)['nu']
                ro1av = self.gas.p_t(P1,T1av)['rho']
                
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33)*(((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst=T11
        try:
            Tfirst=T12
        except:
            Tfirst=T11*0.9
          #  max(Tfirst,T21+5)
        sol = root(T12sved,T11*0.99, method=self.calcmethod, tol=self.calctolerance)
        T12=float(sol.x)
        H12 = self.gas.p_t(P1,T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = H21 + (Q/G2)
        T22 = self.water.p_h(P22,H22)['T']
        return [T12,P12,H12,G1,T22,P22,H22,G2,Q]
        
    
    
    
class evaporVD:
    def __init__(self, stream11, stream12,stream21,stream22,KPD,calctolerance,gas,gas0,water,calcmethod,gas_streams0,water_streams0,gas_streams,water_streams):
        self.KPD=KPD
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas_streams=gas_streams
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.stream11=stream11
        self.stream12=stream12
        self.stream21=stream21
        self.stream22=stream22
        self.gas=gas
        self.water=water
        self.gas0=gas0
        self.H011 = gas_streams0.at[stream11,'H']
        self.H012 = gas_streams0.at[stream12,'H']
        self.H021 = water_streams0.at[stream21,'H']
        self.H022 = water_streams0.at[stream22,'H']
        self.G01  = gas_streams0.at[stream11,'G']
        self.G02  = water_streams0.at[stream21,'G']
        self.P01  = gas_streams0.at[stream11,'P']
        self.P02  = water_streams0.at[stream21,'P']
        self.Q0     = self.G01*(self.H011-self.H012)*self.KPD
        T011   = self.gas0.p_h(self.P01,self.H011)['T']
        T012   = self.gas0.p_h(self.P01,self.H012)['T']
        T021   = self.water.p_h(self.P02,self.H021)['T']
        T022   = self.water.p_h(self.P02,self.H022)['T']
        T0np   = self.water.p_q(self.P02,1)['T']        
        dTmin0 = min(T011-T0np,T012-T0np)
        dTmax0 = max(T011-T0np,T012-T0np)
        self.LMTD0  = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av  = (T011+T012)/2
        self.lambda01av = self.gas0.p_t(self.P01,T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01,T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01,T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01,T01av)['rho']
      
    def calc(self):
        H11  = self.gas_streams.at[self.stream11,'H']
        H21  = self.water_streams.at[self.stream21,'H']
        G1   = self.gas_streams.at[self.stream11,'G']
        P1   = self.gas_streams.at[self.stream11,'P']
        P2   = self.water_streams.at[self.stream21,'P']
        T21 = self.water.p_h(P2,H21)['T']
        T11 = self.gas.p_h(P1,H11)['T']
        T11 = self.gas.p_h(P1,H11)['T']
        T21 = self.water.p_h(P2,H21)['T']
        def T12sved(T12):
            T12=float(T12)
            if T12<T21 or T12>T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1,T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = self.water.p_q(P2,1)['h']  
                T22 = self.water.p_h(P2,H22)['T']
                G2 = Q/(H22-H21)
                dTmin= min(T11-T22,T12-T21)
                dTmax= max(T11-T22,T12-T21)
                if dTmin<0 or dTmax<0 or dTmin==dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))            
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1,T1av)['L']
                Pr1av = self.gas.p_t(P1,T1av)['Prandtl']
                nu1av = self.gas.p_t(P1,T1av)['nu']
                ro1av = self.gas.p_t(P1,T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33)*(((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst=T11
        try:
            Tfirst=T12
        except :
            Tfirst=T11*0.9
        sol = root(T12sved, max(Tfirst,T21+5), method=self.calcmethod, tol=self.calctolerance)
        T12=float(sol.x)
        H12 = self.gas.p_t(P1,T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = self.water.p_q(P2,1)['h']
        T22 = self.water.p_h(P2,H22)['T']
        G2 = Q/(H22-H21)
        return [T12,P1,H12,G1,T22,P2,H22,G2,Q]       
    
    

    
class evaporND:
    def __init__(self, stream11, stream12,stream21,stream22,streamVD, KPD,calctolerance,gas,gas0,water,calcmethod,gas_streams0,water_streams0,gas_streams,water_streams):
        self.KPD=KPD
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas_streams=gas_streams
        self.water_streams=water_streams
        self.calcmethod=calcmethod
        self.calctolerance=calctolerance
        self.stream11=stream11 
        self.stream12=stream12
        self.stream21=stream21
        self.stream22=stream22
        self.streamVD=streamVD
        self.gas=gas
        self.water=water
        self.gas0=gas0
        self.H011 = gas_streams0.at[stream11,'H']
        self.H012 = gas_streams0.at[stream12,'H']
        self.H021 = water_streams0.at[stream21,'H']
        self.H022 = water_streams0.at[stream22,'H']
        self.G01  = gas_streams0.at[stream11,'G']
        self.G02  = water_streams0.at[stream21,'G']
        self.P01  = gas_streams0.at[stream11,'P']
        self.P02  = water_streams0.at[stream21,'P']
        self.Q0     = self.G01*(self.H011-self.H012)*self.KPD
        T011   = self.gas0.p_h(self.P01,self.H011)['T']
        T012   = self.gas0.p_h(self.P01,self.H012)['T']
        T021   = self.water.p_h(self.P02,self.H021)['T']
        T022   = self.water.p_h(self.P02,self.H022)['T']
        T0np   = self.water.p_q(self.P02,1)['T']
        dTmin0 = min(T011-T0np,T012-T0np)
        dTmax0 = max(T011-T0np,T012-T0np)
        self.LMTD0  = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av  = (T011+T012)/2
        self.lambda01av = self.gas0.p_t(self.P01,T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01,T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01,T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01,T01av)['rho']
       
    def calc(self):

        H11  = self.gas_streams.at[self.stream11,'H']
        H21  = self.water_streams.at[self.stream21,'H']
        G1   = self.gas_streams.at[self.stream11,'G']
        P1   = self.gas_streams.at[self.stream11,'P']
        P2   = self.water_streams.at[self.stream21,'P']
        T21 = self.water.p_h(P2,H21)['T']
        T11 = self.gas.p_h(P1,H11)['T']
        Dvd = self.water_streams.at[self.streamVD,'G']
        def T12sved(T12):
            if T12<T21 or T12>T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1,T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = self.water.p_q(P2,1)['h']  
                T22 = self.water.p_h(P2,H22)['T']
                Hvd =self.water.p_q(P2,0)['h']
                G2 = (Q-Dvd*(Hvd-H21))/(H22-H21)
                dTmin= min(T11-T22, T12-T22)
                dTmax= max(T11-T22, T12-T22)
                if dTmin<0 or dTmax<0 or dTmin==dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))            
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1,T1av)['L']
                Pr1av = self.gas.p_t(P1,T1av)['Prandtl']
                nu1av = self.gas.p_t(P1,T1av)['nu']
                ro1av = self.gas.p_t(P1,T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33)*(((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst=T11
        try:
            Tfirst=T12
        except:
            Tfirst=T11*0.9
        sol = root(T12sved, max(Tfirst,T21+5), method=self.calcmethod, tol=self.calctolerance)        
        T12=float(sol.x)
        H12 = self.gas.p_t(P1,T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = self.water.p_q(P2,1)['h']
        T22 = self.water.p_h(P2,H22)['T']
        Hvd =self.water.p_q(P2,0)['h']
        G2 = (Q- Dvd*(Hvd-H21))/(H22-H21)
        Tvd = self.water.p_q(P2,0)['T']
        Pvd = P2
        return [T12,P1,H12,G1,T22,P2,H22,G2,Q,Tvd,Pvd,Hvd,Dvd]