import numpy as n
from scipy.optimize import root
import time
import nasos

class heatex:
    def __init__(self, stream11, stream12, stream21, stream22, KPD, calctolerance, gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.calctolerance = calctolerance
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.gas = gas1
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P021 = water_streams0.at[stream21, 'P']
        self.P022 = water_streams0.at[stream22, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_h(self.P021, self.H021)['T']
        T022 = self.water.p_h(self.P022, self.H022)['T']
        dTmin0 = min(T011-T022, T012-T021)
        dTmax0 = max(T011-T022, T012-T021)
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        T02av = (T021+T022)/2
        P02av = (self.P021+self.P022)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']
        self.ro02av = self.water.p_t(P02av, T02av)['rho']
        self.ro021 = self.water.p_q(self.P021, 1)['rho']

    def calc(self):
        start_timeProp = time.time()
        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        G2 = self.water_streams.at[self.stream21, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P21 = self.water_streams.at[self.stream21, 'P']
        ro21 = self.water.p_q(P21, 1)['rho']
        ddp = (ro21/self.ro021)*((self.G02/G2)**2)
        P22 = P21 - ((self.P021-self.P022)/ddp)
        P12 = P1
        T21 = self.water.p_h(P21, H21)['T']
        T11 = self.gas.p_h(P1, H11)['T']

        print("Prop:--- %s сек. ---" %
              round((time.time() - start_timeProp), 1))

        start_timeIter = time.time()

        def T12sved(T12):
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = H21 + (Q/G2)
                T22 = self.water.p_h(P22, H22)['T']
                dTmin = min(T11-T22, T12-T21)
                dTmax = max(T11-T22, T12-T21)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']

                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
          #  max(Tfirst,T21+5)
        sol = root(T12sved, max(Tfirst, T21+5),
                   method=self.calcmethod, tol=self.calctolerance)
        T12 = float(sol.x)
        print("Iter:--- %s сек. ---" %
              round((time.time() - start_timeIter), 1))
        start_timeElse = time.time()

        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = H21 + (Q/G2)
        T22 = self.water.p_h(P22, H22)['T']
        print("Else:--- %s сек. ---" %
              round((time.time() - start_timeElse), 1))

        return {'Tg': T12, 'Pg': P12, 'Hg': H12, 'Gg': G1, 'Qg': Q, 'Tw': T22, 'Pw': P22, 'Hw': H22, 'Gw': G2}


class heatexPEND:
    def __init__(self, stream11, stream12, stream21, stream22, KPD, calctolerance, gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.calctolerance = calctolerance
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.gas = gas1
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P021 = water_streams0.at[stream21, 'P']
        self.P022 = water_streams0.at[stream22, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_h(self.P021, self.H021)['T']
        T022 = self.water.p_h(self.P022, self.H022)['T']
        dTmin0 = min(T011-T022, T012-T021)
        dTmax0 = max(T011-T022, T012-T021)
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        T02av = (T021+T022)/2
        P02av = (self.P021+self.P022)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']
        self.ro02av = self.water.p_t(P02av, T02av)['rho']
        self.ro021 = self.water.p_q(self.P021, 1)['rho']

    def calc(self):
        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        G2 = self.water_streams.at[self.stream21, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P21 = self.water_streams.at[self.stream21, 'P']
        ro21 = self.water.p_q(P21, 1)['rho']
        ddp = (ro21/self.ro021)*((self.G02/G2)**2)
        P22 = P21 - ((self.P021-self.P022)/ddp)
        P12 = P1
        T21 = self.water.p_h(P21, H21)['T']
        T11 = self.gas.p_h(P1, H11)['T']

        def T12sved(T12):
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = H21 + (Q/G2)
                T22 = self.water.p_h(P22, H22)['T']
                dTmin = min(T11-T22, T12-T21)
                dTmax = max(T11-T22, T12-T21)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']

                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
          #  max(Tfirst,T21+5)
        sol = root(T12sved, T11*0.99, method=self.calcmethod,
                   tol=self.calctolerance)
        T12 = float(sol.x)
        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = H21 + (Q/G2)
        T22 = self.water.p_h(P22, H22)['T']
        return [T12, P12, H12, G1, T22, P22, H22, G2, Q]


class evaporVD:
    def __init__(self, stream11, stream12, stream21, stream22, KPD, calctolerance, gas, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.calctolerance = calctolerance
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.gas = gas
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P02 = water_streams0.at[stream21, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_h(self.P02, self.H021)['T']
        T022 = self.water.p_h(self.P02, self.H022)['T']
        T0np = self.water.p_q(self.P02, 1)['T']
        dTmin0 = min(T011-T0np, T012-T0np)
        dTmax0 = max(T011-T0np, T012-T0np)
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']

    def calc(self):
        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P2 = self.water_streams.at[self.stream21, 'P']
        T21 = self.water.p_h(P2, H21)['T']
        T11 = self.gas.p_h(P1, H11)['T']
        T11 = self.gas.p_h(P1, H11)['T']
        T21 = self.water.p_h(P2, H21)['T']

        def T12sved(T12):
            T12 = float(T12)
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = self.water.p_q(P2, 1)['h']
                T22 = self.water.p_h(P2, H22)['T']
                G2 = Q/(H22-H21)
                dTmin = min(T11-T22, T12-T21)
                dTmax = max(T11-T22, T12-T21)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
        sol = root(T12sved, max(Tfirst, T21+5),
                   method=self.calcmethod, tol=self.calctolerance)
        T12 = float(sol.x)
        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = self.water.p_q(P2, 1)['h']
        T22 = self.water.p_h(P2, H22)['T']
        G2 = Q/(H22-H21)
        return [T12, P1, H12, G1, T22, P2, H22, G2, Q]


class evaporND:
    def __init__(self, stream11, stream12, stream21, stream22, streamVD, KPD, calctolerance, gas, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams):
        self.KPD = KPD
        self.gas_streams0 = gas_streams0
        self.water_streams0 = water_streams0
        self.gas_streams = gas_streams
        self.water_streams = water_streams
        self.calcmethod = calcmethod
        self.calctolerance = calctolerance
        self.stream11 = stream11
        self.stream12 = stream12
        self.stream21 = stream21
        self.stream22 = stream22
        self.streamVD = streamVD
        self.gas = gas
        self.water = water
        self.gas0 = gas0
        self.H011 = gas_streams0.at[stream11, 'H']
        self.H012 = gas_streams0.at[stream12, 'H']
        self.H021 = water_streams0.at[stream21, 'H']
        self.H022 = water_streams0.at[stream22, 'H']
        self.G01 = gas_streams0.at[stream11, 'G']
        self.G02 = water_streams0.at[stream21, 'G']
        self.P01 = gas_streams0.at[stream11, 'P']
        self.P02 = water_streams0.at[stream21, 'P']
        self.Q0 = self.G01*(self.H011-self.H012)*self.KPD
        T011 = self.gas0.p_h(self.P01, self.H011)['T']
        T012 = self.gas0.p_h(self.P01, self.H012)['T']
        T021 = self.water.p_h(self.P02, self.H021)['T']
        T022 = self.water.p_h(self.P02, self.H022)['T']
        T0np = self.water.p_q(self.P02, 1)['T']
        dTmin0 = min(T011-T0np, T012-T0np)
        dTmax0 = max(T011-T0np, T012-T0np)
        self.LMTD0 = (dTmax0 - dTmin0) / (n.log(dTmax0/dTmin0))
        T01av = (T011+T012)/2
        self.lambda01av = self.gas0.p_t(self.P01, T01av)['L']
        self.Pr01av = self.gas0.p_t(self.P01, T01av)['Prandtl']
        self.nu01av = self.gas0.p_t(self.P01, T01av)['nu']
        self.ro01av = self.gas0.p_t(self.P01, T01av)['rho']

    def calc(self):

        H11 = self.gas_streams.at[self.stream11, 'H']
        H21 = self.water_streams.at[self.stream21, 'H']
        G1 = self.gas_streams.at[self.stream11, 'G']
        P1 = self.gas_streams.at[self.stream11, 'P']
        P2 = self.water_streams.at[self.stream21, 'P']
        T21 = self.water.p_h(P2, H21)['T']
        T11 = self.gas.p_h(P1, H11)['T']
        Dvd = self.water_streams.at[self.streamVD, 'G']

        def T12sved(T12):
            if T12 < T21 or T12 > T11:
                return 10**9
            else:
                H12 = self.gas.p_t(P1, T12)['h']
                Q = G1*(H11-H12)*self.KPD
                H22 = self.water.p_q(P2, 1)['h']
                T22 = self.water.p_h(P2, H22)['T']
                Hvd = self.water.p_q(P2, 0)['h']
                G2 = (Q-Dvd*(Hvd-H21))/(H22-H21)
                dTmin = min(T11-T22, T12-T22)
                dTmax = max(T11-T22, T12-T22)
                if dTmin < 0 or dTmax < 0 or dTmin == dTmax:
                    LMTD = (dTmax+dTmin)/2
                else:
                    LMTD = (dTmax - dTmin) / (n.log(dTmax/dTmin))
                dt = self.LMTD0/LMTD
                T1av = (T11+T12)/2
                lambda1av = self.gas.p_t(P1, T1av)['L']
                Pr1av = self.gas.p_t(P1, T1av)['Prandtl']
                nu1av = self.gas.p_t(P1, T1av)['nu']
                ro1av = self.gas.p_t(P1, T1av)['rho']
                kk = (self.lambda01av/lambda1av)*((self.Pr01av/Pr1av)**0.33) * \
                    (((self.G01/G1)*(ro1av/self.ro01av)*(nu1av/self.nu01av))**0.685)
                Qq = self.Q0 / (kk*dt)
                return ((Q-Qq)/Q)*100
        Tfirst = T11
        try:
            Tfirst = T12
        except:
            Tfirst = T11*0.9
        sol = root(T12sved, max(Tfirst, T21+5),
                   method=self.calcmethod, tol=self.calctolerance)
        T12 = float(sol.x)
        H12 = self.gas.p_t(P1, T12)['h']
        Q = G1*(H11-H12)*self.KPD
        H22 = self.water.p_q(P2, 1)['h']
        T22 = self.water.p_h(P2, H22)['T']
        Hvd = self.water.p_q(P2, 0)['h']
        G2 = (Q - Dvd*(Hvd-H21))/(H22-H21)
        Tvd = self.water.p_q(P2, 0)['T']
        Pvd = P2
        return [T12, P1, H12, G1, T22, P2, H22, G2, Q, Tvd, Pvd, Hvd, Dvd]


class cotel_all:
    def __init__(self, KPD,KPDnasos, calctolerance, gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams):
        self.PEVD_obj = heatex('GTU-PEVD', 'PEVD-IVD', 'IVD-PEVD', 'PEVD-DROSVD',
                          KPD, calctolerance, gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams)

        self.IVD_obj = evaporVD('PEVD-IVD', 'IVD-EVD', 'EVD-IVD', 'IVD-PEVD',
                           KPD, calctolerance, gas1, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams)

        self.EVD_obj = heatex('IVD-EVD', 'EVD-PPND', 'PEN-EVD', 'EVD-IVD',
                         KPD, calctolerance, gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams)

        self.PEN_obj = nasos.nasos('BND-PEN', 'PEN-EVD', water, KPDnasos, water_streams)

        self.PPND_obj = heatexPEND('EVD-PPND', 'PPND-IND', 'IND-PPND', 'PPND-DROSND',
                              KPD, calctolerance, gas0, gas1, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams)

        self.IND_obj = evaporND('PPND-IND', 'IND-GPK', 'GPK-IND', 'IND-PPND',  'PEVD-DROSVD',
                           KPD, calctolerance, gas1, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams)

        self.GPK_obj = heatex('IND-GPK', 'GPK-out', 'X-GPK', 'GPK-IND',
                         KPD, calctolerance, gas1, gas0, water, calcmethod, gas_streams0, water_streams0, gas_streams, water_streams)
        self.KPD=KPD
        self.calctolerance=calctolerance
        self.water_streams=water_streams
        self.gas_streams=gas_streams
        self.gas_streams0=gas_streams0
        self.water_streams0=water_streams0
        self.gas1=gas1

       

    def calc(self,maxiterations=50):

        it = maxiterations
        start_time = time.time()

        for i in range(it):
            # Связвка высокого давления
            for j in range(it):

                start_timePEVD = time.time()

                # Расчёт ПЕВД
                PEVD = self.PEVD_obj.calc()
                self.gas_streams.loc['PEVD-IVD', 'T':'G'] = [PEVD['Tg'],
                                                        PEVD['Pg'], PEVD['Hg'], PEVD['Gg']]
                self.water_streams.loc['PEVD-DROSVD', 'T':'G'] = [PEVD['Tw'],
                                                     PEVD['Pw'], PEVD['Hw'], PEVD['Gw']]
                
                print("PEVD:--- %s сек. ---" %
                      round((time.time() - start_timePEVD), 1))

                start_timeIVD = time.time()

                # Расчёт ИВД
                IVD = self.IVD_obj.calc()
                self.gas_streams.loc['IVD-EVD', 'T':'G'] = [IVD[0], IVD[1], IVD[2], IVD[3]]
                self.water_streams.loc['IVD-PEVD','T':'G'] = [IVD[4], IVD[5], IVD[6], IVD[7]]

                print("IVD:--- %s сек. ---" % round((time.time() - start_timeIVD), 1))

                # Переопределение расхода в ВД
                self.water_streams.loc['PEVD-DROSVD':'PEN-EVD', 'G'] = IVD[7]
                self.water_streams.loc['BND-PEN', 'G'] = IVD[7]

                start_timeEVD = time.time()
                
                # Расчёт ЭВД
                EVD = self.EVD_obj.calc()
                print("EVD:--- %s сек. ---" % round((time.time() - start_timePEVD), 1))
                
                self.gas_streams.loc['EVD-PPND', 'T':'G'] = [EVD['Tg'],
                                                EVD['Pg'], EVD['Hg'], EVD['Gg']]
                self.water_streams.loc['EVD-IVD', 'T':'G'] = [EVD['Tw'],
                                                 EVD['Pw'], EVD['Hw'], EVD['Gw']]

                # Баланс ПЕВ+ИВД+ЭВД
                Qgas1 = self.KPD*self.gas_streams.at['GTU-PEVD', 'G'] * \
            (self.gas_streams.at['GTU-PEVD', 'H']-self.gas_streams.at['IVD-EVD', 'H'])
                Qwat1 = self.water_streams.at['PEVD-DROSVD', 'G']*(\
            self.water_streams.at['PEVD-DROSVD', 'H']-self.water_streams.at['EVD-IVD', 'H'])
                print('dQ/Q ПЕВД+ИВД+ЭВД', ((Qgas1-Qwat1)/Qgas1*100))
                if abs((Qgas1-Qwat1)/Qgas1*100) < self.calctolerance:
                    break
            print("ВД: --- %s сек. ---" % round((time.time() - start_time), 1))
            # Для сходимости
            if i == 0:
                self.gas_streams.loc['PPND-IND', 'T'] = self.gas_streams.loc['EVD-PPND', 'T'] - 1
                self.gas_streams.loc['PPND-IND', 'H'] = self.gas1.p_t(self.gas_streams.loc['PPND-IND', 'P'], self.gas_streams.loc['PPND-IND', 'T'])['h']
                
            # Связка низкого давления
            for j in range(it):
                # Расчёт ППНД
                PPND = self.PPND_obj.calc()
                self.gas_streams.loc['PPND-IND', 'T':'G'] = [PPND[0],
                                                PPND[1], PPND[2], PPND[3]]
                self.water_streams.loc['PPND-DROSND',
                          'T':'G'] = [PPND[4], PPND[5], PPND[6], PPND[7]]

                # Расчёт ИНД
                IND = self.IND_obj.calc()
                self.gas_streams.loc['IND-GPK', 'T':'G'] = [IND[0], IND[1], IND[2], IND[3]]
                self.water_streams.loc['IND-PPND',
                          'T':'G'] = [IND[4], IND[5], IND[6], IND[7]]

                # Переопределение расхода в НД
                self.water_streams.loc['PPND-DROSND':'IND-PPND', 'G'] = IND[7]

                # ПЭН
                self.water_streams.loc['BND-PEN', 'T':'G'] = [IND[9],
                                                 IND[10], IND[11], IND[12]]
                PEN = self.PEN_obj.calc()
                self.water_streams.loc['PEN-EVD',
                          'T':'G'] = [PEN[0], PEN[1], PEN[2], PEN[3]]

                # Баланс ППНД+ИНД
                Qgas = self.KPD*self.gas_streams.at['EVD-PPND', 'G'] * \
            (self.gas_streams.at['EVD-PPND', 'H']-self.gas_streams.at['IND-GPK', 'H'])
                Qwat = self.water_streams.at['IND-PPND', 'G']*(self.water_streams.at['IND-PPND', 'H']-self.water_streams.at['GPK-IND', 'H']) +\
            self.water_streams.at['PPND-DROSND', 'G']*(self.water_streams.at['PPND-DROSND', 'H'] -\
                                                  self.water_streams.at['IND-PPND', 'H'])+\
        self.water_streams.at['BND-PEN', 'G']*(self.water_streams.at['BND-PEN', 'H']-self.water_streams.at['GPK-IND', 'H'])
                
                # Расчёт расхода в ГПК
                self.water_streams.loc['X-GPK', 'G'] = self.water_streams.at['PPND-DROSND',
                                                           'G']+self.water_streams.at['PEVD-DROSVD', 'G']
                self.water_streams.loc['GPK-IND', 'G'] = self.water_streams.at['X-GPK', 'G']
                
                start_timeGPK = time.time()
                # Расчёт ГПК
                GPK = self.GPK_obj.calc()
                self.gas_streams.loc['GPK-out', 'T':'G'] = [GPK['Tg'],
                                               GPK['Pg'], GPK['Hg'], GPK['Gg']]
                self.water_streams.loc['GPK-IND', 'T':'G'] = [GPK['Tw'],
                                                 GPK['Pw'], GPK['Hw'], GPK['Gw']]
                print("GPK:--- %s сек. ---" % round((time.time() - start_timeGPK), 1))

                # Баланс ППНД+ИНД+ГПК
                Qgas1 = self.KPD*self.gas_streams.at['EVD-PPND', 'G'] * \
            (self.gas_streams.at['EVD-PPND', 'H']-self.gas_streams.at['GPK-out', 'H'])
                Qwat1 = self.water_streams.at['GPK-IND', 'G']*(self.water_streams.at['GPK-IND', 'H']-self.water_streams.at['X-GPK', 'H']) +\
            self.water_streams.at['IND-PPND', 'G']*(self.water_streams.at['IND-PPND', 'H']-self.water_streams.at['GPK-IND', 'H']) +\
            self.water_streams.at['PPND-DROSND', 'G']*(self.water_streams.at['PPND-DROSND', 'H']-self.water_streams.at['IND-PPND', 'H']) +\
            self.water_streams.at['BND-PEN', 'G'] * \
            (self.water_streams.at['BND-PEN', 'H']-self.water_streams.at['GPK-IND', 'H'])
                print('dQ/Q ППНД+ИНД+ГПК', (Qgas1-Qwat1)/Qgas1*100)
                if abs((Qgas1-Qwat1)/Qgas1*100) < self.calctolerance:
                    break
            print("НД+ --- %s сек. ---" % round((time.time() - start_time), 1))

            # Баланс общий
            Qgasall = self.KPD*self.gas_streams.at['GTU-PEVD', 'G'] * \
        (self.gas_streams.at['GTU-PEVD', 'H']-self.gas_streams.at['GPK-out', 'H'])
            Qwatall = self.water_streams.at['PPND-DROSND', 'G']*(self.water_streams.at['PPND-DROSND', 'H']-self.water_streams.at['X-GPK', 'H'])+self.water_streams.at['PEVD-DROSVD', 'G']*(
        self.water_streams.at['PEVD-DROSVD', 'H']-self.water_streams.at['X-GPK', 'H'])-self.water_streams.at['BND-PEN', 'G']*(self.water_streams.at['PEN-EVD', 'H']-self.water_streams.at['BND-PEN', 'H'])
            print('dQ/Qsumm', (Qgasall-Qwatall)/Qgasall*100)
            if abs((Qgasall-Qwatall)/Qgasall) < self.calctolerance:
                break
