
# Импорт библиотек
import matplotlib
import matplotlib.pyplot as plt
import time
import cotel
import Turboustanovka
import mat_properties as prop

# Основные константы
Calcmethod = "hybr"
KPD_PN = 0.8074
KPD_KN = 0.75
KPD_to = 0.99
KPD_SP = 0.9
Calctolerance = 10**-2
Teplo = 1
Maxiterations_cotel = 10
Maxiterations_KU_TU = 20
Maxiterations_turbine = 30


class ku_tu:
    def __init__(self, gas0, gas1,  water, gas_streams0, gas_streams, water_streams0, water_streams, heaters, electric, Calcmethod, KPD_SP, KPKN, KPD_to, KPDPN):

        self.TU = Turboustanovka.turboustanovka(
            water, water_streams0, water_streams, heaters, electric, KPD_SP, KPKN)
        self.Whole_cotel = cotel.cotel_all(KPD_to, KPDPN,  gas0, gas1, water, Calcmethod,
                                           gas_streams0, water_streams0, gas_streams, water_streams, heaters, electric)

        self.Vvd0 = 1/water.p_h(water_streams0.at["PEVD-DROSVD", 'P'],
                           water_streams0.at["PEVD-DROSVD", 'H'])['rho']
        self.Vnd0 = 1/water.p_h(water_streams0.at["PPND-DROSND", 'P'],
                           water_streams0.at["PPND-DROSND", 'H'])['rho']
        self.Gvd0 = water_streams0.at["PEVD-DROSVD", 'G']
        self.Gnd0 = water_streams0.at["PPND-DROSND", 'G']
        self.dPvd0 = water_streams0.at["PEVD-DROSVD", 'P'] - \
            water_streams0.at["DROSVD-TURBVD", 'P']
        self.dPnd0 = water_streams0.at["PPND-DROSND", 'P'] - \
            water_streams0.at["DROSND-TURBND", 'P']
        self.water_streams0 = water_streams0
        self.water_streams = water_streams
        self.gas_streams0 = gas_streams0
        self.gas_streams = gas_streams
        self.gas0 = gas0
        self.gas1 = gas1
        self.water = water

    def calculate(self, Teplo, Calctolerance, Maxiterations_KU_TU, Maxiterations_cotel, Maxiterations_turbine):
        Pnd_it = []
        Pvd_it = []
        start_time = time.time()
        Max_error_P = 2
        Max_error_G = 2
        Max_error = 2
        Error_nd_P = 2
        Error_vd_P = 2
        Calctolerance_new = 10**-1
        Teplo_overflow = 0

        # Первое приближение по давлению
        G_gas1 = self.gas_streams.at["GTU-PEVD", "G"]
        G_gas0 = self.gas_streams0.at["GTU-PEVD", "G"]
        g_gas = G_gas1/G_gas0
        Pvd_1 = 0.4538+7.7601*g_gas
        Pnd_1 = -0.0189+0.6885*g_gas
        self.water_streams.at["PPND-DROSND", 'P'] = Pnd_1
        self.water_streams.at["PEVD-DROSVD", 'P'] = Pvd_1

        for i in range(Maxiterations_KU_TU):

            for j in range(Maxiterations_cotel):

                # Расчет котла
                Cotel_result = self.Whole_cotel.calc(
                    Calctolerance_new, Maxiterations_cotel)

                Gnd1 = self.water_streams.at["PPND-DROSND", "G"]
                Gnd2 = self.water_streams.at["DROSND-TURBND", "G"]
                Gvd1 = self.water_streams.at["PEVD-DROSVD", "G"]
                Gvd2 = self.water_streams.at["DROSVD-TURBVD", "G"]

                # Перекидываем расходы
                self.water_streams.at["DROSVD-TURBVD", 'G'] = Gvd1
                self.water_streams.at["DROSND-TURBND", 'G'] = Gnd1

                # Перекидываем энтальпии
                self.water_streams.at["DROSVD-TURBVD",
                                      'H'] = self.water_streams.at["PEVD-DROSVD", 'H']
                self.water_streams.at["DROSND-TURBND",
                                      'H'] = self.water_streams.at["PPND-DROSND", 'H']

                # В первм приближении для упрощения считаем конденсационный режим
                if i > 0 and Teplo == 1:
                    teplofikacia = 1
                else:
                    teplofikacia = 0

                # Расчет турбины
                TU_res = self.TU.calculate(
                    teplofikacia, calcmethod="hybr", calctolerance=Calctolerance_new, maxiterations=Maxiterations_turbine
                )

              
                Calctolerance_new = Calctolerance*10
                if i > 2:
                    Calctolerance_new = Calctolerance
                    if i==3 and j==0:
                        print('Переход к оригинальной точности расчета', Calctolerance)
            
                Error_water_G = abs((self.water_streams.at["SMESHOD-REC", "G"] -
                                     self.water_streams.at["GPK-IND", "G"])/self.water_streams.at["GPK-IND", "G"]*100)
                Error_nd_G = abs((Gnd1 - Gnd2)/Gnd1*100)
                Error_vd_G = abs((Gvd1 - Gvd2)/Gvd1*100)
                Max_error_G = max(Error_water_G, Error_nd_G, Error_vd_G)
                Max_error = max(Error_water_G, Error_nd_G,
                                Error_vd_G, Error_nd_P, Error_vd_P)
               
                if Error_water_G > 20:
                    Teplo_overflow = 1
                if Error_water_G > 1:
                    print(self.water_streams)
                if abs(Max_error_G) < Calctolerance_new:
                    print(
                        "Максимальная погрешность определения расхода в КУ+ПТУ", Max_error_G)
                    break
                if j == Maxiterations_cotel - 1:
                    print("Достигнуто максимальное количество итераций расхода КУ+ПТУ")

            # Переписываю давления
            P_turb_vd = self.water_streams.at["DROSVD-TURBVD", 'P']
            P_turb_nd = self.water_streams.at["DROSND-TURBND", 'P']
            Vvd1 = 1/self.water.p_h(self.water_streams.at["PEVD-DROSVD", 'P'],
                                    self.water_streams.at["PEVD-DROSVD", 'H'])['rho']
            Vnd1 = 1/self.water.p_h(self.water_streams.at["PPND-DROSND", 'P'],
                                    self.water_streams.at["PPND-DROSND", 'H'])['rho']
            dPvd = self.dPvd0 * (Gvd1*Vvd1/self.Gvd0/self.Vvd0)**2*self.Vvd0/Vvd1
            dPnd = self.dPnd0*(Gnd1*Vnd1/self.Gnd0/self.Vnd0)**2*self.Vnd0/Vnd1
            P_kotel_vd = self.water_streams.at["PEVD-DROSVD", 'P']
            P_kotel_nd = self.water_streams.at["PPND-DROSND", 'P']
            P_kotel_nd_new = P_turb_nd+dPnd
            P_kotel_vd_new = P_turb_vd+dPvd
            self.water_streams.at["PPND-DROSND",
                                  'P'] = (P_kotel_nd_new+P_kotel_nd)/2
            self.water_streams.at["PEVD-DROSVD",
                                  'P'] = (P_kotel_vd_new+P_kotel_vd)/2

            # Закидываю давления в массив
            Pnd_it.append(round(self.water_streams.loc['PPND-DROSND', 'P'], 5))
            Pvd_it.append(round(self.water_streams.loc['PEVD-DROSVD', 'P'], 5))
            print('Pnd_it', Pnd_it)
            print('Pvd_it', Pvd_it)

            # Ошибки расчета
            Error_nd_P = abs((P_kotel_nd - P_kotel_nd_new)/P_kotel_nd*100)
            Error_vd_P = abs((P_kotel_vd - P_kotel_vd_new)/P_kotel_vd*100)
            Max_error = max(Error_water_G, Error_nd_G,
                            Error_vd_G, Error_nd_P, Error_vd_P)
            # print('Error_nd_P', Error_nd_P)
            # print('Error_vd_P', Error_vd_P)
            Max_error_P = max(Error_nd_P, Error_vd_P)
            # print('Max_error', Max_error)
            print(
                f"Время {i+1} итерации расчета КУ+ТУ:--- %s сек. --- {round((time.time() - start_time), 1)}")
            if Teplo_overflow == 1:
                print('Слишком большая теплофикационная мощность, расчет окончен.')
                print(
                    'Для правильного расчета необходимо повысить мощность ГТУ или уменьшить мощность теплофикации.')
            if abs(Max_error) < Calctolerance and Calctolerance_new == Calctolerance:
                print(
                    "Максимальная погрешность определения расходов", Error_water_G)
                break
            if i == Maxiterations_KU_TU - 1:
                print("Достигнуто максимальное количество итераций давления КУ+ПТУ", i+1)
