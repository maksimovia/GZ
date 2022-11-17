
# Импорт библиотек
import matplotlib
import time
import cotel
import Turboustanovka
import mat_properties as prop

# Основные константы
# Calcmethod = "hybr"
# KPD_PN = 0.8074
# KPD_KN = 0.75
# KPD_to = 0.99
# KPD_SP = 0.9
# Calctolerance = 10**-2
# Teplo = 1
# Maxiterations_cotel = 10
# Maxiterations_KU_TU = 20
# Maxiterations_turbine = 30


class ku_tu:
    def __init__(self, gas0, gas1,  water, gas_streams0, gas_streams, water_streams0, water_streams, heaters, electric, streamKU_VD, streamKU_ND, streamST_VD, streamST_ND, Calcmethod, KPD_SP, KPKN, KPD_to, KPDPN, steamVD_fraction_to_turbine=1, steamVD_to_turbine=0):

        self.TU = Turboustanovka.turboustanovka(
            water, water_streams0, water_streams, heaters, electric, KPD_SP, KPKN)
        self.Whole_cotel = cotel.cotel_all(KPD_to, KPDPN,  gas0, gas1, water, Calcmethod,
                                           gas_streams0, water_streams0, gas_streams, water_streams, heaters, electric)
        self.steamVD_to_turbine = steamVD_to_turbine

        self.streamKU_VD = streamKU_VD
        self.streamKU_ND = streamKU_ND
        self.streamST_VD = streamST_VD
        self.streamST_ND = streamST_ND
        self.steamVD_fraction_to_turbine = steamVD_fraction_to_turbine

        self.Vvd0 = 1/water.p_h(water_streams0.at[self.streamKU_VD, 'P'],
                                water_streams0.at[self.streamKU_VD, 'H'])['rho']
        self.Vnd0 = 1/water.p_h(water_streams0.at[self.streamKU_ND, 'P'],
                                water_streams0.at[self.streamKU_ND, 'H'])['rho']
        self.Gvd0 = water_streams0.at[self.streamKU_VD, 'G']
        self.Gnd0 = water_streams0.at[self.streamKU_ND, 'G']
        self.dPvd0 = water_streams0.at[self.streamKU_VD, 'P'] - \
            water_streams0.at[self.streamST_VD, 'P']
        self.dPnd0 = water_streams0.at[self.streamKU_ND, 'P'] - \
            water_streams0.at[self.streamST_ND, 'P']
        self.water_streams0 = water_streams0
        self.water_streams = water_streams
        self.gas_streams0 = gas_streams0
        self.gas_streams = gas_streams
        self.gas0 = gas0
        self.gas1 = gas1
        self.water = water
        

        # print(self.water_streams)
        # Первое приближение по давлению
        

        if self.water_streams0.at[self.streamKU_ND, 'P'] == self.water_streams.at[self.streamKU_ND, 'P']:
            G_gas1 = self.gas_streams.at["GTU-PEVD", "G"]
            G_gas0 = self.gas_streams0.at["GTU-PEVD", "G"]
            if isinstance(G_gas1, float):
                print(self.gas_streams)
                # print(G_gas1)
                # G_gas1=G_gas0
            g_gas = G_gas1/G_gas0
            Pvd_1 = 0.4538+7.7601*g_gas
            Pnd_1 = -0.0189+0.6885*g_gas
            self.water_streams.at[self.streamKU_ND, 'P'] = Pnd_1
            self.water_streams.at[self.streamKU_VD, 'P'] = Pvd_1
        # print(self.water_streams)

    def calculate(self, Teplo, Calctolerance, Maxiterations_KU_TU, Maxiterations_cotel, Maxiterations_turbine):
        Pnd_it = []
        Pvd_it = []
        start_time = time.time()
        Max_error_P = 2
        Max_error_G = 2
        Max_error = 2
        Error_nd_P = 2
        Error_vd_P = 2
        Calctolerance_new = max(10**-1, Calctolerance)
        Teplo_overflow = 0
        Maxiterations_cotel_new = min(3,Maxiterations_cotel)
        Maxiterations_cotel_tu_rashod_new=min(3,Maxiterations_KU_TU)
        # print('Teplo',Teplo)

        for i in range(Maxiterations_KU_TU):
            Maxiterations_cotel_tu_rashod=Maxiterations_KU_TU
            

            for j in range(Maxiterations_cotel_tu_rashod_new):

                # Расчет котла

                Cotel_result = self.Whole_cotel.calc(
                    Calctolerance_new, Maxiterations_cotel_new)
                
                Gnd1 = self.water_streams.at[self.streamKU_ND, "G"]
                Gnd2 = self.water_streams.at[self.streamST_ND, "G"]
                Gvd1 = self.water_streams.at[self.streamKU_VD, "G"]
                Gvd2 = self.water_streams.at[self.streamST_VD, "G"]

                # Перекидываем расходы с учетом расхода в другие места
                if self.steamVD_to_turbine == 0:
                    new_VD_massflow = Gvd1*self.steamVD_fraction_to_turbine
                else:
                    new_VD_massflow = self.steamVD_to_turbine

                self.water_streams.at[self.streamST_VD, 'G'] = new_VD_massflow
                self.water_streams.at[self.streamST_ND, 'G'] = Gnd1

                # Перекидываем энтальпии
                self.water_streams.at[self.streamST_VD,
                                      'H'] = self.water_streams.at[self.streamKU_VD, 'H']
                self.water_streams.at[self.streamST_ND,
                                      'H'] = self.water_streams.at[self.streamKU_ND, 'H']

                # В первм приближении для упрощения считаем конденсационный режим
                if i>0 and Teplo == 1:
                    teplofikacia = 1
                else:
                    teplofikacia = 0
                    
                   
                # print('Teplo',Teplo)
                # print('teplofikacia',teplofikacia)
                
                # Расчет турбины
                TU_res = self.TU.calculate(
                    teplofikacia, calcmethod="hybr", calctolerance=Calctolerance_new, maxiterations=Maxiterations_turbine
                )
                

                # Calctolerance_new = Calctolerance
                if i > 2:
                    Calctolerance_new = Calctolerance
                    Maxiterations_cotel_new = Maxiterations_cotel
                if i>5:
                    Maxiterations_cotel_tu_rashod_new= Maxiterations_KU_TU
#                     if i == 3 and j == 0:
#                         print('Переход к оригинальной точности расчета',
#                               Calctolerance)
#                         print('Переход к оригинальному количетсву итераций',
#                               Maxiterations_cotel)

                # точка смешения на входе в ГПК с ПКМ

                if self.water_streams.at["ST-GPK", "G"] > 0:
                    G_smesh_od = self.water_streams.at["SMESHOD-REC", "G"]
                    H_smesh_od = self.water_streams.at["SMESHOD-REC", "H"]
                    G_smesh_PKM = self.water_streams.at["ST-GPK", "G"]
                    H_smesh_PKM = self.water_streams.at["ST-GPK", "H"]
                    G_v_GPK = G_smesh_od+G_smesh_PKM
                    P_v_GPK = self.water_streams.at["SMESHOD-REC", "P"]
                    H_v_GPK = (G_smesh_od*H_smesh_od +
                               G_smesh_PKM*H_smesh_PKM)/G_v_GPK
                    T_v_GPK = self.water.p_h(P_v_GPK, H_v_GPK)["T"]
                    self.water_streams.loc["SMESH-GPK",
                                           "T":"G"] = T_v_GPK, P_v_GPK, H_v_GPK, G_v_GPK
                else:
                    self.water_streams.loc["SMESH-GPK",
                                           "T":"G"] = self.water_streams.loc["SMESHOD-REC", "T":"G"]

                ##################

                G_turb = self.water_streams.at["SMESHOD-REC", "G"]
                G_ku = self.water_streams.at["GPK-IND",
                                             "G"]-(Gvd1-new_VD_massflow)
                Error_water_G = abs((self.water_streams.at["SMESHOD-REC", "G"] -
                                     (self.water_streams.at["GPK-IND", "G"]-(Gvd1-new_VD_massflow)))/(self.water_streams.at["GPK-IND", "G"]-(Gvd1-new_VD_massflow))*100)
                Error_nd_G = abs((Gnd1 - Gnd2)/Gnd1*100)
                Error_vd_G = abs((new_VD_massflow - Gvd2) /
                                 (new_VD_massflow)*100)
                Max_error_G = max(Error_water_G, Error_nd_G, Error_vd_G)
                Max_error = max(Error_water_G, Error_nd_G,
                                Error_vd_G, Error_nd_P, Error_vd_P)

                if Error_water_G > 20 and j != Maxiterations_cotel_tu_rashod - 1:
                    Teplo_overflow = 1
                    print(f"Расход из турбины G: {G_turb}")
                    print(f"Расход в ГПК G: {G_ku}")
#                 if Error_water_G > 1 and Error_water_G < 20:
#                     print("Погрешность определения расхода выше допустимой")
#                     print(f"Расход из турбины: {G_turb}")
#                     print(f"Расход в ГПК: {G_ku}")
                # if abs(Max_error_G) < Calctolerance_new:
                    #                     print(
                    #                         "Максимальная погрешность определения расхода в КУ+ПТУ", Max_error_G)
                    # break
                if j == Maxiterations_cotel_tu_rashod - 1:
                    print(f"Достигнуто максимальное количество итераций расхода КУ+ПТУ: {Maxiterations_cotel_tu_rashod}")
                    print(f"Error_water_G: {Error_water_G}, Error_nd_G: {Error_nd_G}, Error_vd_G: {Error_vd_G}")

            # Переписываю давления
            P_turb_vd = self.water_streams.at[self.streamST_VD, 'P']
            P_turb_nd = self.water_streams.at[self.streamST_ND, 'P']
            Vvd1 = 1/self.water.p_h(self.water_streams.at[self.streamKU_VD, 'P'],
                                    self.water_streams.at[self.streamKU_VD, 'H'])['rho']
            Vnd1 = 1/self.water.p_h(self.water_streams.at[self.streamKU_ND, 'P'],
                                    self.water_streams.at[self.streamKU_ND, 'H'])['rho']
            dPvd = self.dPvd0 * (Gvd1*Vvd1/self.Gvd0 /
                                 self.Vvd0)**2*self.Vvd0/Vvd1
            dPnd = self.dPnd0*(Gnd1*Vnd1/self.Gnd0/self.Vnd0)**2*self.Vnd0/Vnd1
            P_kotel_vd = self.water_streams.at[self.streamKU_VD, 'P']
            P_kotel_nd = self.water_streams.at[self.streamKU_ND, 'P']
            P_kotel_nd_new = P_turb_nd+dPnd
            P_kotel_vd_new = P_turb_vd+dPvd
            self.water_streams.at[self.streamKU_ND,
                                  'P'] = (P_kotel_nd_new+P_kotel_nd)/2
            self.water_streams.at[self.streamKU_VD,
                                  'P'] = (P_kotel_vd_new+P_kotel_vd)/2

            # Закидываю давления в массив
            Pnd_it.append(
                round(self.water_streams.loc[self.streamKU_ND, 'P'], 5))
            Pvd_it.append(
                round(self.water_streams.loc[self.streamKU_VD, 'P'], 5))

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
                f"Время {i+1} итерации расчета КУ+ТУ:---  {round((time.time() - start_time), 1)} сек. ---")
            if Teplo_overflow == 1:
                print('Слишком большая теплофикационная мощность, расчет окончен.')
                print(
                    'Для правильного расчета необходимо повысить мощность ГТУ или уменьшить мощность теплофикации.')
            if abs(Max_error) < Calctolerance and Calctolerance_new == Calctolerance:
                print('Расчет КУ+ПТУ окончен.')
                print("Максимальная погрешность определения расходов при расчете КУ+ПТУ", Error_water_G)
                print('Pnd_it', Pnd_it)
                print('Pvd_it', Pvd_it)
                break

            if i == Maxiterations_KU_TU - 1:
                print("Достигнуто максимальное количество итераций давления КУ+ПТУ:", i+1)
                print('Pnd_it', Pnd_it)
                print('Pvd_it', Pvd_it)
