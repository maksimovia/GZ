import mat_properties as prop
import numpy as n
import pandas as pd
import SP
import Turbine
from scipy.optimize import root
import nasos
import time


class turboustanovka:
    def __init__(self, water, water_streams0, water_streams, heaters, electric, KPD_SP=0.99, KPDnasos=0.8):
        self.Tepl_systema = SP.teplofik_systema(
            KPD_SP, water, water_streams0, water_streams
        )

        self.Turb = Turbine.turbine(
            water,
            "DROSVD-TURBVD",
            "ENDOFVD",
            "DROSND-TURBND",
            "SMESHEND",
            "DOOTB2",
            "DOOTB1",
            "INCND",
            "INKOND",
            water_streams0,
            water_streams,
        )
        self.water_streams = water_streams
        self.KPD_SP = KPD_SP
        self.heaters = heaters
        self.electric = electric
        self.water = water

        self.KN = nasos.nasos('KOND-KN', 'KN-GPK', water,
                              KPDnasos, water_streams, water_streams0)

    def Find_Potb2(self, Diafragma):
        Diafragma_it = [0]
        G_sp2_it = 1
        G_sp1_it = 1
        P_turb_it = 1
        P_tepl_it = 1
        for i in range(self.maxiterations):
            Tepl_systema_res = self.Tepl_systema.calculate(
                self.maxiterations, calctolerance
            )
            G_sp2 = Tepl_systema_res["SP2"]["Gotb"]
            P_sp2 = Tepl_systema_res["SP2"]["p_otb"]
            G_sp1 = Tepl_systema_res["SP1"]["Gotb"]
            G_sp2 = max(0.01, G_sp2)
            G_sp1 = max(0.01, G_sp1)
            G_ots3 = max(
                0.01, self.water_streams.at["DOOTB2", "G"] - G_sp2)
            self.water_streams.at["DOOTB1", "G"] = G_ots3
            G_CND = max(
                0.001, self.water_streams.at["DOOTB1", "G"] - G_sp1)
            self.water_streams.loc["INCND":"INKOND", "G"] = G_CND
            self.water_streams.loc["OTB1-SP1", "T":"H"] = self.water_streams.loc[
                "DOOTB1", "T":"H"
            ]
            self.water_streams.at["OTB2-SP2", "T"] = self.water_streams.at[
                "DOOTB2", "T"
            ]
            self.water_streams.at["OTB2-SP2", "H"] = self.water_streams.at[
                "DOOTB2", "H"
            ]
            Turb_res = self.Turb.calculate(
                Diafragma, self.maxiterations, calctolerance)
            Potb2_turb = self.water_streams.at["DOOTB2", "P"]
            Potb2_teplof = self.water_streams.at["OTB2-SP2", "P"]

            Error_Gsp2 = (G_sp2_it-G_sp2)/G_sp2_it*100
            Error_Gsp1 = (G_sp1_it-G_sp1)/G_sp1_it*100
            Error_P_turb = (P_turb_it-Potb2_turb)/P_turb_it*100
            Error_P_tepl = (P_tepl_it-Potb2_teplof)/P_tepl_it*100
            G_sp2_it = G_sp2
            G_sp1_it = G_sp1
            P_turb_it = Potb2_turb
            P_tepl_it = Potb2_teplof

            Max_error = max(Error_Gsp2, Error_Gsp1, Error_P_turb, Error_P_tepl)
            # print("Potb2_turb", Potb2_turb)
            print("G_sp2", G_sp2)
            print("G_sp2_it", G_sp2_it)
            print("G_sp1", G_sp1)
            print("G_sp1_it", G_sp1_it)

            # print("Diafragma", Diafragma)
            if abs(Max_error) < calctolerance:
                print("G_sp2", G_sp2, "G_sp1", G_sp1, "P_sp2", P_sp2)
                print(
                    "Выход из цикла сведения турбины и теплофикации",
                    Max_error,
                )
                break
            if i == self.maxiterations - 1:
                print(
                    "Достигнуто максимальное количество итераций нахождения давления в верхнем отборе")

        Error_p = (Potb2_turb - Potb2_teplof) / Potb2_teplof * 100
        Diafragma_it.append(Diafragma)
        print("Error_p", Error_p)
        print("Diafragma", Diafragma)
        return Error_p

    def Find_Potb2_it(self, maxiterations, calctolerance):
        Diafragma = 0.05
        Diafragma_it = [0]
        for i in range(maxiterations):
            G_sp2_it = 1
            G_sp1_it = 1
            P_turb_it = 1
            P_tepl_it = 1
            for j in range(maxiterations):
                Tepl_systema_res = self.Tepl_systema.calculate(
                    maxiterations, calctolerance
                )
                G_sp2 = Tepl_systema_res["SP2"]["Gotb"]
                P_sp2 = Tepl_systema_res["SP2"]["p_otb"]
                G_sp1 = Tepl_systema_res["SP1"]["Gotb"]
                # print("G_sp2", G_sp2, "G_sp1", G_sp1, "P_sp2", P_sp2)
                G_sp2 = max(0.01, G_sp2)
                G_sp1 = max(0.01, G_sp1)
                G_ots3 = max(
                    0.001, self.water_streams.at["DOOTB2", "G"] - G_sp2)
                self.water_streams.at["DOOTB1", "G"] = G_ots3
                G_CND = max(
                    0.001, self.water_streams.at["DOOTB1", "G"] - G_sp1)
                self.water_streams.loc["INCND":"INKOND", "G"] = G_CND
                self.water_streams.loc["OTB1-SP1", "T":"H"] = self.water_streams.loc[
                    "DOOTB1", "T":"H"
                ]
                self.water_streams.at["OTB2-SP2", "T"] = self.water_streams.at[
                    "DOOTB2", "T"
                ]
                self.water_streams.at["OTB2-SP2", "H"] = self.water_streams.at[
                    "DOOTB2", "H"
                ]
                Turb_res = self.Turb.calculate(
                    Diafragma, maxiterations, calctolerance)
                Potb2_turb = self.water_streams.at["DOOTB2", "P"]
                Potb2_teplof = self.water_streams.at["OTB2-SP2", "P"]

                Error_Gsp2 = (G_sp2_it-G_sp2)/G_sp2_it*100
                Error_Gsp1 = (G_sp1_it-G_sp1)/G_sp1_it*100
                Error_P_turb = (P_turb_it-Potb2_turb)/P_turb_it*100
                Error_P_tepl = (P_tepl_it-Potb2_teplof)/P_tepl_it*100
                G_sp2_it = G_sp2
                G_sp1_it = G_sp1
                P_turb_it = Potb2_turb
                P_tepl_it = Potb2_teplof

                Max_error = max(Error_Gsp2, Error_Gsp1,
                                Error_P_turb, Error_P_tepl)
                # print("G_sp2", G_sp2)
                # print("G_sp2_it", G_sp2_it)
                # print("G_sp1", G_sp1)
                # print("G_sp1_it", G_sp1_it)

                if abs(Max_error) < calctolerance:
                    # print(
                    #     "Выход из цикла сведения турбины и теплофикации",
                    #     Max_error,
                    # )
                    break
                if i == maxiterations - 1:
                    print("Достигнуто максимальное количество итераций расход и давлений в турбине при работе с теплофикацией")

            Error_p = (Potb2_turb - Potb2_teplof) / Potb2_teplof * 100
            Diafragma = max(0, Diafragma - Error_p / 2500)
            Diafragma_it.append(Diafragma)
            # print("Diafragma", Diafragma)
            # print("Potb2_turb", Potb2_turb)
            if abs(Error_p) < calctolerance/10:
                print(
                    "Максимальная погрешность определения давления в верхнем отборе",
                    Error_p,
                )
                break
            if i == maxiterations - 1:
                print(
                    "Достигнуто максимальное количество итераций давления верхнего отбора"
                )
        # print(Error_p)
        # print(Diafragma_it)
        return Diafragma

    def calculate_t_rejim(self, calcmethod, calctolerance, maxiterations):

        Diafragma = self.Find_Potb2_it(maxiterations, calctolerance)
        # self.maxiterations=maxiterations
        # Solution_Potb2 = root(self.Find_Potb2, 0.08,
        #                                       method=calcmethod, tol=calctolerance)
        # Diafragma = float(Solution_Potb2.x)

        Turb_res = self.Turb.calculate(Diafragma, maxiterations, calctolerance)
        Tepl_systema_res = self.Tepl_systema.calculate(
            maxiterations, calctolerance)

        Result = {
            "Turb_res": self.Turb.calculate_power(),
            "Tepl_systema_res": Tepl_systema_res,
            "Diafragma": Diafragma,
        }
        return Result

    def calculate(
        self,
        teplofikacia=0,
        calcmethod="hybr",
        calctolerance=10**-3,
        maxiterations=20,
    ):

        start_time = time.time()

        if teplofikacia == 1:
            Result = self.calculate_t_rejim(
                calcmethod, calctolerance, maxiterations)
            self.water_streams.at['OTB2-SP2','G']=Result["Tepl_systema_res"]['SP2']['Gotb']
            self.water_streams.at['OTB1-SP1','G']=Result["Tepl_systema_res"]['SP1']['Gotb']
            Qsp1 = Result['Tepl_systema_res']['SP1']['Qw']
            Qsp2 = Result['Tepl_systema_res']['SP2']['Qw']
            Qod = Result['Tepl_systema_res']['OD']['Qw']
            self.heaters.at["SP2", "Qw"] = Qsp2
            self.heaters.at["SP1", "Qw"] = Qsp1
            self.heaters.at["OD", "Qw"] = Qod
            self.heaters.at["SP2", "KPD"] = self.KPD_SP
            self.heaters.at["SP1", "KPD"] = self.KPD_SP
            self.heaters.at["OD", "KPD"] = self.KPD_SP

        if teplofikacia == 0:
            diafragma = 0
            self.water_streams.loc["DOOTB1":"INKOND", "G"] = self.water_streams.at["DROSVD-TURBVD",
                                                                                   'G']+self.water_streams.at["DROSND-TURBND", 'G']
            Turb_res = self.Turb.calculate(
                diafragma, maxiterations, calctolerance)
            Result = {"Turb_res": self.Turb.calculate_power()}

        # запись данных в таблицу блоков

        Ni = Result['Turb_res']['Ni']
        Ni1 = Result['Turb_res']['Ni1']
        Ni2 = Result['Turb_res']['Ni2']
        Ni3 = Result['Turb_res']['Ni3']
        Ni4 = Result['Turb_res']['Ni4']

        self.electric.at["Turbine", "Ni"] = Ni
        self.electric.loc["Tots1":"Tots4", "Ni"] = [Ni1, Ni2, Ni3, Ni4]

        # Расчет насоса конденсатного (расход должен меняться из расчета - добавочная вода в конденсатор)
        pk = self.water_streams.loc["INKOND", 'P']
        KONDout = self.water.p_q(pk, 0)
        T_KNin = KONDout['T']
        H_KNin = KONDout['h']
        self.water_streams.loc["KOND-KN", "T":'H'] = [T_KNin, pk, H_KNin]
        self.water_streams.at["KOND-KN",
                              "G"] = self.water_streams.at["INKOND", 'G']
        self.water_streams.at["KN-GPK",
                              "P"] = self.water_streams.at["SMESHOD-REC", "P"]
        KN_res = self.KN.calc()
        Result['KN'] = KN_res
        self.water_streams.loc["KN-GPK", 'T':'G'] = [KN_res['T2'],
                                                     KN_res['P2'], KN_res['h2real'], KN_res['G1']]
        self.electric.loc['KN', 'Ni':'KPD'] = [KN_res['Ni'],
                                               KN_res['Ngm'], KN_res['KPDm'], KN_res['KPD']]
        # Расчет смешения перед ГПК
        h_kn = self.water_streams.at["KN-GPK", "H"]
        G_kn = self.water_streams.at["KN-GPK", "G"]

        h_od = self.water_streams.at["OD-GPK", "H"]
        G_od = self.water_streams.at["OD-GPK", "G"]
        if teplofikacia == 0:
            G_od = 0
            h_od = 0
        G_smeshod = G_kn+G_od
        h_smeshod = (h_od*G_od+h_kn*G_kn)/G_smeshod
        p_smeshod = self.water_streams.at["SMESHOD-REC", "P"]
        t_smeshod = self.water.p_h(p_smeshod, h_smeshod)['T']
        self.water_streams.at["SMESHOD-REC", "T"] = t_smeshod
        self.water_streams.at["SMESHOD-REC", "H"] = h_smeshod
        self.water_streams.at["SMESHOD-REC", "G"] = G_smeshod
        # print('G_smeshod',G_smeshod,'h_smeshod',h_smeshod,'t_smeshod',t_smeshod,)
        print("Fin турбоустановка:--- %s сек. ---" %
              round((time.time() - start_time), 2))

        return Result
