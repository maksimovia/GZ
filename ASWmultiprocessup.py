def ParallelCompute(args):
              # Импорт библиотек
 
    import os
    import time
    import ASW
    import GTU
    import KU_TU
    import mat_properties as prop
    import numpy as n
    import pandas as pd
    import SP
    from scipy.optimize import root
    import calculate_all as calc
    # таблица номинального режима
    gas_streams0 = pd.read_excel("streams0.xlsx", sheet_name="gas", index_col=0)
    water_streams0 = pd.read_excel("streams0.xlsx", sheet_name="water", index_col=0)
    GTU_ISO = pd.read_excel("input.xlsx", sheet_name="ISO", index_col=0)
    GTU_input = pd.read_excel("input.xlsx", sheet_name="GTU_input", index_col=0)
    # рабочая таблица (=номинал в 1 итерации)
    gas_streams = pd.read_excel("streams.xlsx", sheet_name="gas", index_col=0)
    water_streams = pd.read_excel("streams.xlsx", sheet_name="water", index_col=0)
    # рабочая таблица показателей блоков
    heaters = pd.read_excel("blocks.xlsx", sheet_name="heaters", index_col=0)
    electric = pd.read_excel("blocks.xlsx", sheet_name="electric", index_col=0)
    ############################################################
    # Теплосеть и перекидка температуры воздуха
    gas_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    water_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    Tnv = gas_streams.at["AIR", "T"]
    water_streams.at["SWINB", "T"] = SP.Tset(Tnv)[1]
    water_streams.at["SWOUT", "T"] = SP.Tset(Tnv)[0]
    water_streams.at["SWIN-TURB", "T"] = water_streams.at["SWIN", "T"]
    water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"]
    water_streams.at["SP2-WOUT", "T"] = water_streams.at["SWOUT", "T"]
    ############################################################


    # Основные эффективности оборудования
    KPD_PN = 0.8074
    KPD_KN = 0.75
    KPD_to = 0.99
    KPD_SP = 0.99


    # Параметры, отвечающие за процесс расчета
    Calcmethod = "hybr"
    Calctolerance = 10**-2
    Maxiterations_KU_TU = 20
    Maxiterations_cotel = 4
    Maxiterations_turbine = 30

    # Параметры режима работы ПГУ
    # Расчет для минимума нагрузки
    Сalculate_minimum = True
    # Расчет для работы с теплофикацией
    Teplo = int(True)
    # Отбор пара высокго давления или доля или кг/с
    steamVD_fraction_to_turbine = 1


    # Состав газов в номинале в ГТУ
    gasmix = "Nitrogen*Oxygen*CO2*Water*Argon"
    # Считывание рефпропа
    RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")
    fractiongas0 = list(gas_streams0.loc["GTU-PEVD", "N2":"Ar"])
    gas0 = prop.Materials_prop(
        gasmix,
        fractiongas0,
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
    )
    water = prop.Materials_prop(
        "water",
        [1.0, 0, 0, 0, 0],
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
    )

    #####################
    ######Максимов#######
    fractionwaterMethane = (0.833372660622383, 0.166627339377617, 0, 0, 0)
    waterMethanemix = "Water*METHANE"

    waterMethane = prop.Materials_prop(
        waterMethanemix,
        fractionwaterMethane,
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
    )
    Methane = prop.Materials_prop(
        "METHANE",
        [1.0, 0, 0, 0, 0],
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
    )
    gas_KU_PKM = prop.Materials_prop(
        gasmix,
        (
        0.710320591016015,
        0.00996710270335893,
        0.090538556815177,
        0.180531273012258,
        0.00864247645319178,
    ),
        prop.REFPROP_h_s,
        prop.REFPROP_p_t,
        prop.REFPROP_p_h,
        prop.REFPROP_p_s,
        prop.REFPROP_p_q,
        prop.REFPROP_t_q,
        prop.REFPROP_p_rho,
        prop.REFPROP_s_q,
        RP=RP,
    )

    # Вкл выкл заряд по ПКМ
    PKM_zaryad = False
    #####################
    ASW_bull = 1

    # Задание энтальпий газа в номинальном режиме
    Temperatures = gas_streams0.loc["GTU-KU":"GPK-out", "T"]
    Pressure = gas_streams0.loc["GTU-KU", "P"]
    Enthalpies = list(map(lambda x: gas0.p_t(Pressure, x)["h"], Temperatures))


    gas_streams0.loc["GTU-KU":"GPK-out", "H"] = Enthalpies

    ############################################################
    # Задание ГТУ
    Gas_turbine = GTU.gtu(GTU_ISO, "GTU-KU")

    # Расчет всей ПГУ вместе
    def calculate_CCGT(
        Iterations_KU_TU,
        Iterations_cotel,
        Iterations_turbine,):

        # Расчет ГТУ
        Gas_turbine_res = Gas_turbine.calc(GTU_input)

        # Запись данных об электричестве
        electric.at["GTU", "N"] = Gas_turbine_res["N"]
        electric.at["GTU", "KPD"] = Gas_turbine_res["eff"]
        electric.at["DK", "N"] = Gas_turbine_res["Ndk"]

        # Запись данных о газе на выходе из ГТУ
        gas_streams = pd.read_excel("streams.xlsx", sheet_name="gas", index_col=0)
        gas_streams.at["GTU-KU", "T"] = Gas_turbine_res["T"]
        gas_streams.at["GTU-KU", "G"] = Gas_turbine_res["G"]
        gas_streams.at["GTU-KU", "P"] = 0.1
        gas_streams.at["GTU-KU", "H"] = gas0.p_t(
            gas_streams.at["GTU-KU", "P"], gas_streams.at["GTU-KU", "T"]
        )["h"]
        Gas_turbine_composition = pd.read_excel(
            "input.xlsx", sheet_name="Gas_composition0", index_col=0
        )

        #####################пкм
        ######Максимов#######
        if PKM_zaryad == True:

            steamVD_to_turbine = 0.25*water_streams0.at["PEVD-DROSVD", "G"]   ##!!!максимов
            from PKM import steam_transformer
            # Пар в паротрансформатор
            water_streams.loc["DROSVD-ST", "T":"H"] = water_streams.loc["PEVD-DROSVD", "T":"H"]
            water_streams.loc["DROSVD-ST", "G"] = water_streams.at["PEVD-DROSVD", "G"] - steamVD_to_turbine  ##!!!!
            # паротрансформатор
            ST = steam_transformer(
                stream11="DROSVD-ST",
                water=water,
                water_streams=water_streams,
                heaters=heaters,
                Pdr1=2,Pdr2=0.8,P2=2,dT=15,dTmin=5,Tdec=10)
            steam_trans = ST.calc()

            # Ввод в табл выходов из паротрансформатора
            water_streams.loc["ST-GPK", "T":"G"] = [steam_trans["T17"],steam_trans["P17"],steam_trans["H17"],steam_trans["G1"]]
            water_streams.loc["ST-PKM", "T":"G"] = [steam_trans["T24"],steam_trans["P2"],steam_trans["H24"],steam_trans["G2"]]
            heaters.at["Strans", "Qw"] = steam_trans["Q"]
            heaters.at["Strans_Qcool", "Qw"] = steam_trans['Qcool80']

            print(heaters.at["Strans_Qcool", "Qw"])
            # реформер
            from PKM import reformer

            ref = reformer(
                stream11="ST-PKM",
                water=water,
                gas_KU=gas_KU_PKM,
                Methane=Methane,
                waterMethane=waterMethane,
                water_streams=water_streams,
                heaters=heaters,
                Tref=700,
                Pref=2,
                T1gas=1968.58395330148,
                T2gas=800,
            )
            reform = ref.calc()
            # Газы реформера
            gas_streams.loc["AIR-REF", "T":"G"] = [15, 0.1, 414.38, reform["Gair"]]
            gas_streams.loc["CH4-REF", "T":"G"] = [15, 0.7, 881.50, reform["Gch4"]]
            gas_streams.loc["REF-SMESH", "T":"G"] = [800, 0.1, reform["H2gas"], reform["Ggas"]]
            gas_streams.loc["REF-SMESH", "N2":"Ar"] = list(reform["Gasfrac"].values())

            # Смешение
            gas_streams.loc["GTU-PEVD", "G"] = (gas_streams.at["REF-SMESH", "G"] + gas_streams.at["GTU-KU", "G"])
            gas_streams.loc["GTU-PEVD", "H"] = (gas_streams.at["REF-SMESH", "G"] * gas_streams.at["REF-SMESH", "H"]
                                                + gas_streams.at["GTU-KU", "G"] * gas_streams.at["GTU-KU", "H"]) / gas_streams.loc["GTU-PEVD", "G"]
            gas_streams.loc["GTU-PEVD", "P"] = 0.1

            from PKM import mixing_gases_molar

            mixing_gases_molar("GTU-KU", "REF-SMESH", "GTU-PEVD", gas_streams)
            for stream in gas_streams.index[4:10]:
                gas_streams.loc[stream, "N2":"Ar"] = gas_streams.loc["GTU-PEVD", "N2":"Ar"]
        else:
            gas_streams.loc["GTU-PEVD", "T":"Ar"] = gas_streams.loc["GTU-KU", "T":"Ar"]
            water_streams.loc["ST-GPK", "T":"G"] = [80,2,320,0]
            steamVD_to_turbine=0
    # ----------------------------- конец Максимова 
    # ----------------------------- Начало Максима
        ASW = ASW.Accum(water,water_streams,accumulation,
                         stream12 = "ASW-WOUT",
                         stream11 = "SP2-ASW",
                         stream_obratnoi_setevoi_vody = "SWIN-TURB",
                         stream_pryamoi_setevoi_vody = "SP2-WOUT",
                         T_nar_vozd = water_streams.at["AIR", "T"])

        ASW.set_construct(Diametr=Diametr,
                           kolichestvo = kolichestvo,
                           Visota = Visota,
                           lambda_min_vata = lambda_min_vata,
                           delta_min_vata = delta_min_vata)

        start_time = time.time()



            # # если Зарядка 
        if ASWbul == 1:
            G_ASW_zarydka = ASW.zaryadka(vremya)['G']
            water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka
            water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] + G_ASW_zarydka
                ## если разрядка
        if ASWbul == 2:
            ASW.zaryadka(vremya)
            ASW.jdat(vremyaojidaniya)
            G_ASW_razryadka = ASW.razryadka(vremya)['G']
            water_streams.at["SWIN-TURB", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
            water_streams.at["SP2-WOUT", "G"] = water_streams.at["SWIN", "G"] - G_ASW_razryadka
            water_streams.at["SWOUT", "H"] = water.p_t(water_streams.at["SWOUT","P"],water_streams.at["SWOUT", "T"])['h']
            water_streams.at["SP2-WOUT", "H"] = (water_streams.at["SWOUT", "H"]*water_streams.at["SWOUT", "G"] -
                                                water_streams.at["ASW-WOUT", "H"]*G_ASW_razryadka)/(water_streams.at["SWOUT", "G"]-G_ASW_razryadka)
            water_streams.at["SP2-WOUT", "T"] = water.p_h(water_streams.at["SP2-WOUT", "P"],water_streams.at["SP2-WOUT", "H"])['T']


    #     # Параметры газа на входе в КУ
    #     gas_streams.loc["GTU-PEVD", "T":"P"] = gas_streams.loc["GTU-KU", "T":"P"]
    #     gas_streams.at["GTU-PEVD", "G"] = gas_streams.loc["GTU-KU", "G"]
    #     gas_streams.loc["GTU-PEVD", "N2":"Ar"] = Gas_turbine_composition.loc["Fraction", "N2":"Ar"]

        # Состав газов при частичной нагрузке
        fractiongas = list(gas_streams.loc["GTU-PEVD", "N2":"Ar"])

        gas1 = prop.Materials_prop(
            gasmix,
            fractiongas,
            prop.REFPROP_h_s,
            prop.REFPROP_p_t,
            prop.REFPROP_p_h,
            prop.REFPROP_p_s,
            prop.REFPROP_p_q,
            prop.REFPROP_t_q,
            prop.REFPROP_p_rho,
            prop.REFPROP_s_q,
            RP=RP,
        )

        # Инициализаця KU+TU, она здесь потому что нжно менять состав газа на входе в КУ

        KU_and_TU = KU_TU.ku_tu(
            gas0,
            gas1,
            water,
            gas_streams0,
            gas_streams,
            water_streams0,
            water_streams,
            heaters,
            electric,
            "PEVD-DROSVD",
            "PPND-DROSND",
            "DROSVD-TURBVD",
            "DROSND-TURBND",
            Calcmethod,
            KPD_SP,
            KPD_KN,
            KPD_to,
            KPD_PN,
            steamVD_fraction_to_turbine,
            steamVD_to_turbine,
        )

        start_time = time.time()

        # Расчет КУ и ТУ
        KU_and_TU.calculate(
            Teplo,
            Calctolerance,
            Iterations_KU_TU,
            Iterations_cotel,
            Iterations_turbine,
        )

        print(f"fin КУ и ТУ:--- {round((time.time() - start_time), 1)} сек. ---")


    calculate_CCGT(Maxiterations_KU_TU, Maxiterations_cotel, Maxiterations_turbine)
    GTU_input = pd.read_excel("input.xlsx", sheet_name="GTU_input", index_col=0)

    Max_iterations_minimum = 10

    if Сalculate_minimum == True:
        n_GTU = GTU_input.at["n", 1]
        start_time = time.time()
        for i in range(Max_iterations_minimum):
            print(f"Началась {i+1} итерация расчета ПГУ")
            if i < round(Max_iterations_minimum/2,1):
                New_iterations_KU_TU, New_iterations_cotel, New_iterations_turbine = (
                    2,
                    2,
                    15,
                )
            else:
                New_iterations_KU_TU, New_iterations_cotel, New_iterations_turbine = (
                    Maxiterations_KU_TU,
                    Maxiterations_cotel,
                    Maxiterations_turbine,
                )

            Delt_Gcnd = 100
            Delt_Nturb = 100
            if Teplo == 1:
                Delt_Gcnd = (water_streams.at["INKOND", "G"] - 4.44) / 4.44
                Delt_Nturb = (electric.at["Turbine", "Ni"] - 17.6) / 17.6
            Delt_Gvd = (
                water_streams.at["PEVD-DROSVD", "G"] / water_streams0.at["PEVD-DROSVD", "G"]
                - 0.25
            ) / 0.25
            Delt_Gnd = (
                water_streams.at["PPND-DROSND", "G"] / water_streams0.at["PPND-DROSND", "G"]
                - 0.5
            ) / 0.5
            Delta_min = min(Delt_Gcnd, Delt_Nturb, Delt_Gvd, Delt_Gnd)
            if n_GTU == 1 and Delta_min < 0:
                print("Мощность ГТУ 100% и расход пара все еще слишком мал")
            n_GTU = n_GTU - Delta_min / 5
            GTU_input.at["n", 1] = n_GTU

            calculate_CCGT(
                New_iterations_KU_TU, New_iterations_cotel, New_iterations_turbine
            )
            print(f"Отклонение от ограничения минимальное равно {Delta_min}")
            Delta_n_GTU = abs((n_GTU_it[-1] - n_GTU_it[-2]) / n_GTU_it[-1] * 100)
            if abs(Delta_min) < Calctolerance and Delta_n_GTU < Calctolerance:
                calculate_CCGT(
                New_iterations_KU_TU, New_iterations_cotel, New_iterations_turbine
            )
                print(f"Отклонение от ограничения минимальное равно {Delta_min}")
                print(f"Относительная мощность ГТУ равна {n_GTU}")
                print(
                    f"fin минимальная мощность ПГУ:--- {round((time.time() - start_time), 1)} сек. ---"
                )
                break
            if i == Max_iterations_minimum - 1:
                print("Достигнуто максимальное количество итераций минимального расхода в ПГУ", i+1)
    
#     Перепад на дифрагму для проверки параметров работы программы при данных исходных данных
    P_pered_diafragma=water_streams.at["DOOTB1", "P"]
    P_posle_diafragma=water_streams.at["INCND", "P"]
    Delta_P_Diafragma=P_pered_diafragma-P_posle_diafragma
    
    result = {
    "T_air":round(t_air_list,2),
    "n_GTU":round(nagr,2),
    "GTU": round(electric.at["GTU", "N"], 4),
    "GTU_KPD": round(electric.at["GTU", "KPD"], 4),
    "Turbine": round(electric.at["Turbine", "Ni"], 4),
    "KN": round(electric.at["KN", "Ni"], 4),
    "DK": round(electric.at["DK", "N"], 4),
    "PEN": round(electric.at["PEN", "Ni"], 4),
    "Turbine_Qt":round(heaters.at["SP2", "Qw"]+heaters.at["SP1", "Qw"]+heaters.at["OD", "Qw"], 4),
    "ASW_Qt":round(accumulation.at["ASW", "Qw"]/(vremya*3600), 4),
    "Delta_P_Diafragma":round(Delta_P_Diafragma, 4)
    }


    return result
    

