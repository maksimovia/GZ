# Импорт библиотек

def ParallelCompute_PKM(air_temperature):
    import os
    import time

    import GTU
    import KU_TU
    import mat_properties as prop
    import numpy as n
    import pandas as pd
    import SP
    from calculate_CCGT_PKM import Calculate_CCGT_PKM_iter
    from scipy.optimize import root

    # таблица номинального режима
    gas_streams0 = pd.read_excel(
        "streams0.xlsx", sheet_name="gas", index_col=0)
    water_streams0 = pd.read_excel(
        "streams0.xlsx", sheet_name="water", index_col=0)
    GTU_ISO = pd.read_excel("input.xlsx", sheet_name="ISO", index_col=0)
    GTU_input = pd.read_excel(
        "input.xlsx", sheet_name="GTU_input", index_col=0)
    # рабочая таблица (=номинал в 1 итерации)
    GTU_input.at["tair", 1] = air_temperature
    gas_streams = pd.read_excel("streams.xlsx", sheet_name="gas", index_col=0)
    water_streams = pd.read_excel(
        "streams.xlsx", sheet_name="water", index_col=0)
    # рабочая таблица показателей блоков
    heaters = pd.read_excel("blocks.xlsx", sheet_name="heaters", index_col=0)
    electric = pd.read_excel("blocks.xlsx", sheet_name="electric", index_col=0)
    accumulation = pd.read_excel(
        "blocks.xlsx", sheet_name="accumulation", index_col=0)
    ############################################################
    # Теплосеть и перекидка температуры воздуха
    gas_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    water_streams.loc["AIR", "T":"P"] = [GTU_input.loc["tair", 1], 0.1]
    Tnv = gas_streams.at["AIR", "T"]
    water_streams.at["SWIN", "T"] = SP.Tset(Tnv)[1]
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
    Maxiterations_KU_TU = 5
    Maxiterations_cotel = 3
    Maxiterations_turbine = 15
    Iter_pkm = 8

    # Параметры режима работы ПГУ
    # Расчет для минимума нагрузки
    Сalculate_minimum = False

    # Расчет для работы с теплофикацией
    # Teplo = int(False)
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
    # Задание энтальпий газа в номинальном режиме
    Temperatures = gas_streams0.loc["GTU-KU":"GPK-out", "T"]
    Pressure = gas_streams0.loc["GTU-KU", "P"]
    gas_streams0.loc["GTU-KU":"GPK-out", "H"] = list(
        map(lambda x: gas0.p_t(Pressure, x)["h"], Temperatures)
    )

    #####################Максимов#####################
    water_streams.at["SWIN", "H"] = water.p_t(
        1, water_streams.at["SWIN", "T"])["h"]
    water_streams.at["SWOUT", "H"] = water.p_t(
        1, water_streams.at["SWOUT", "T"])["h"]

    # Время зарядки/разрядки, часы
    time_ac = 4
    # Время ожидания, часы
    time_jdat = 12
    # Конструкция аккумулятора
    constr = {
        "Diametr": 20,
        "kolichestvo": 8,
        "Visota": 20,
        "lambda_min_vata": 0.045,
        "delta_min_vata": 0.01,
    }
    PKM_zaryad = True
    PKM_razryad = False
    syngas_streams = pd.read_excel(
        "streams.xlsx", sheet_name="syngas", index_col=0)
    ##################################################

    ############################################################
    # Задание ГТУ
    Gas_turbine = GTU.gtu(GTU_ISO, "GTU-KU")

    pkm_pgu_tol = 10**-2

    arguments_all_it = [
        Maxiterations_KU_TU,
        Maxiterations_cotel,
        Maxiterations_turbine,
        gas_streams0,
        water_streams0,
        GTU_ISO,
        GTU_input,
        gas_streams,
        water_streams,
        heaters,
        electric,
        Gas_turbine,
        gas0,
        water,
        PKM_zaryad,
        PKM_razryad,
        syngas_streams,
        Calcmethod,
        Calctolerance,
        KPD_PN,
        KPD_KN,
        KPD_to,
        KPD_SP,
        steamVD_fraction_to_turbine,
        accumulation,
        time_ac,
        constr,
        time_jdat,
    ]

    ########################ОГРАНИЧЕНИЯ НА РАБОТУ ПГУ########################
    # ЕСЛИ ТЕПЛОФИКАЦИЯ (из документов по ПГУ-220Т)
    # Gк_мин=4,44 кг/с вроде
    # Nтурбины мин = 17,6 МВт (примерно 25%)
    # Максимальноа давление в отборах СП 0,245, 0,198 МПа (Теплофикационная паровая турбина Т-63/76-8.8 для серии ПГУ-230)
    # ДЛЯ ВСЕХ РЕЖИМОВ (Из Трухния по ПГУ-450Т)
    # Gвд_мин=25% от номинала
    # Gнд_мин = 50% от номинала

    ########РАСЧЕТ расхода пара на ПКМ из условия заполнения хранилища################
    ########РАСЧЕТ мощности ГТУ из условия тепловой мощности и ################

    ################Расчет минимальной нагрузки ГТУ при остальных нормальных условиях############

    Max_iterations_minimum = 20
    
    start_time_all = time.time()

    n_GTU_it = [0.5]
    Delta_n_GTU = 100
    coeficient_PGU = 5
    if Сalculate_minimum == True:
        gas_streams.loc["GTU-PEVD", "G"] = gas_streams.loc["GTU-KU", "G"]
        n_GTU = GTU_input.at["n", 1]
        start_time = time.time()
        Delta_min = 0
        n_GTU_it.append(round(n_GTU, 5))
        for i in range(Max_iterations_minimum):
            print("n_GTU:", n_GTU_it)
            print("Delta_n_GTU: ", Delta_n_GTU)
            print("Delta_min: ", Delta_min)
            if i<7 :#Delta_n_GTU > 1 :
                (
                    New_iterations_KU_TU,
                    New_iterations_cotel,
                    New_iterations_turbine,
                    New_Iter_pkm,
                    New_coeficient_PGU
                ) = (3, 2, 15, 4, 3)
            else:
                # print("Delta_n_GTU: ", Delta_n_GTU)
                (
                    New_iterations_KU_TU,
                    New_iterations_cotel,
                    New_iterations_turbine,
                    New_Iter_pkm,
                    New_coeficient_PGU
                ) = (
                    Maxiterations_KU_TU,
                    Maxiterations_cotel,
                    Maxiterations_turbine,
                    Iter_pkm,
                    coeficient_PGU
                )

            gas_streams = Calculate_CCGT_PKM_iter(
                arguments_all_it, New_Iter_pkm, pkm_pgu_tol
            )
            Delt_Gcnd = 100
            Delt_Nturb = 100
            Delt_Gcnd = (water_streams.at["INKOND", "G"] - 4.44) / 4.44
            Delt_Nturb = (electric.at["Turbine", "Ni"] - 17.6) / 17.6
            Delt_Gvd = (
                water_streams.at["PEVD-DROSVD", "G"] /
                water_streams0.at["PEVD-DROSVD", "G"]
                - 0.25
            ) / 0.25
            Delt_Gnd = (
                water_streams.at["PPND-DROSND", "G"] /
                water_streams0.at["PPND-DROSND", "G"]
                - 0.5
            ) / 0.5

            arguments_all_it[0], arguments_all_it[1], arguments_all_it[2] = (
                New_iterations_KU_TU,
                New_iterations_cotel,
                New_iterations_turbine,
            )

            Delta_min = min(Delt_Gcnd, Delt_Nturb, Delt_Gvd, Delt_Gnd)
            n_GTU = n_GTU - Delta_min / New_coeficient_PGU
            if n_GTU>1:
                print("Мощность больше 1 у ГТУ")
                n_GTU=min(n_GTU,1.1)
                
            if n_GTU<0:
                print("Мощность меньше 0 у ГТУ")
                n_GTU=0.6
            Delta_n_GTU = abs(
                (n_GTU_it[-1] - n_GTU_it[-2]) / n_GTU_it[-1] * 100)
            GTU_input.at["n", 1] = n_GTU
            n_GTU_it.append(round(n_GTU, 5))

            print(
                f"Время {i+1} итерации расчета мощности ГТУ при ПГУ с ПКМ: --- {round((time.time() - start_time), 1)} сек. ---"
            )
            print(f"Отклонение от ограничения минимальное равно {Delta_min}")

            if n_GTU == 1 and Delta_min < 0:
                print("Мощность ГТУ 100% и расход пара все еще слишком мал")

            if abs(Delta_min) < Calctolerance and Delta_n_GTU < Calctolerance/10:
                arguments_all_it[0], arguments_all_it[1], arguments_all_it[2] = (
                    Maxiterations_KU_TU,
                    Maxiterations_cotel,
                    Maxiterations_turbine,
                )
                n_GTU_it.append(round(n_GTU, 5))
                gas_streams = Calculate_CCGT_PKM_iter(
                    arguments_all_it, New_Iter_pkm, pkm_pgu_tol
                )
                print(
                    f"Отклонение от ограничения минимальное равно {Delta_min}")
                print(f"Относительная мощность ГТУ равна {n_GTU}")
                print(
                    f"fin минимальная мощность ПГУ:--- {round((time.time() - start_time), 1)} сек. ---"
                )
                print("n_GTU_it", n_GTU_it)
                break
            if i == Max_iterations_minimum - 1:
                print(
                    "Достигнуто максимальное количество итераций минимального расхода в ПГУ",
                    i + 1,
                )
                print("n_GTU_it", n_GTU_it)

    gas_streams_zaryad = gas_streams
    water_streams_zaryad = water_streams
    syngas_streams_zaryad = syngas_streams
    electric_zaryad = electric
    heaters_zaryad = heaters
    accumulation_zaryad = accumulation

    #############ПРОЦЕСС РАЗРЯДКИ####################
    PKM_zaryad = False
    PKM_razryad = True

    Maxiterations_KU_TU = 20
    Maxiterations_cotel = 5
    Maxiterations_turbine = 20

    n_GTU = 1
    GTU_input.at["n", 1] = n_GTU

    arguments_all_it = [
        Maxiterations_KU_TU,
        Maxiterations_cotel,
        Maxiterations_turbine,
        gas_streams0,
        water_streams0,
        GTU_ISO,
        GTU_input,
        gas_streams,
        water_streams,
        heaters,
        electric,
        Gas_turbine,
        gas0,
        water,
        PKM_zaryad,
        PKM_razryad,
        syngas_streams,
        Calcmethod,
        Calctolerance,
        KPD_PN,
        KPD_KN,
        KPD_to,
        KPD_SP,
        steamVD_fraction_to_turbine,
        accumulation,
        time_ac,
        constr,
        time_jdat,
    ]
    New_Iter_pkm = 20
    pkm_pgu_tol = 10**-2
    CCGT = Calculate_CCGT_PKM_iter(arguments_all_it, New_Iter_pkm, pkm_pgu_tol)
    print(gas_streams)
    print(water_streams)
    print(syngas_streams)
    print(electric)
    print(heaters)
    print(accumulation)

    gas_streams_razryad = gas_streams
    water_streams_razryad = water_streams
    syngas_streams_razryad = syngas_streams
    electric_razryad = electric
    heaters_razryad = heaters
    accumulation_razryad = accumulation
    Potb2_turb = water_streams.at["DOOTB2", "P"]
    Potb2_teplof = water_streams.at["OTB2-SP2", "P"]
    Delta_P_Diafragma=Potb2_teplof-Potb2_turb
    
    time_all=time.time() -start_time_all
    print(f"Расчет окончен для температуры {air_temperature}:--- {round((time_all), 1)} сек. ---")

    result = {
        "T_air": round(air_temperature, 2),
        "n_GTU": round(GTU_input.at["n", 1], 5),
        "GTU": round(electric.at["GTU", "N"], 4),
        "GTU_KPD": round(electric.at["GTU", "KPD"], 4),
        "Turbine": round(electric.at["Turbine", "Ni"], 4),
        "KN": round(electric.at["KN", "Ni"], 4),
        "DK": round(electric.at["DK", "N"], 4),
        "PEN": round(electric.at["PEN", "Ni"], 4),
        "Turbine_Qt": round(heaters.at["SP2", "Qw"]+heaters.at["SP1", "Qw"]+heaters.at["OD", "Qw"], 4),
        # "PKM_Qt": round(accumulation.at["ASW", "Qw"]/(vremya*3600), 4),
        # "Delta_P_Diafragma":round(water_streams.at["DOOTB1", "P"]-water_streams.at["INCND", "P"],4),
        "INKOND": round(water_streams.at["INKOND", "G"], 4),
        # "Calculate_minimum": Сalculate_minimum,
        "Delta_P_Diafragma": round(Delta_P_Diafragma, 4),
        "T_accum": round(accumulation.at["PKM", "T"], 4),
        "time to calc":time_all,
        
        "gas_streams_zaryad": gas_streams_zaryad,
        "water_streams_zaryad": water_streams_zaryad,
        "syngas_streams_zaryad": syngas_streams_zaryad,
        "electric_zaryad": electric_zaryad,
        "heaters_zaryad": heaters_zaryad,
        "accumulation_zaryad": accumulation_zaryad,
        "gas_streams_razryad":   gas_streams_razryad,
        "water_streams_razryad": water_streams_razryad,
        "syngas_streams_razryad":   syngas_streams_razryad,
        "electric_razryad": electric_razryad,
        "heaters_razryad":  heaters_razryad,
        "accumulation_razryad":    accumulation_razryad,
    }
    print(result)
    return result
