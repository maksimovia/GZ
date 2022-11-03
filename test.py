import os
import time

import GTU
import KU_TU
import mat_properties as prop
import numpy as n
import pandas as pd
import SP
from scipy.optimize import root



def calculate_CCGT(Iterations_KU_TU, Iterations_cotel, Iterations_turbine,
                    gas_streams0,water_streams0,
                    GTU_ISO,GTU_input,
                    gas_streams,water_streams,
                    heaters,electric,
                    Gas_turbine,gas0,
                    water,PKM_zaryad,PKM_razryad,
                    syngas_streams,Calcmethod,
                    Calctolerance,KPD_PN,KPD_KN,KPD_to,KPD_SP,
                    steamVD_fraction_to_turbine,Teplo,accumulation,time_ac):
    
    gasmix = "Nitrogen*Oxygen*CO2*Water*Argon"
    RP = prop.init_REFPROP(r"C:\Program Files (x86)\REFPROP")

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

    #####################Максимов#####################
    from PKM import accum
    if PKM_zaryad:
        Accumulator = accum.zaryad(time_ac,accumulation,gas_streams,syngas_streams,water_streams,water_streams0,heaters,electric)
        steamVD_to_turbine = Accumulator['steamVD_to_turbine']
        #Мощность осн ГТУ????
    elif PKM_razryad:
        Accumulator = accum.razryad(time_ac,accumulation,gas_streams,syngas_streams,water_streams,water_streams0,heaters,electric)
        #......
        steamVD_to_turbine = Accumulator['steamVD_to_turbine']
    ##################################################

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
    
    # Инициализаця KU+TU, она здесь потому что нужно менять состав газа на входе в КУ

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
#     print(gas_streams)
    print(f"fin КУ и ТУ:--- {round((time.time() - start_time), 1)} сек. ---")
    return gas_streams