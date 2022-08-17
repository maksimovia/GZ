


# Импорт библиотек
import numpy as n
import pandas as pd

def n_gtu_tatm(x):
    if x >= -2:
        Result = -0.006*x+1.091
    else:
        Result = 1.1
    n_gtu_tatm = Result
    return n_gtu_tatm


def n_gtu_pin(x):
    Result = -0.0015*x+1
    n_gtu_pin = Result
    return n_gtu_pin


def n_gtu_pout(x):
    Result = -0.00055*x+1
    n_gtu_pout = Result
    return n_gtu_pout

def gkt_gtu_tatm(x):  
    if x>=-2:
        Result=0.0000004251 * x ** 3 - 0.0000460133 * x ** 2 - 0.0024923115 * x + 1.0457022925
    else:
        Result =(0.000012573) * x * x + 0.0015678 * x + 1.0545
    gkt_gtu_tatm = Result
    return gkt_gtu_tatm

def gkt_gtu_n(x):  
    if x <= 0.53:
        Result= 0.019 * x + 0.685
    else:
        if x > 0.53:
            Result= 0.649 * x + 0.351
        else:
            if x >= 1:
                Result= 1
    gkt_gtu_n = Result
    return gkt_gtu_n

def gkt_gtu_pin(x):  
    Result= -0.001 * x + 1
    gkt_gtu_pin = Result
    return gkt_gtu_pin

def tkt_gtu_n(x):  
    if x<=0.53:
        Result=1.019 * x + 0.46
    else:
        Result =1
    tkt_gtu_n = Result
    return tkt_gtu_n

def tkt_gtu_tatm(x):  
    Result=0.000000057 * x ** 3 + 0.00000898 * x ** 2 + 0.000693897 * x + 0.987213678
    tkt_gtu_tatm = Result
    return tkt_gtu_tatm

def tkt_gtu_pin(x):  
    Result=0.000250404 * x + 1
    tkt_gtu_pin = Result
    return tkt_gtu_pin

def tkt_gtu_pout(x):  
    Result=2.67 * 10 ** -4 * x + 1
    tkt_gtu_pout = Result
    return tkt_gtu_pout

def eff_GTU_tatm(x):  
    if x>=0:
        Result=0.000000176 * x ** 3 - 0.000027204 * x ** 2 - 0.000959283 * x + 1.020148468
    else:
        Result = -0.000000254 * x ** 3 - 0.000018365 * x ** 2 - 0.000386071 * x + 1.020114209
    eff_GTU_tatm = Result
    return eff_GTU_tatm

def eff_gtu_pin(x):  
    Result=(20 - 0.012 * x) / 20
    eff_gtu_pin = Result
    return eff_gtu_pin

def eff_gtu_pout(x):  
    Result=-0.00055 * x + 1
    eff_gtu_pout = Result
    return eff_gtu_pout

def eff_gtu_n(x):  
    Result=-3.352525 * x ** 5 + 6.645291 * x ** 4 - 1.744334 * x ** 3 - 4.031808 * x ** 2 + 3.485841 * x - 0.00455
    eff_gtu_n = Result
    return eff_gtu_n

def n_pgu15(x):  
    Result=0.0007 * x * x + 0.19852 * x + 21.92782
    n_pgu15 = Result
    return n_pgu15

def n_pgu37(x):  
    Result=0.00022 * x * x + 0.36184 * x + 12.70505
    n_pgu37 = Result
    return n_pgu37

def n_pgu3(x):  
    Result=-0.00038 * x * x + 0.31242 * x + 17.7701
    n_pgu3 = Result
    return n_pgu3

def n_pgu10(x):  
    Result=-0.0005 * x * x + 0.31659 * x + 16.4482
    n_pgu10 = Result
    return n_pgu10

def n_pgu28(x):  
    Result=-0.00006 * x * x + 0.17351 * x + 19.75728
    n_pgu28 = Result
    return n_pgu28

def n_pgu40(x):  
    Result=-0.00032 * x * x + 0.24488 * x + 14.15046
    n_pgu40 = Result
    return n_pgu40

def g_pgu40(x):  
    Result=2.2047 * x + 116.5489
    g_pgu40 = Result
    return g_pgu40

def g_pgu28(x):  
    Result=2.2104 * x + 115.9821
    g_pgu28 = Result
    return g_pgu28

def g_pgu10(x):  
    Result=2.2201 * x + 115.6876
    g_pgu10 = Result
    return g_pgu10

def g_pgu3(x):  
    Result=2.220859 * x + 114.436408
    g_pgu3 = Result
    return g_pgu3

def g_pgu15(x):  
    Result=2.2705 * x + 105.7854
    g_pgu15 = Result
    return g_pgu15

def g_pgu37(x):  
    Result=2.3584 * x + 95.8445
    g_pgu37 = Result
    return g_pgu37

def q_pgu15(x):  
    if x>=111.9:
        Result=0.13747 * x + 96.1876
    else:
        Result = 111.6
    q_pgu15 = Result
    return q_pgu15

def q_pgu37(x):  
    if x>=96.9:
        Result=0.1571 * x + 96.44722
    else:
        Result = 111.6
    q_pgu37 = Result
    return q_pgu37

def q_pgu3(x):  
    Result=0.00123 * x * x + 0.4586 * x + 78.51066
    q_pgu3 = Result
    return q_pgu3

def q_pgu10(x):  
    Result=0.00127 * x * x + 0.42754 * x + 81.32958
    q_pgu10 = Result
    return q_pgu10

def q_pgu28(x):  
    Result=0.00112 * x * x + 0.43102 * x + 78.42656
    q_pgu28 = Result
    return q_pgu28

def q_pgu40(x):  
    Result=0.00077 * x * x + 0.51375 * x + 70.85735
    q_pgu40 = Result
    return q_pgu40

def t_ex37(x):  
    if x<= 110:
        Result= -0.0108876 * x * x + 2.3614201 * x + 64.983432
    else:
        Result = -0.027027 * x + 195.972973
    t_ex37 = Result
    return t_ex37

def t_ex15(x):  
    if x<= 123.6:
        Result= 0.1540541 * x + 172.7459459
    else:
        Result = 0.0359116 * x + 185.9779006
    t_ex15 = Result
    return t_ex15

def t_ex3(x):  
    if x<= 123.6:
        Result= 0.135468 * x + 171.7561576
    else:
        Result = 0.0827068 * x + 178.2774436
    t_ex3 = Result
    return t_ex3

def t_ex10(x):  
    if x<= 127.4:
        Result= 0.1375291 * x + 170.9787879
    else:
        Result = 0.0783848 * x + 178.5137767
    t_ex10 = Result
    return t_ex10

def t_ex28(x):  
    if x<= 127.4:
        Result= 0.1375291 * x + 169.9787879
    else:
        Result = 0.0783848 * x + 178.5137767
    t_ex28 = Result
    return t_ex28

def t_ex40(x):  
    if x<= 127.4:
        Result= 0.1445221 * x + 167.5878788
    else:
        Result = 0.0700935 * x + 177.0700935
    t_ex40 = Result
    return t_ex40



def gtu_raschet(ISO, n, T, Pin, Pout):
    N_out=ISO['N']*n_gtu_tatm(T)*n_gtu_pin(Pin)*n_gtu_pout(Pout)*n
    eff_out=ISO['eff']*eff_GTU_tatm(T)*eff_gtu_n(n)*eff_gtu_pin(Pin)*eff_gtu_pout(Pout)
    G_out=ISO['G']*gkt_gtu_tatm(T)*gkt_gtu_n(n)*gkt_gtu_pin(Pin)
    T_out=ISO['T']*tkt_gtu_tatm(T)*tkt_gtu_n(n)*tkt_gtu_pin(Pin)*tkt_gtu_pout(Pout) 
    N_dk=(2.121779534*N_out+5949.887022724)/1000
    return {'N':N_out,'eff':eff_out, 'G':G_out, 'T':T_out, 'Ndk':N_dk}


class gtu:
    def __init__(self, ISO_input, GTU_input, streamout):
        self.streamout = streamout
        self.GTU_input = GTU_input
        self.ISO = {
    "N": ISO_input.loc["Value", "N"],
    "eff": ISO_input.loc["Value", "eff"],
    "G": ISO_input.loc["Value", "G"],
    "T": ISO_input.loc["Value", "T"],
}
       

    def calc(self):
        T = self.GTU_input.at["tair", 1]
        n = self.GTU_input.at['n', 1]
        Pin = self.GTU_input.at['Pin', 1]
        Pout = self.GTU_input.at['Pout', 1]
        GTU_res = gtu_raschet(self.ISO, n, T, Pin , Pout)
        return GTU_res


