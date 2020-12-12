from math import sqrt, log10, log2, pow, erfc, e
from scipy.special import erfcinv
from matplotlib import pyplot as plt

import argparse

#Parâmetros de referência: Link Budget de Uplink
uplink_pot_tx = 26 #dBm
uplink_Gtx = 1.2
uplink_tx_loss = 0
uplink_SNR_req = 0
uplink_sens_req_rx = -101.5 #dBm
uplink_Grx = 19 #dBi
uplink_rx_loss = 3 #dB
uplink_multipath = 3 #dB
uplink_shadow_margin = 4 #dB

def start_mcs():
    # 0 -> bit rate
    # 1 -> code rate

    mcs = {
        "QPSK 1/3":[2, 0.333333, 2342010, 14.0051, -0.577897],
        "16-QAM 1/2":[4, 0.5, 47613.1, 0.0926275, -0.29583],
        "64-QAM 3/4":[6, 0.75, 26405.8, 0.0220186, -0.24491]}

    return mcs

def count_hex(area, radius):

    hex_area_m = 1.5*sqrt(3)*pow(radius, 2)
    hex_area_km = hex_area_m/1000000
    total_hex = area//hex_area_km

    return hex_area_km, total_hex

def link_budget(pot_tx, pot_rx, Ms, Gt, Pt, decrease):
    Lmax = pot_tx - pot_rx + Gt - Pt - Ms - decrease
    return Lmax

def max_radius(max_loss, freq, modulation):
    alfa = 4
    beta = 10.2
    gama = 2.36

    log_R_max = (max_loss - beta - 10*gama*log10(freq))/(10*alfa)
    R_max = pow(10, log_R_max)
    #print(modulation, "log:", log_R_max, "normal:", R_max)

    return R_max

def loss_model(dist, freq):
    alfa = 4
    beta = 10.2
    gama = 2.36

    loss = 10*alfa*log10(dist) + beta + 10*gama*log10(freq)

    return loss

def outdoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs):
    #prob_cobertura_borda = float()
    Ms = 4*n/sigma - 3
    Q = 1 - (erfc(Ms/(sigma*sqrt(2))))/2
    #prob_cobertura_borda = 1 - Q
     
    sigma_sir = sqrt(2*sigma*(1-ro))
    Q_inv = sqrt(2)*erfcinv(2*Q)
    M_in = -1*Q_inv*sigma_sir

    m_in = pow(10, M_in/10)
    D_in = 10*log10(m_in*gama*eta)

    setor_gain = 4.77
    rate_rb = float()
    rate_total = float()
    temp_rate = float()
    temp_sinr = float()
    
    for modulation in mcs:
        for sinr in reversed(range(0, 34)):
            rate_rb = mcs[modulation][2]/(mcs[modulation][3] + pow(e, mcs[modulation][4]*sinr))
            rate_total = rate_rb*(5*bw)

            if rate_total > desired_throughput:
                temp_rate = rate_total
                temp_sinr = sinr
            elif rate_total == desired_throughput:
                temp_rate = rate_total
                temp_sinr = sinr
                break
            else:
                break

        mcs[modulation].append(temp_rate)
        mcs[modulation].append(temp_sinr)

        sir_mean = m_in*gama*temp_sinr + setor_gain
        #print("SIR mínima:", temp_sinr, "SIR média:", sir_mean)
        mcs[modulation].append(sir_mean)

        sens = temp_sinr + figura + 10*log10(180000) - 174 + D_in
        loss = link_budget(uplink_pot_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss,
                            0)
        mcs[modulation].append(max_radius(loss, 2.6, modulation))

    return mcs

def indoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs):
    sigma = sqrt(64 + pow(sigma, 2))

    Ms = 4*n/sigma - 3
    Q = (erfc(Ms/(sigma*sqrt(2))))/2
    prob_cobertura_borda = 1 - Q

    sigma_sir = sqrt(2*sigma*(1-ro))
    Q_inv = sqrt(2)*erfcinv(2*prob_cobertura_borda)
    M_in = -1*Q_inv*sigma_sir

    m_in = pow(10, M_in/10)
    D_in = 10*log10(m_in*gama*eta)

    rate_rb = float()
    rate_total = float()
    temp_rate = float()
    temp_sinr = float()

    for modulation in mcs:
        for sinr in reversed(range(0, 34)):
            rate_rb = mcs[modulation][2]/(mcs[modulation][3] + pow(e, mcs[modulation][4]*sinr))
            rate_total = rate_rb*(5*bw)

            if rate_total > desired_throughput:
                temp_rate = rate_total
                temp_sinr = sinr
            elif rate_total == desired_throughput:
                temp = rate_total
                temp_sinr = sinr
                break
            else:
                break

        #mcs[modulation].append(temp_rate)
        #mcs[modulation].append(temp_sinr)

        sens = temp_sinr + figura + 10*log10(180000) - 174 + D_in
        loss = link_budget(uplink_pot_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss,
                            20)
        mcs[modulation].append(max_radius(loss, 2.6, modulation))
    
    return mcs

def incar_radius(sigma, n, ro, gama, eta, bw, figura, mcs):
    sigma = sqrt(9 + pow(sigma, 2))

    Ms = 4*n/sigma - 3
    Q = (erfc(Ms/(sigma*sqrt(2))))/2
    prob_cobertura_borda = 1 - Q

    sigma_sir = sqrt(2*sigma*(1-ro))
    Q_inv = sqrt(2)*erfcinv(2*prob_cobertura_borda)
    M_in = -1*Q_inv*sigma_sir

    m_in = pow(10, M_in/10)
    D_in = 10*log10(m_in*gama*eta)

    rate_rb = float()
    rate_total = float()
    temp_rate = float()
    temp_sinr = float()

    for modulation in mcs:
        for sinr in reversed(range(0, 34)):
            rate_rb = mcs[modulation][2]/(mcs[modulation][3] + pow(e, mcs[modulation][4]*sinr))
            rate_total = rate_rb*(5*bw)

            if rate_total > desired_throughput:
                temp_rate = rate_total
                temp_sinr = sinr
            elif rate_total == desired_throughput:
                temp = rate_total
                temp_sinr = sinr
                break
            else:
                break

        #mcs[modulation].append(temp_rate)
        #mcs[modulation].append(temp_sinr)

        sens = temp_sinr + figura + 10*log10(180000) - 174 + D_in
        loss = link_budget(uplink_pot_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss,
                            10)
        mcs[modulation].append(max_radius(loss, 2.6, modulation))
    
    return mcs

if __name__ == "__main__":

    inputs = {
        "bandwidth": 0,
        "area": 0,
        "throughput": 0,
        "ro": 0,
        "gama": 0,
        "eta": 0,
        "prob_outdoor": 0,
        "prob_indoor": 0,
        "prob_incar": 0,
        "figura": 0
    }

    #inputs
    inputs["bandwidth"] = float(input("Please select Bandwidth: "))
    inputs["area"] = float(input("Select desired area to be covered: "))
    inputs["throughput"] = float(input("Select desired throughput: "))
    inputs["ro"] = float(input("Select ro: "))
    inputs["gama"] = float(input("Select gama: "))
    inputs["eta"] = float(input("Select eta: "))
    inputs["prob_outdoor"] = float(input("Select probability of user being outdoors: "))
    inputs["prob_indoor"] = float(input("Select probability of user being indoors: "))
    inputs["prob_incar"] = float(input("Select probability of user being in a car: "))
    inputs["figura"] = float(input("Select noise figure: "))

    bw = inputs["bandwidth"]
    area = inputs["area"]
    desired_throughput = inputs["throughput"]
    ro = inputs["ro"]
    gama = inputs["gama"]
    eta = inputs["eta"]
    prob_outdoor = inputs["prob_outdoor"]
    prob_indoor = inputs["prob_indoor"]
    prob_incar = inputs["prob_incar"]
    figura = inputs["figura"]
    mcs = start_mcs()
    n = 4
    sigma = 7.6
    prob_cobertura_celula = 0.9
    
    #substituir as entradas ->
    mcs = outdoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs)
    mcs = indoor_radius(sigma, n, ro, gama, eta, bw, figura, mcs)
    mcs = incar_radius(sigma, n, ro, gama, eta, bw, figura, mcs)

    rate_list = list()

    for modulation in mcs:
        rate_list.append(mcs[modulation][5])

        radius = prob_outdoor*mcs[modulation][8] + prob_indoor*mcs[modulation][9] + prob_incar*mcs[modulation][10]
        cell_area, cell_quant = count_hex(area, radius*0.95)
        mcs[modulation].append(cell_area)
        mcs[modulation].append(cell_quant)
        print(cell_area, cell_quant)
