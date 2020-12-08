from math import sqrt, log10, log2, pow, erfc, e
from scipy.special import erfcinv
from matplotlib import pyplot as plt

import argparse

#Parâmetros de referência: Link Budget de Uplink
uplink_pot_tx = 23 #dBm
uplink_Gtx = 0
uplink_tx_loss = 0
uplink_SNR_req = 0
uplink_sens_req_rx = -101.5 #dBm
uplink_Grx = 18 #dBi
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

    tri_area = radius*(radius*sqrt(3)/2)/2
    hex_area = 6*tri_area
    total_hex = area/hex_area

    return hex_area, total_hex

def link_budget(pot_tx, pot_rx, Ms, Gt, Pt, decrease):
    Lmax = pot_tx - pot_rx + Gt - Pt - Ms - decrease
    return Lmax

def max_radius(max_loss, freq):
    alfa = 4
    beta = 10.2
    gama = 2.36

    log_R_max = (max_loss - beta - 10*gama*log10(freq))/10*alfa
    R_max = pow(10, log_R_max)
    R_max = R_max/pow(10, 3)

    return R_max

def loss_model(dist, freq):
    alfa = 4
    beta = 10.2
    gama = 2.36

    loss = 10*alfa*log10(dist) + beta + 10*gama*log10(freq)

    return loss

def outdoor_radius(sigma, n, ro, gama, eta, bw, figura)

    Ms = 4*n/sigma - 3
    Q = (erfc(Ms/(sigma*sqrt(2)))/2
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
                temp = rate_total
                temp_sinr = sinr
            elif rate_total == desired_throughput:
                temp = rate_total
                temp_sinr = sinr
                break
            else:
                break

        mcs[modulation].append(temp_rate)
        mcs[modulation].append(temp_sinr)

        sens = temp_sinr + figura + 10*log10(180000) - 174 + D_in
        loss = link_budget(P_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss)
        return max_radius(loss, 2000)

def indoor_radius(sigma, n, ro, gama, eta, bw, figura)
    sigma = sqrt(64 + pow(sigma, 2))

    Ms = 4*n/sigma - 3
    Q = (erfc(Ms/(sigma*sqrt(2)))/2
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
                temp = rate_total
                temp_sinr = sinr
            elif rate_total == desired_throughput:
                temp = rate_total
                temp_sinr = sinr
                break
            else:
                break

        mcs[modulation].append(temp_rate)
        mcs[modulation].append(temp_sinr)

        sens = temp_sinr + figura + 10*log10(180000) - 174 + D_in
        loss = link_budget(P_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss)
        return max_radius(loss, 2000)

def incar_radius(sigma, n, ro, gama, eta, bw, figura)
    sigma = sqrt(36 + pow(sigma, 2))

    Ms = 4*n/sigma - 3
    Q = (erfc(Ms/(sigma*sqrt(2)))/2
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
                temp = rate_total
                temp_sinr = sinr
            elif rate_total == desired_throughput:
                temp = rate_total
                temp_sinr = sinr
                break
            else:
                break

        mcs[modulation].append(temp_rate)
        mcs[modulation].append(temp_sinr)

        sens = temp_sinr + figura + 10*log10(180000) - 174 + D_in
        loss = link_budget(P_tx,
                            sens,
                            Ms,
                            uplink_multipath + uplink_Gtx + uplink_Grx,
                            uplink_rx_loss + uplink_tx_loss)
        return max_radius(loss, 2000)

if __name__ == "__main__":

    inputs = {
        "bandwidth" = 0,
        "area" = 0,
        "throughput" = 0,
        "ro" = 0,
        "gama" = 0,
        "eta" = 0,
        "prob_outdoor" = 0,
        "prob_indoor" = 0,
        "prob_incar" = 0,
        "figura" = 0
    }

    #inputs
    inputs.bandwidth = input("Please select Bandwidth: ")
    inputs.area = input("Select desired area to be covered: ")
    inputs.throughput = input("Select desired throughput: ")
    inputs.ro = input("Select ro: ")
    inputs.gama = input("Select gama: ")
    inputs.eta = input("Select eta: ")
    inputs.prob_outdoor = input("Select probability of user being outdoors: ")
    inputs.prob_indoor = input("Select probability of user being indoors: ")
    inputs.prob_incar = input("Select probability of user being in a car: ")


    bw = inputs.bandwidth
    area = inputs.area
    desired_throughput = inputs.throughput
    ro = inputs.ro
    gama = inputs.gama
    eta = inputs.eta
    prob_outdoor = inputs.prob_outdoor
    prob_indoor = inputs.prob_indoor
    prob_incar = inputs.prob_incar


    n = 4
    sigma = 7.6
    prob_cobertura_celula = 0.9
    
    #substituir as entradas ->
    out_r = outdoor_radius(sigma, n, ro, gama, eta, bw, figura)
    indoor_r = indoor_radius(sigma, n, ro, gama, eta, bw, figura)
    incar_r = incar_radius(sigma, n, ro, gama, eta, bw, figura)

    radius = prob_outdoor*out_r + prob_indoor*indoor_r + prob_incar*incar_r

    mcs = start_mcs()
