# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using HDF5

const c_to_mps = 299792458
const mps_to_c = 3.335641e-09

const lsol_to_erg= 3.828e33
const erg_to_lsol = 2.61233e-34

const pc_to_m = 3.08568e+16
const m_to_pc = 3.240777e-17

const kpc_to_m = 3.08568e+19
const m_to_kpc = 3.240777e-20

const mpc_to_m = 3.08568e+22
const m_to_mpc = 3.240777e-23

const mpc_to_cm = 3.08568e+24
const cm_to_mpc = 3.240777e-25

const msol_to_kg = 1.989e30
const kg_to_msol = 5.027652e-31

const lsol_to_w = 3.828e26
const w_to_lsol = 2.61233e-27

const jansky_to_cgs = 1e-23
const cgs_to_jansky = 1e23

const jansky_to_si = 1e-26
const si_to_jansky = 1e26

const cgs_to_si = 1e-3
const si_to_cgs = 1e3

const BC03lr = bc_data()
