import mesa_reader as mr
import pickle

prefix = '/Users/ajermyn/Dropbox/Active_Projects/LowMassMagnetic/output/runs/'
DIR = prefix + 'magnetic_low_mass_Z_MW_profiles_time_2022_03_08_14_44_28_sha_f557' + '/runs/' # The directory where you unpacked the data

mods = '1.0 1.025 1.05 1.075 1.1 1.125 1.15 1.175 1.2 1.225 1.25 1.275 1.3 1.325 1.35 1.375 1.4 1.425 1.45 1.475 1.5 1.525 1.55 1.575 1.6 1.625 1.65 1.675 1.7 1.725 1.75 1.775 1.8 1.825 1.85 1.875 1.9 1.925 1.95 1.975'
mods = list(map(float,mods.split(' ')))

hs = list(mr.MesaData(DIR+str(j)+'/LOGS/history.data') for j in mods)
pf = list(mr.MesaData(DIR+str(j)+'/LOGS/profile_mid_MS.data') for j in mods)

pickle.dump([prefix,DIR,mods,hs,pf],open('parsed.data','wb'))

