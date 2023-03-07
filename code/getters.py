import numpy as np
from numpy import log10
from functions import concat, rotation_law

def g(h,CZ,property):
	return getattr(h, concat(CZ, property))

def logB_shutoff_tau_3(CZ,h): # CZ is a dummy argument
	return log10(h.B_shutoff_down_to_tau_3)

def logB_shutoff_tau_10(CZ,h): # CZ is a dummy argument
	return log10(h.B_shutoff_down_to_tau_10)

def logB_shutoff_tau_100(CZ,h): # CZ is a dummy argument
	return log10(h.B_shutoff_down_to_tau_100)

def logB_shutoff(CZ, h):
	return log10(g(h,CZ,'B_shutoff'))

def logPr(CZ, h):
	return log10(g(h,CZ,'viscosity'))-log10(g(h,CZ,'thermal_diffusivity'))

def logPr_local(CZ, h):
	return log10(g(h,CZ,'prandtl'))

def logPm(CZ, h):
	return log10(g(h,CZ,'viscosity'))-log10(g(h,CZ,'magnetic_diffusivity'))

def logBeta(CZ, h):
	return log10(g(h,CZ,'Prad_div_P'))

def logA(CZ, h):
	dR = g(h,CZ,'cz_top_r') - g(h,CZ,'cz_bottom_r')
	return log10(g(h,CZ,'cz_top_r')) - log10(dR)

def logD(CZ, h):
	return log10(g(h,CZ,'density_ratio'))

def logLcdivL(CZ, h):
	return log10(g(h,CZ,'L_conv_div_L'))

def logTurnover(CZ, h):
	return log10(g(h,CZ,'turnover_time') / (24*3600))

def logRa(CZ, h):
	dGrad = g(h,CZ,'gradr_sub_grada')
	dR = g(h,CZ,'cz_top_r') - g(h,CZ,'cz_bottom_r')
	grav = g(h,CZ,'gravity')
	nu = g(h,CZ,'viscosity')
	alpha = g(h,CZ,'thermal_diffusivity')
	hp = g(h,CZ,'pressure_scale_height')
	radiation_term = 1 - g(h,CZ,'Prad_div_P') # Gives Pgas/P = beta_gas
	radiation_term = (4 - 3 * radiation_term) / radiation_term
	return log10(dGrad) + 4*log10(dR) - log10(hp) + log10(grav) - log10(nu) - log10(alpha) + log10(radiation_term)

def logTauCZ(CZ, h):
	return log10(g(h,CZ,'cz_dtau'))

def logTauSurf(CZ, h):
	return log10(g(h,CZ,'cz_top_tau'))

def logGammaEdd(CZ, h):
	return log10(g(h,CZ,'L_div_Ledd'))

def logRe(CZ, h):
	dR = g(h,CZ,'cz_top_r') - g(h,CZ,'cz_bottom_r')
	return log10(g(h,CZ,'conv_vel')) + log10(dR) - log10(g(h,CZ,'viscosity'))

def logPe(CZ, h):
	dR = g(h,CZ,'cz_top_r') - g(h,CZ,'cz_bottom_r')
	return log10(g(h,CZ,'conv_vel')) + log10(dR) - log10(g(h,CZ,'thermal_diffusivity'))

def logMad(CZ, h):
	return log10(g(h,CZ,'conv_vel')) - log10(g(h,CZ,'sound_speed_adiabatic'))

def logMiso(CZ, h):
	return log10(g(h,CZ,'conv_vel')) - log10(g(h,CZ,'sound_speed_isothermal'))

def logNu(CZ, h):
	return log10(g(h,CZ,'Nusselt'))

def logStop(CZ, h):
	return log10(g(h,CZ,'stiffness_top'))

def logSbot(CZ, h):
	return log10(g(h,CZ,'stiffness_bottom'))

def logEk(CZ, h):
	dR = g(h,CZ,'cz_top_r') - g(h,CZ,'cz_bottom_r')
	return log10(g(h,CZ,'viscosity') / (2 * rotation_law(h) * dR**2))

def logRo(CZ, h):
	dR = g(h,CZ,'cz_top_r') - g(h,CZ,'cz_bottom_r')
	return log10(g(h,CZ,'conv_vel') / (2 * rotation_law(h) * dR))