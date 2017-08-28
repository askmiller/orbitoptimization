from scipy.special import ellipk
from scipy.special import ellipe
import math

G = 6.673e-11 # m^3 /kg /s^2, the gravitational constant

####################################################################################
#
# Each member, p, of the population must be a list of pairs (mass, radius)
#
####################################################################################
#
# If p has length n, fitness_vector returns a list of length n + 2
# whose first n entries are stability value checks which should be non-negative,
# whose last two entries are the angular momentum L, and the total energy.
#
####################################################################################

def fitness_vector(p):
		assert p[0][1] == 0.0
		vector = []
		L_sum = 0.0
		GE_KE_sum = 0.0
		for mj,rj in p:
			L_inner_sum = 0.0
			for mi,ri in p:
				if ri == rj:
					continue
				elliptic_m_value = 4*ri * rj/(ri + rj)/(ri + rj)
				K = ellipk(elliptic_m_value)
				E = ellipe(elliptic_m_value)
				L_inner_summand =  K/(rj + ri) + E/(rj - ri)
				L_inner_sum = L_inner_sum + L_inner_summand
				GE_KE_summand =   - K/(rj + ri) + E/(rj - ri)
				GE_KE_sum = GE_KE_sum + GE_KE_summand 

			vector.append(L_inner_sum)
			if L_inner_sum >= 0.0:
				L_sum = L_sum + L_inner_sum**0.5*mj*rj
					#Tolerate invalid values
					#Rely on fitness penalty to get rid of them


		L = (G/math.pi)**0.5*L_sum
		vector.append(L)
		GE_KE = (G/math.pi/2.0)* GE_KE_sum
		vector.append(GE_KE)
		return vector
